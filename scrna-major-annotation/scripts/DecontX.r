#!/usr/bin/env Rscript
#
# 功能：对 seurat-qc-metrics.R 产生的项目级 Seurat RDS 按 orig.ident 拆分，以单样本/10X channel 为建模单位并行运行 DecontX。
# 输入：
#   1. --input：seurat-qc-metrics.R 实际生成的项目级 RDS 路径（当前命名规则为 <project>.qc-metrics.rds，不是字面固定文件名）；必须保留 RNA 原始非负整数 UMI counts，orig.ident 必须非空。
#   2. --project 与 --output-dir：项目名和输出根目录；必须与 workflow_state.json 中的 qc_metrics 事件一致。
#   3. 可选 --background：仅单样本时使用的同 channel 10X raw 目录/H5 或 counts-bearing RDS。
#   4. 可选 --background-sheet：CSV/TSV，前两列必须依次为 sample_id,background_path，每个 orig.ident 恰好一行。
#   5. 可选 --background-is-empty-only：仅当 background 已预先严格去除所有 called-cell barcodes 时使用。
# 输出：
#   1. <output-dir>/count-versions/original/<project>.decontx-original.rds：RNA 为原始 counts，并新增 decontX_contamination/decontX_clusters metadata。
#   2. <output-dir>/count-versions/corrected/<project>.decontx-corrected.rds：RNA 为 DecontX 矫正后整数 counts，重算 count 相关 QC/标准化/细胞周期，并保留同一污染与原始-count scDblFinder metadata。
#   3. <output-dir>/<project>.decontx-cell-metadata.csv 与 <project>.decontx-sample-summary.csv：逐细胞和逐样本污染摘要。
#   4. <output-dir>/<project>.decontx-marker-plots.<input_md5>/：每样本污染校正前/后 marker detection PDF 和高分辨率 PNG。
#   5. <output-dir>/<project>.count-version-manifest.csv、<project>.count-version-review.csv：双分支清单和第一轮比较后的人工选择门。
#   6. 根目录及两个 count-version 分支各自的 workflow_state.json：分支独立运行第一轮，避免结果文件互相覆盖。
# 示例命令：
#   Rscript scripts/DecontX.r --input output/prostate.qc-metrics.rds --project prostate --output-dir output
#   Rscript scripts/DecontX.r --input output/prostate.qc-metrics.rds --project prostate --output-dir output --background-sheet backgrounds.csv
# 注意事项：
#   1. DecontX 内置 broad clustering 作为 z；官方教程指出 cluster 选择会直接影响污染估计，因此结果必须结合 marker 校正前/后图复核。
#   2. background 必须与样本来自同一 channel，基因集必须完全一致；默认要求 raw background 完整包含本样本 called-cell barcodes，并按精确 barcode 排除。
#   3. 两个对象第一轮使用完全相同的细胞；scDblFinder 只在原始 observed counts 上计算并复制到 corrected 对象，不在矫正 counts 上重跑。
#   4. 两个分支均完成第一轮聚类、QC 和 provisional annotation 后，人工在 count-version-review.csv 选择一个版本；未选择的分支不得进入过滤及后续轮次。

script_arg <- grep("^--file=", commandArgs(), value = TRUE)
script_dir <- if (length(script_arg)) dirname(normalizePath(sub("^--file=", "", script_arg[[1]]))) else getwd()
source(file.path(script_dir, "common.R"))
source(file.path(script_dir, "workflow-state.R"))

as_dgC_matrix <- function(x) {
  x <- methods::as(x, "generalMatrix")
  x <- methods::as(x, "CsparseMatrix")
  methods::as(x, "dMatrix")
}

validate_count_matrix <- function(x, label) {
  if (is.null(dim(x)) || length(dim(x)) != 2L || !nrow(x) || !ncol(x)) {
    stop(label, " must be a non-empty gene-by-cell count matrix.")
  }
  if (is.null(rownames(x)) || is.null(colnames(x)) || anyNA(rownames(x)) || anyNA(colnames(x)) ||
      any(!nzchar(rownames(x))) || any(!nzchar(colnames(x))) || anyDuplicated(rownames(x)) || anyDuplicated(colnames(x))) {
    stop(label, " must have non-empty, unique gene and cell names.")
  }
  x <- as_dgC_matrix(x)
  values <- x@x
  if (any(!is.finite(values)) || any(values < 0) || any(abs(values - round(values)) > 1e-8)) {
    stop(label, " must contain finite, non-negative integer UMI values.")
  }
  Matrix::drop0(x)
}

resolve_background_path <- function(path) {
  path <- normalizePath(path, mustWork = TRUE)
  if (dir.exists(path) && basename(path) == "outs") {
    candidates <- c(file.path(path, "raw_feature_bc_matrix.h5"), file.path(path, "raw_feature_bc_matrix"))
    found <- candidates[file.exists(candidates)]
    if (!length(found)) stop("Cell Ranger outs directory lacks raw_feature_bc_matrix.h5 or raw_feature_bc_matrix/.")
    path <- found[[1]]
  }
  normalizePath(path, mustWork = TRUE)
}

background_fingerprint <- function(path) {
  path <- resolve_background_path(path)
  files <- if (!dir.exists(path)) {
    path
  } else {
    choose_file <- function(candidates, label) {
      found <- file.path(path, candidates)
      found <- found[file.exists(found)]
      if (!length(found)) stop("Background directory lacks ", label, ": ", path)
      found[[1]]
    }
    c(
      choose_file(c("matrix.mtx.gz", "matrix.mtx"), "matrix.mtx[.gz]"),
      choose_file(c("features.tsv.gz", "features.tsv", "genes.tsv.gz", "genes.tsv"), "features/genes TSV"),
      choose_file(c("barcodes.tsv.gz", "barcodes.tsv"), "barcodes.tsv[.gz]")
    )
  }
  files <- normalizePath(files, mustWork = TRUE)
  hashes <- unname(as.character(tools::md5sum(files)))
  list(
    resolved_path = path,
    files = files,
    md5 = paste(paste0(basename(files), "=", hashes), collapse = ";")
  )
}

read_gene_expression <- function(path) {
  path <- resolve_background_path(path)
  value <- if (dir.exists(path)) {
    Seurat::Read10X(path)
  } else if (tolower(tools::file_ext(path)) %in% c("h5", "hdf5")) {
    Seurat::Read10X_h5(path)
  } else {
    readRDS(path)
  }
  if (inherits(value, "Seurat")) {
    value <- join_rna_layers(value)
    value <- rna_layer(value, "counts")
  } else if (inherits(value, "SingleCellExperiment")) {
    if (!"counts" %in% SummarizedExperiment::assayNames(value)) stop("Background SCE lacks a counts assay.")
    value <- SummarizedExperiment::assay(value, "counts")
  } else if (is.list(value)) {
    if (!"Gene Expression" %in% names(value)) stop("Background contains multiple modalities but no 'Gene Expression' matrix.")
    value <- value[["Gene Expression"]]
  }
  validate_count_matrix(value, paste0("Background counts from ", path))
}

read_background_sheet <- function(path, sample_ids) {
  path <- normalizePath(path, mustWork = TRUE)
  ext <- tolower(tools::file_ext(path))
  sheet <- if (ext == "csv") {
    utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  } else if (ext %in% c("tsv", "txt")) {
    utils::read.delim(path, stringsAsFactors = FALSE, check.names = FALSE)
  } else {
    stop("--background-sheet must be .csv, .tsv, or .txt.")
  }
  if (ncol(sheet) < 2L || !identical(colnames(sheet)[1:2], c("sample_id", "background_path"))) {
    stop("The first two background-sheet columns must be ordered as sample_id, background_path.")
  }
  sheet$sample_id <- trimws(as.character(sheet$sample_id))
  sheet$background_path <- trimws(as.character(sheet$background_path))
  if (!nrow(sheet) || anyNA(sheet$sample_id) || anyNA(sheet$background_path) ||
      any(!nzchar(sheet$sample_id)) || any(!nzchar(sheet$background_path)) || anyDuplicated(sheet$sample_id)) {
    stop("background-sheet sample_id and background_path must be non-empty; sample_id must be unique.")
  }
  if (!setequal(sheet$sample_id, sample_ids)) {
    stop("background-sheet sample_id values must exactly match orig.ident values.")
  }
  relative <- !grepl("^/|^[A-Za-z]:[/\\\\]", sheet$background_path)
  sheet$background_path[relative] <- file.path(dirname(path), sheet$background_path[relative])
  paths <- setNames(sheet$background_path, sheet$sample_id)
  paths[sample_ids]
}

called_barcode_candidates <- function(sample_id, cells) {
  prefix <- paste0(sample_id, "_")
  stripped <- ifelse(startsWith(cells, prefix), substring(cells, nchar(prefix) + 1L), cells)
  unique(c(cells, stripped))
}

prepare_background <- function(path, sample_id, sample_counts, empty_only = FALSE) {
  fingerprint <- background_fingerprint(path)
  background <- read_gene_expression(fingerprint$resolved_path)
  if (!setequal(rownames(background), rownames(sample_counts))) {
    stop("Background genes do not exactly match filtered-cell genes for sample ", sample_id, ".")
  }
  background <- background[rownames(sample_counts), , drop = FALSE]
  raw_n <- ncol(background)
  called_ids <- called_barcode_candidates(sample_id, colnames(sample_counts))
  called <- colnames(background) %in% called_ids
  if (empty_only && any(called)) {
    stop("--background-is-empty-only was set, but called-cell barcodes were found for sample ", sample_id, ".")
  }
  if (!empty_only && sum(called) != ncol(sample_counts)) {
    stop(
      "Raw background must contain every filtered-cell barcode for sample ", sample_id,
      "; matched ", sum(called), " of ", ncol(sample_counts),
      ". Confirm that background comes from the same channel, or use --background-is-empty-only for a prefiltered empty-droplet matrix."
    )
  }
  background <- background[, !called, drop = FALSE]
  if (!ncol(background)) stop("No empty droplets remain after removing called cells for sample ", sample_id, ".")
  list(
    counts = background, raw_droplets = raw_n, called_removed = sum(called),
    empty_droplets = ncol(background), resolved_path = fingerprint$resolved_path,
    fingerprint = fingerprint$md5
  )
}

validate_decontx_result <- function(corrected, contamination, expected_counts, sample_id) {
  corrected <- as_dgC_matrix(corrected)
  if (!identical(dim(corrected), dim(expected_counts)) ||
      !identical(rownames(corrected), rownames(expected_counts)) ||
      !identical(colnames(corrected), colnames(expected_counts))) {
    stop("DecontX corrected-count dimensions or dimnames changed for sample ", sample_id, ".")
  }
  values <- corrected@x
  if (any(!is.finite(values)) || any(values < 0) || any(abs(values - round(values)) > 1e-8)) {
    stop("DecontX corrected counts are not finite, non-negative integer values for sample ", sample_id, ".")
  }
  if (length(contamination) != ncol(expected_counts) || any(!is.finite(contamination)) ||
      any(contamination < 0 | contamination > 1)) {
    stop("Invalid DecontX contamination estimates for sample ", sample_id, ".")
  }
  Matrix::drop0(corrected)
}

marker_sets <- function(features) {
  requested <- list(
    T_cells = c("CD3D", "CD3E", "TRAC"),
    B_cells = c("MS4A1", "CD79A", "CD79B"),
    Monocytes = c("LYZ", "S100A8", "S100A9"),
    NK_cells = c("GNLY", "NKG7"),
    epithelial = c("EPCAM", "KRT8", "KRT18"),
    fibroblast = c("DCN", "LUM", "COL1A1")
  )
  used <- lapply(requested, intersect, y = features)
  used <- used[lengths(used) > 0L]
  if (!length(used)) stop("No bundled T/B/Monocyte/NK/Epithelial/Fibroblast marker genes were found; marker QC plots cannot be generated.")
  list(requested = requested, used = used)
}

safe_file_tag <- function(x) {
  tag <- gsub("[^A-Za-z0-9._-]+", "_", x)
  ifelse(nzchar(tag), tag, "sample")
}

save_marker_plots <- function(results, original_counts, corrected_counts, markers, output_dir) {
  ensure_dir(output_dir)
  tags <- make.unique(vapply(names(results), safe_file_tag, character(1)), sep = "_")
  files <- vector("list", length(results))
  for (i in seq_along(results)) {
    result <- results[[i]]
    cells <- result$cells
    sce <- SingleCellExperiment::SingleCellExperiment(
      assays = list(
        counts = original_counts[, cells, drop = FALSE],
        decontXcounts = corrected_counts[, cells, drop = FALSE]
      ),
      colData = S4Vectors::DataFrame(
        decontX_clusters = result$clusters,
        row.names = cells
      )
    )
    plot <- decontX::plotDecontXMarkerPercentage(
      sce, markers = markers, assayName = c("counts", "decontXcounts")
    ) + ggplot2::ggtitle(paste0(result$sample_id, ": marker detection before/after DecontX"))
    prefix <- file.path(output_dir, paste0(tags[[i]], ".marker-before-after"))
    save_plot_pair(plot, prefix, width = 12, height = 7, dpi = 300)
    files[[i]] <- c(paste0(prefix, ".pdf"), paste0(prefix, ".png"))
  }
  unname(unlist(files, use.names = FALSE))
}

atomic_write_lines <- function(x, path) {
  path <- absolute_path(path, must_work = FALSE)
  tmp <- tempfile(pattern = paste0(".", basename(path), "."), tmpdir = dirname(path), fileext = ".tmp")
  on.exit(if (file.exists(tmp)) unlink(tmp), add = TRUE)
  writeLines(x, tmp, useBytes = TRUE)
  replace_file(tmp, path)
}

build_corrected_object <- function(original, corrected_counts, project) {
  derived <- c(
    "nCount_RNA", "nFeature_RNA", "percent_mito", "percent_ribo", "percent_hb",
    "S.Score", "G2M.Score", "Phase"
  )
  metadata <- original[[]]
  metadata <- metadata[, setdiff(colnames(metadata), derived), drop = FALSE]
  corrected <- Seurat::CreateSeuratObject(
    counts = corrected_counts, assay = "RNA", project = project,
    meta.data = metadata, min.cells = 0, min.features = 0
  )
  if (!identical(colnames(corrected), colnames(original)) || !identical(rownames(corrected), rownames(original))) {
    stop("Corrected Seurat object changed the input cell or gene order.")
  }
  stored_counts <- rna_layer(corrected, "counts")
  if (!identical(as_dgC_matrix(stored_counts), as_dgC_matrix(corrected_counts))) {
    stop("Corrected Seurat RNA counts do not exactly match the validated DecontX count matrix.")
  }
  corrected$nCount_RNA <- as.numeric(Matrix::colSums(corrected_counts))
  corrected$nFeature_RNA <- as.numeric(Matrix::colSums(corrected_counts > 0))
  assert_positive_qc_counts(corrected)
  corrected <- Seurat::PercentageFeatureSet(corrected, pattern = "^MT-", col.name = "percent_mito")
  corrected <- Seurat::PercentageFeatureSet(corrected, pattern = "^RP[SL]", col.name = "percent_ribo")
  corrected <- Seurat::PercentageFeatureSet(corrected, pattern = "^HB[ABDEGMQZ]", col.name = "percent_hb")
  corrected <- Seurat::NormalizeData(
    corrected, normalization.method = "LogNormalize", scale.factor = 10000, verbose = FALSE
  )
  cc <- Seurat::cc.genes.updated.2019
  s_genes <- intersect(cc$s.genes, rownames(corrected))
  g2m_genes <- intersect(cc$g2m.genes, rownames(corrected))
  if (length(s_genes) < 5L || length(g2m_genes) < 5L) {
    stop("Too few human cell-cycle genes remain in corrected counts for reliable scoring.")
  }
  corrected <- Seurat::CellCycleScoring(
    corrected, s.features = s_genes, g2m.features = g2m_genes, set.ident = FALSE
  )
  corrected@misc <- original@misc
  corrected
}

count_version_review_template <- function(project, input_file, input_md5, shared_cell_set_md5,
                                          original_dir, corrected_dir, original_file, corrected_file) {
  data.frame(
    project = project,
    parent_qc_rds = input_file,
    parent_qc_md5 = input_md5,
    shared_cell_set_md5 = shared_cell_set_md5,
    original_branch_dir = original_dir,
    corrected_branch_dir = corrected_dir,
    original_rds = original_file,
    original_rds_md5 = file_md5(original_file),
    corrected_rds = corrected_file,
    corrected_rds_md5 = file_md5(corrected_file),
    agent_version = "", agent_reason = "", agent_evidence = "",
    human_version = "", human_reason = "", approved = "0",
    reviewer = "", reviewed_at = "", stringsAsFactors = FALSE
  )
}

opts <- parse_cli_args()
input_file <- normalizePath(required_opt(opts, "input"), mustWork = TRUE)
project <- required_opt(opts, "project")
output_dir <- ensure_dir(required_opt(opts, "output_dir"))
state_file <- file.path(output_dir, "workflow_state.json")
input_md5 <- file_md5(input_file)
count_versions_dir <- ensure_dir(file.path(output_dir, "count-versions"))
original_dir <- ensure_dir(file.path(count_versions_dir, "original"))
corrected_dir <- ensure_dir(file.path(count_versions_dir, "corrected"))
original_file <- file.path(original_dir, paste0(project, ".decontx-original.rds"))
corrected_file <- file.path(corrected_dir, paste0(project, ".decontx-corrected.rds"))
original_state_file <- file.path(original_dir, "workflow_state.json")
corrected_state_file <- file.path(corrected_dir, "workflow_state.json")
cell_metadata_file <- file.path(output_dir, paste0(project, ".decontx-cell-metadata.csv"))
sample_summary_file <- file.path(output_dir, paste0(project, ".decontx-sample-summary.csv"))
marker_dir <- file.path(output_dir, paste0(project, ".decontx-marker-plots.", substr(input_md5, 1L, 12L)))
session_file <- file.path(output_dir, paste0(project, ".decontx-sessionInfo.txt"))
manifest_file <- file.path(output_dir, paste0(project, ".count-version-manifest.csv"))
count_version_review_file <- file.path(output_dir, paste0(project, ".count-version-review.csv"))
assert_output_not_input(input_file, original_file)
assert_output_not_input(input_file, corrected_file)

seed <- if (is.null(opts$seed)) WORKFLOW_SEED else parse_integer(opts$seed, "seed", 1L, .Machine$integer.max - 10000L)
background <- opts$background %||% NULL
background_sheet <- opts$background_sheet %||% NULL
background_is_empty_only <- isTRUE(opts$background_is_empty_only)
if (!is.null(background) && !is.null(background_sheet)) stop("Use only one of --background or --background-sheet.")
if (background_is_empty_only && is.null(background) && is.null(background_sheet)) {
  stop("--background-is-empty-only requires --background or --background-sheet.")
}

packages <- c(
  "Seurat", "SeuratObject", "decontX", "SingleCellExperiment", "SummarizedExperiment",
  "S4Vectors", "Matrix", "future", "future.apply", "parallelly", "ggplot2", "jsonlite"
)

tryCatch({
  require_packages(packages, minimum = list(Seurat = "5.0.0", SeuratObject = "5.0.0"))
  assert_decontx_input(state_file, project, input_file)
  set.seed(seed)
  object <- readRDS(input_file)
  assert_seurat_input(object)
  object <- join_rna_layers(object)
  counts <- validate_count_matrix(rna_layer(object, "counts"), "RNA counts")
  sample_vector <- as.character(object$orig.ident)
  sample_ids <- sort(unique(sample_vector))
  if (!is.null(background) && length(sample_ids) != 1L) {
    stop("--background is only valid for one orig.ident; use --background-sheet for multiple samples.")
  }
  background_paths <- setNames(rep(NA_character_, length(sample_ids)), sample_ids)
  if (!is.null(background)) background_paths[[sample_ids[[1]]]] <- background
  if (!is.null(background_sheet)) background_paths <- read_background_sheet(background_sheet, sample_ids)

  workers <- min(worker_count(), length(sample_ids))
  backend <- "sequential"
  old_plan <- future::plan()
  old_max <- getOption("future.globals.maxSize")
  on.exit({ future::plan(old_plan); options(future.globals.maxSize = old_max) }, add = TRUE)
  options(future.globals.maxSize = future_max_size_for(object))
  if (workers > 1L) {
    future_parallel_plan(workers)
    backend <- if (.Platform$OS.type == "unix" && Sys.info()[["sysname"]] != "Windows") "multicore" else "multisession"
  } else {
    future::plan(future::sequential)
  }

  run_sample <- function(i) {
    sample_id <- sample_ids[[i]]
    cells <- colnames(counts)[sample_vector == sample_id]
    sample_counts <- counts[, cells, drop = FALSE]
    sample_seed <- as.integer(seed + i - 1L)
    background_info <- NULL
    background_sce <- NULL
    background_path <- background_paths[[sample_id]]
    if (!is.na(background_path) && nzchar(background_path)) {
      background_info <- prepare_background(
        background_path, sample_id, sample_counts,
        empty_only = background_is_empty_only
      )
      background_sce <- SingleCellExperiment::SingleCellExperiment(
        assays = list(counts = background_info$counts)
      )
    }
    sce <- SingleCellExperiment::SingleCellExperiment(assays = list(counts = sample_counts))
    fit <- if (is.null(background_sce)) {
      decontX::decontX(sce, seed = sample_seed, verbose = FALSE)
    } else {
      decontX::decontX(sce, background = background_sce, seed = sample_seed, verbose = FALSE)
    }
    if (!identical(colnames(fit), cells)) stop("DecontX changed cell order for sample ", sample_id, ".")
    contamination <- as.numeric(SummarizedExperiment::colData(fit)$decontX_contamination)
    corrected <- round(decontX::decontXcounts(fit))
    corrected <- validate_decontx_result(corrected, contamination, sample_counts, sample_id)
    clusters <- as.character(SummarizedExperiment::colData(fit)$decontX_clusters)
    if (length(clusters) != length(cells) || anyNA(clusters) || any(!nzchar(clusters))) {
      stop("Invalid DecontX broad-cluster labels for sample ", sample_id, ".")
    }
    list(
      sample_id = sample_id, cells = cells, corrected = corrected,
      contamination = contamination, clusters = clusters, seed = sample_seed,
      background_path = if (is.null(background_info)) "" else background_info$resolved_path,
      background_md5 = if (is.null(background_info)) "" else background_info$fingerprint,
      raw_droplets = if (is.null(background_info)) NA_integer_ else background_info$raw_droplets,
      called_removed = if (is.null(background_info)) NA_integer_ else background_info$called_removed,
      empty_droplets = if (is.null(background_info)) NA_integer_ else background_info$empty_droplets
    )
  }

  indices <- seq_along(sample_ids)
  results <- if (workers > 1L) {
    future.apply::future_lapply(indices, run_sample, future.seed = TRUE)
  } else {
    lapply(indices, run_sample)
  }
  names(results) <- sample_ids
  future::plan(future::sequential)

  corrected_counts <- do.call(cbind, lapply(results, `[[`, "corrected"))
  if (anyDuplicated(colnames(corrected_counts)) || !setequal(colnames(corrected_counts), colnames(object))) {
    stop("Per-sample DecontX outputs do not contain each input cell exactly once.")
  }
  corrected_counts <- corrected_counts[, colnames(object), drop = FALSE]
  corrected_counts <- validate_decontx_result(corrected_counts, rep(0, ncol(object)), counts, "merged object")

  contamination <- setNames(rep(NA_real_, ncol(object)), colnames(object))
  clusters <- setNames(rep(NA_character_, ncol(object)), colnames(object))
  for (result in results) {
    contamination[result$cells] <- result$contamination
    clusters[result$cells] <- paste(result$sample_id, result$clusters, sep = ":")
  }
  if (anyNA(contamination) || anyNA(clusters)) stop("Failed to map all DecontX results back to input cells.")

  object$decontX_contamination <- unname(contamination[colnames(object)])
  object$decontX_clusters <- unname(clusters[colnames(object)])
  SeuratObject::DefaultAssay(object) <- "RNA"

  sample_summary <- do.call(rbind, lapply(results, function(result) {
    x <- result$contamination
    data.frame(
      orig.ident = result$sample_id, n_cells = length(x),
      contamination_min = min(x), contamination_q05 = as.numeric(stats::quantile(x, 0.05)),
      contamination_median = stats::median(x), contamination_q95 = as.numeric(stats::quantile(x, 0.95)),
      contamination_max = max(x), seed = result$seed,
      background_path = result$background_path, background_md5 = result$background_md5,
      background_is_empty_only = background_is_empty_only, raw_droplets = result$raw_droplets,
      called_cells_removed_from_background = result$called_removed,
      empty_droplets_used = result$empty_droplets, stringsAsFactors = FALSE
    )
  }))
  marker_info <- marker_sets(rownames(counts))
  marker_files <- save_marker_plots(results, counts, corrected_counts, marker_info$used, marker_dir)

  decontx_provenance <- list(
    package_version = as.character(utils::packageVersion("decontX")),
    seed = seed, modeling_unit = "orig.ident", cluster_strategy = "DecontX internal broad clustering",
    corrected_counts_rounded = TRUE,
    parallel_backend = backend, workers = workers,
    requested_markers = marker_info$requested, used_markers = marker_info$used,
    per_sample = sample_summary
  )
  object@misc$decontX <- decontx_provenance
  object@misc$count_version <- list(
    version = "original", rna_counts_source = "observed_filtered_UMI_counts",
    parent_qc_rds = input_file, parent_qc_md5 = input_md5,
    scDblFinder_source = "original_observed_counts"
  )
  corrected_object <- build_corrected_object(object, corrected_counts, project)
  corrected_object@misc$decontX <- decontx_provenance
  corrected_object@misc$count_version <- list(
    version = "corrected", rna_counts_source = "rounded_decontX_corrected_counts",
    parent_qc_rds = input_file, parent_qc_md5 = input_md5,
    scDblFinder_source = "copied_from_original_observed_counts",
    count_dependent_qc_recalculated = c(
      "nCount_RNA", "nFeature_RNA", "percent_mito", "percent_ribo", "percent_hb",
      "S.Score", "G2M.Score", "Phase"
    )
  )
  if (!identical(object$scDblFinder.score, corrected_object$scDblFinder.score) ||
      !identical(as.character(object$scDblFinder.class), as.character(corrected_object$scDblFinder.class)) ||
      !identical(object$decontX_contamination, corrected_object$decontX_contamination) ||
      !identical(as.character(object$decontX_clusters), as.character(corrected_object$decontX_clusters))) {
    stop("Corrected branch did not preserve contamination or raw-count scDblFinder metadata exactly.")
  }
  cell_metadata <- data.frame(
    cell = colnames(object), orig.ident = as.character(object$orig.ident),
    decontX_contamination = object$decontX_contamination,
    decontX_clusters = object$decontX_clusters, stringsAsFactors = FALSE
  )

  atomic_save_rds(object, original_file)
  atomic_save_rds(corrected_object, corrected_file)
  atomic_write_csv(cell_metadata, cell_metadata_file)
  atomic_write_csv(sample_summary, sample_summary_file)
  atomic_write_lines(capture.output(sessionInfo()), session_file)
  manifest <- data.frame(
    count_version = c("original", "corrected"),
    branch_dir = c(original_dir, corrected_dir),
    rds = c(original_file, corrected_file),
    rds_md5 = c(file_md5(original_file), file_md5(corrected_file)),
    rna_counts_source = c("observed_filtered_UMI_counts", "rounded_decontX_corrected_counts"),
    n_cells = c(ncol(object), ncol(corrected_object)),
    n_genes = c(nrow(object), nrow(corrected_object)),
    cell_set_md5 = rep(cell_set_md5(colnames(object)), 2L),
    stringsAsFactors = FALSE
  )
  atomic_write_csv(manifest, manifest_file)
  count_version_review <- count_version_review_template(
    project, input_file, input_md5, cell_set_md5(colnames(object)),
    original_dir, corrected_dir, original_file, corrected_file
  )
  atomic_write_csv(count_version_review, count_version_review_file)
  record_workflow_event(
    state_file, project, stage = "decontx_compare", status = "waiting_count_version_comparison",
    script = "DecontX.r",
    inputs = list(
      qc_metrics_rds = input_file, qc_metrics_md5 = input_md5,
      background = background,
      background_sheet = background_sheet,
      background_sheet_md5 = if (is.null(background_sheet)) NULL else file_md5(background_sheet),
      background_is_empty_only = background_is_empty_only
    ),
    outputs = list(
      original_rds = original_file, original_md5 = file_md5(original_file),
      corrected_rds = corrected_file, corrected_md5 = file_md5(corrected_file),
      original_branch_dir = original_dir, corrected_branch_dir = corrected_dir,
      count_version_manifest = manifest_file, count_version_review = count_version_review_file,
      cell_metadata = cell_metadata_file, sample_summary = sample_summary_file,
      marker_plots = marker_files, session_info = session_file
    ),
    metrics = list(
      n_cells = ncol(object), n_genes = nrow(object), n_samples = length(sample_ids),
      modeling_unit = "orig.ident", parallel_backend = backend, workers = workers,
      corrected_counts_rounded = TRUE, first_round_requires_both_count_versions = TRUE,
      per_sample = sample_summary
    ),
    packages = packages
  )
  root_state <- read_workflow_state(state_file)
  write_workflow_state(root_state, original_state_file)
  write_workflow_state(root_state, corrected_state_file)
  for (branch in list(
    list(version = "original", file = original_file, state = original_state_file),
    list(version = "corrected", file = corrected_file, state = corrected_state_file)
  )) {
    record_workflow_event(
      branch$state, project, stage = "decontx", status = "ready_for_round1_clustering",
      script = "DecontX.r",
      inputs = list(qc_metrics_rds = input_file, qc_metrics_md5 = input_md5),
      outputs = list(
        decontx_rds = branch$file, decontx_md5 = file_md5(branch$file),
        original_rds = original_file, original_md5 = file_md5(original_file),
        corrected_rds = corrected_file, corrected_md5 = file_md5(corrected_file),
        count_version_review = count_version_review_file,
        cell_set_md5 = cell_set_md5(colnames(object))
      ),
      metrics = list(
        count_version = branch$version, n_cells = ncol(object), n_genes = nrow(object),
        scDblFinder_source = "original_observed_counts"
      ),
      packages = packages
    )
  }
  message("Created original-count branch: ", original_file)
  message("Created corrected-count branch: ", corrected_file)
  message("Complete Round 1 in both branch directories, then approve: ", count_version_review_file)
  message("Created per-cell metadata: ", cell_metadata_file)
  message("Created per-sample summary: ", sample_summary_file)
  message("Created marker plots in: ", marker_dir)
}, error = function(e) {
  record_failure(
    state_file, project, "decontx", NULL, "DecontX.r", e,
    inputs = list(
      qc_metrics_rds = input_file, background = background,
      background_sheet = background_sheet, background_is_empty_only = background_is_empty_only
    ),
    packages = packages
  )
  stop(e)
})
