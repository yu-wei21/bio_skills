#!/usr/bin/env Rscript
#
# 功能：将一个或多个 Cell Ranger Gene Expression 输出合并为未归一化的 Seurat v5 RDS。
# 输入：
#   1. --sample-sheet：CSV/TSV；前两列必须依次为 sample_id、matrix_path，每行一个独立 10x capture。
#   2. matrix_path：支持 filtered_feature_bc_matrix/、filtered_feature_bc_matrix.h5 或 Cell Ranger outs/。
#   3. 特征必须为 human/GRCh38 gene symbols；不支持 raw matrix、Cell Hashing、Flex 或 lane-split 配对。
#   4. --project 与 --output-dir：项目名和输出根目录；sample sheet 额外列会复制到细胞 metadata。
# 输出：
#   1. <output-dir>/<project>.raw.rds：仅含原始 counts；orig.ident 固定为 sample_id。
#   2. <output-dir>/workflow_state.json：自动维护的运行状态与 provenance。
# 示例命令：
#   Rscript scripts/create-seurat-object.R --sample-sheet samples.csv --project prostate --output-dir output

script_arg <- grep("^--file=", commandArgs(), value = TRUE)
script_dir <- if (length(script_arg)) dirname(normalizePath(sub("^--file=", "", script_arg[[1]]))) else getwd()
source(file.path(script_dir, "common.R"))
source(file.path(script_dir, "workflow-state.R"))

read_sample_sheet <- function(path) {
  path <- normalizePath(path, mustWork = TRUE)
  ext <- tolower(tools::file_ext(path))
  if (ext == "csv") {
    sheet <- utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  } else if (ext %in% c("tsv", "txt")) {
    sheet <- utils::read.delim(path, stringsAsFactors = FALSE, check.names = FALSE)
  } else {
    stop("Sample sheet must be .csv, .tsv, or .txt.")
  }
  if (ncol(sheet) < 2L || !identical(colnames(sheet)[1:2], c("sample_id", "matrix_path"))) {
    stop("The first two sample-sheet columns must be ordered as sample_id, matrix_path.")
  }
  if (!nrow(sheet) || anyNA(sheet$sample_id) || any(!nzchar(trimws(sheet$sample_id))) || anyDuplicated(sheet$sample_id)) {
    stop("sample_id must be non-empty and unique.")
  }
  if (anyNA(sheet$matrix_path) || any(!nzchar(trimws(sheet$matrix_path)))) stop("matrix_path cannot be empty.")
  reserved <- intersect(colnames(sheet)[-(1:2)], c("orig.ident", "nCount_RNA", "nFeature_RNA"))
  if (length(reserved)) stop("Sample-sheet extra columns cannot overwrite reserved metadata: ", paste(reserved, collapse = ", "))
  sheet_dir <- dirname(path)
  is_absolute <- grepl("^/|^[A-Za-z]:[/\\\\]", sheet$matrix_path)
  sheet$matrix_path[!is_absolute] <- file.path(sheet_dir, sheet$matrix_path[!is_absolute])
  sheet
}

resolve_matrix_path <- function(path) {
  path <- normalizePath(path, mustWork = TRUE)
  if (file.info(path)$isdir && basename(path) == "outs") {
    h5 <- file.path(path, "filtered_feature_bc_matrix.h5")
    dir <- file.path(path, "filtered_feature_bc_matrix")
    if (file.exists(h5)) return(h5)
    if (dir.exists(dir)) return(dir)
    stop("Cell Ranger outs/ lacks filtered_feature_bc_matrix.h5 or filtered_feature_bc_matrix/.")
  }
  path
}

read_gene_expression <- function(path) {
  value <- if (file.info(path)$isdir) Seurat::Read10X(path) else Seurat::Read10X_h5(path)
  if (is.list(value)) {
    if (!"Gene Expression" %in% names(value)) stop("10x file contains multiple modalities but no 'Gene Expression' matrix.")
    value <- value[["Gene Expression"]]
  }
  if (!inherits(value, "Matrix")) value <- methods::as(value, "dgCMatrix")
  value
}

opts <- parse_cli_args()
sample_sheet_file <- required_opt(opts, "sample_sheet")
project <- required_opt(opts, "project")
output_dir <- ensure_dir(required_opt(opts, "output_dir"))
output_file <- file.path(output_dir, paste0(project, ".raw.rds"))
state_file <- file.path(output_dir, "workflow_state.json")
packages <- c("Seurat", "SeuratObject", "Matrix", "jsonlite")

tryCatch({
  require_packages(packages, minimum = list(Seurat = "5.0.0", SeuratObject = "5.0.0"))
  sheet <- read_sample_sheet(sample_sheet_file)
  objects <- vector("list", nrow(sheet))
  raw_counts <- vector("list", nrow(sheet))

  for (i in seq_len(nrow(sheet))) {
    sample_id <- trimws(sheet$sample_id[[i]])
    matrix_path <- resolve_matrix_path(sheet$matrix_path[[i]])
    counts <- read_gene_expression(matrix_path)
    assert_human_gene_symbols(rownames(counts))
    raw_counts[[i]] <- c(n_genes = nrow(counts), n_cells = ncol(counts))
    colnames(counts) <- paste0(sample_id, "_", colnames(counts))
    object <- Seurat::CreateSeuratObject(
      counts = counts, project = project, assay = "RNA", min.cells = 0, min.features = 0
    )
    object$orig.ident <- sample_id
    extra <- setdiff(colnames(sheet), c("sample_id", "matrix_path"))
    for (column in extra) object[[column]] <- sheet[[column]][[i]]
    objects[[i]] <- object
  }

  object <- if (length(objects) == 1L) objects[[1]] else merge(objects[[1]], y = objects[-1], project = project)
  assert_seurat_input(object)
  if (!setequal(unique(as.character(object$orig.ident)), sheet$sample_id)) stop("Merged orig.ident values do not match the sample sheet.")
  for (sample_id in sheet$sample_id) {
    observed <- sum(as.character(object$orig.ident) == sample_id)
    expected <- raw_counts[[match(sample_id, sheet$sample_id)]][["n_cells"]]
    if (observed != expected) stop("Cell-count mismatch after merge for sample ", sample_id, ".")
  }

  atomic_save_rds(object, output_file)
  record_workflow_event(
    state_file, project, stage = "create", status = "ready_for_qc_metrics",
    script = "create-seurat-object.R",
    inputs = list(sample_sheet = normalizePath(sample_sheet_file, mustWork = TRUE), sample_sheet_md5 = file_md5(sample_sheet_file)),
    outputs = list(raw_rds = output_file, raw_rds_md5 = file_md5(output_file)),
    metrics = list(n_samples = nrow(sheet), n_cells = ncol(object), n_genes = nrow(object), per_sample = as.list(table(object$orig.ident))),
    packages = packages
  )
  message("Created: ", output_file)
}, error = function(e) {
  record_failure(state_file, project, "create", NULL, "create-seurat-object.R", e,
                 inputs = list(sample_sheet = sample_sheet_file), packages = packages)
  stop(e)
})
