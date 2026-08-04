#!/usr/bin/env Rscript
#
# 功能：以固定 LogNormalize 流程重新降维聚类；每轮过滤后必须重新运行本脚本。
# 输入：
#   1. --input：DecontX.r 产生的 original 或 corrected count-version Seurat v5 RDS；RNA assay 必须是该分支实际用于分析的 counts，且含完整 QC metadata。
#   2. --round：1、2 或 3；--project 与 --output-dir 指定项目名和输出根目录。
#   3. 固定参数：seed=666、2000 HVGs、30 PCs、dims 1:30、Louvain algorithm=1、resolution 0.1–1.5/0.1。
# 输出：
#   1. <output-dir>/roundN/<project>.roundN.clustered.rds：15 个 resolution 列及 UMAP；多 orig.ident 在标准化/整合前按样本 split layers，Harmony 后 JoinLayers。
#   2. <output-dir>/workflow_state.json：记录本轮输入输出、Harmony 使用情况和包版本。
# 示例命令：
#   Rscript scripts/seurat-cluster.R --input output/count-versions/original/prostate.decontx-original.rds --project prostate --round 1 --output-dir output/count-versions/original

script_arg <- grep("^--file=", commandArgs(), value = TRUE)
script_dir <- if (length(script_arg)) dirname(normalizePath(sub("^--file=", "", script_arg[[1]]))) else getwd()
source(file.path(script_dir, "common.R"))
source(file.path(script_dir, "workflow-state.R"))

opts <- parse_cli_args()
input_file <- normalizePath(required_opt(opts, "input"), mustWork = TRUE)
project <- required_opt(opts, "project")
round <- parse_integer(required_opt(opts, "round"), "round", 1L, 3L)
paths <- project_paths(required_opt(opts, "output_dir"), project, round)
output_file <- paths$clustered
assert_output_not_input(input_file, output_file)
base_packages <- c("Seurat", "SeuratObject", "future", "future.apply", "parallelly", "Matrix", "jsonlite")

tryCatch({
  require_packages(base_packages, minimum = list(Seurat = "5.0.0", SeuratObject = "5.0.0"))
  assert_round_clustering_input(paths$state, project, round, input_file)
  set.seed(WORKFLOW_SEED)
  object <- readRDS(input_file)
  assert_seurat_input(object)
  input_cell_set_md5 <- cell_set_md5(colnames(object))
  needed <- c("S.Score", "G2M.Score", "Phase", "percent_mito", "scDblFinder.class", "decontX_contamination")
  missing <- setdiff(needed, colnames(object[[]]))
  if (length(missing)) stop("Run seurat-qc-metrics.R first; missing metadata: ", paste(missing, collapse = ", "))
  count_version <- object@misc$count_version$version %||% ""
  if (!count_version %in% c("original", "corrected")) {
    stop("Run DecontX.r first; object@misc$count_version$version must be original or corrected.")
  }
  if (anyNA(object$decontX_contamination) || any(!is.finite(object$decontX_contamination)) ||
      any(object$decontX_contamination < 0 | object$decontX_contamination > 1)) {
    stop("decontX_contamination must contain finite values in [0, 1].")
  }
  object <- join_rna_layers(object, require_data = TRUE)
  assert_positive_qc_counts(object)
  object <- clean_for_recluster(object)
  object@misc$filter_provenance <- NULL
  SeuratObject::DefaultAssay(object) <- "RNA"

  samples <- as.character(object$orig.ident)
  use_harmony <- length(unique(samples)) > 1L
  packages <- base_packages
  if (use_harmony) {
    require_packages("harmony")
    packages <- c(packages, "harmony")
    object[["RNA"]] <- split(object[["RNA"]], f = samples)
    assert_split_rna_layers(object, samples, require_data = TRUE)
  }

  old_plan <- future::plan()
  old_max <- getOption("future.globals.maxSize")
  on.exit({ future::plan(old_plan); options(future.globals.maxSize = old_max) }, add = TRUE)
  options(future.globals.maxSize = future_max_size_for(object))
  future_parallel_plan()
  object <- Seurat::NormalizeData(object, normalization.method = "LogNormalize", scale.factor = 10000, verbose = FALSE)
  future::plan(future::sequential)
  object <- Seurat::FindVariableFeatures(object, selection.method = "vst", nfeatures = WORKFLOW_NFEATURES, verbose = FALSE)
  Seurat::VariableFeatures(object) <- setdiff(Seurat::VariableFeatures(object), cell_cycle_genes())
  if (length(Seurat::VariableFeatures(object)) < 30L) stop("Fewer than 30 variable features remain after cell-cycle gene removal.")
  future_parallel_plan()
  object <- Seurat::ScaleData(object, features = Seurat::VariableFeatures(object), vars.to.regress = c("S.Score", "G2M.Score"), verbose = FALSE)
  future::plan(future::sequential)
  object <- Seurat::RunPCA(object, features = Seurat::VariableFeatures(object), npcs = 30, seed.use = WORKFLOW_SEED, verbose = FALSE)

  reduction <- "pca"
  if (use_harmony) {
    object <- Seurat::IntegrateLayers(
      object = object, method = Seurat::HarmonyIntegration,
      orig.reduction = "pca", new.reduction = "harmony", verbose = FALSE
    )
    reduction <- "harmony"
  }
  object <- Seurat::RunUMAP(object, reduction = reduction, dims = WORKFLOW_DIMS,
                            seed.use = WORKFLOW_SEED, verbose = FALSE)
  object <- Seurat::FindNeighbors(object, reduction = reduction, dims = WORKFLOW_DIMS, verbose = FALSE)
  future_parallel_plan()
  object <- Seurat::FindClusters(
    object, resolution = WORKFLOW_RESOLUTIONS, algorithm = 1,
    random.seed = WORKFLOW_SEED, verbose = FALSE
  )
  future::plan(future::sequential)
  object <- join_rna_layers(object, require_data = TRUE)
  object@misc$major_annotation <- list(
    round = round, seed = WORKFLOW_SEED, normalization = "LogNormalize",
    nfeatures = WORKFLOW_NFEATURES, dims = WORKFLOW_DIMS,
    resolutions = WORKFLOW_RESOLUTIONS, reduction = reduction,
    cell_cycle_regressed = c("S.Score", "G2M.Score"),
    cell_cycle_genes_removed_from_hvg = TRUE
  )

  atomic_save_rds(object, output_file)
  record_workflow_event(
    paths$state, project, stage = "clustering", status = "waiting_resolution_review", round = round,
    script = "seurat-cluster.R",
    inputs = list(input_rds = input_file, input_md5 = file_md5(input_file), input_cells_md5 = input_cell_set_md5),
    outputs = list(clustered_rds = output_file, clustered_md5 = file_md5(output_file), cell_set_md5 = cell_set_md5(colnames(object))),
    metrics = list(
      n_cells = ncol(object), n_genes = nrow(object), n_samples = length(unique(samples)),
      count_version = count_version, harmony = use_harmony,
      rna_layers_joined_before_save = TRUE
    ),
    packages = packages
  )
  message("Created: ", output_file)
}, error = function(e) {
  record_failure(paths$state, project, "clustering", round, "seurat-cluster.R", e,
                 inputs = list(input_rds = input_file), packages = base_packages)
  stop(e)
})
