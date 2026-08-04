#!/usr/bin/env Rscript
#
# 功能：从 RNA counts 重算基础 QC，计算细胞周期，并按 orig.ident 运行 scDblFinder；输出交给 DecontX.r。
# 输入：
#   1. --input：含 RNA counts 的 Seurat v5 RDS；orig.ident 必须非空，物种为 human/GRCh38 gene symbols。
#   2. --project 与 --output-dir：项目名和输出根目录；输入可来自创建脚本或外部 counts-bearing Seurat RDS。
# 输出：
#   1. <output-dir>/<project>.qc-metrics.rds：metadata 重算 nCount_RNA、nFeature_RNA，并新增 percent_mito、percent_ribo、percent_hb、S.Score、G2M.Score、Phase、scDblFinder.score、scDblFinder.class。
#   2. <output-dir>/workflow_state.json：记录参数、包版本和每样本 doublet 摘要，状态转为 ready_for_decontx；不生成独立 doublet 文件或图。
# 示例命令：
#   Rscript scripts/seurat-qc-metrics.R --input output/prostate.raw.rds --project prostate --output-dir output

script_arg <- grep("^--file=", commandArgs(), value = TRUE)
script_dir <- if (length(script_arg)) dirname(normalizePath(sub("^--file=", "", script_arg[[1]]))) else getwd()
source(file.path(script_dir, "common.R"))
source(file.path(script_dir, "workflow-state.R"))
source(file.path(script_dir, "scdblfinder-utils.R"))

opts <- parse_cli_args()
input_file <- normalizePath(required_opt(opts, "input"), mustWork = TRUE)
project <- required_opt(opts, "project")
output_dir <- ensure_dir(required_opt(opts, "output_dir"))
output_file <- file.path(output_dir, paste0(project, ".qc-metrics.rds"))
state_file <- file.path(output_dir, "workflow_state.json")
assert_output_not_input(input_file, output_file)
packages <- c("Seurat", "SeuratObject", "future", "future.apply", "parallelly", "scDblFinder", "SingleCellExperiment", "BiocParallel", "S4Vectors", "SummarizedExperiment", "Matrix", "jsonlite")

tryCatch({
  require_packages(packages, minimum = list(Seurat = "5.0.0", SeuratObject = "5.0.0"))
  set.seed(WORKFLOW_SEED)
  object <- readRDS(input_file)
  assert_seurat_input(object)
  object <- join_rna_layers(object)
  counts <- rna_layer(object, "counts")
  object$nCount_RNA <- as.numeric(Matrix::colSums(counts))
  object$nFeature_RNA <- as.numeric(Matrix::colSums(counts > 0))
  assert_positive_qc_counts(object)
  SeuratObject::DefaultAssay(object) <- "RNA"

  object <- Seurat::PercentageFeatureSet(object, pattern = "^MT-", col.name = "percent_mito")
  object <- Seurat::PercentageFeatureSet(object, pattern = "^RP[SL]", col.name = "percent_ribo")
  object <- Seurat::PercentageFeatureSet(object, pattern = "^HB[ABDEGMQZ]", col.name = "percent_hb")

  old_plan <- future::plan()
  old_max <- getOption("future.globals.maxSize")
  on.exit({ future::plan(old_plan); options(future.globals.maxSize = old_max) }, add = TRUE)
  options(future.globals.maxSize = future_max_size_for(object))
  future_parallel_plan()
  object <- Seurat::NormalizeData(object, normalization.method = "LogNormalize", scale.factor = 10000, verbose = FALSE)
  future::plan(future::sequential)

  cc <- Seurat::cc.genes.updated.2019
  s_genes <- intersect(cc$s.genes, rownames(object))
  g2m_genes <- intersect(cc$g2m.genes, rownames(object))
  if (length(s_genes) < 5L || length(g2m_genes) < 5L) stop("Too few human cell-cycle genes were found for reliable scoring.")
  object <- Seurat::CellCycleScoring(object, s.features = s_genes, g2m.features = g2m_genes, set.ident = FALSE)
  object <- run_scdblfinder(object)

  keep <- c("scDblFinder.score", "scDblFinder.class")
  unexpected <- grep("^scDblFinder\\.", colnames(object[[]]), value = TRUE)
  for (column in setdiff(unexpected, keep)) object[[column]] <- NULL
  atomic_save_rds(object, output_file)
  record_workflow_event(
    state_file, project, stage = "qc_metrics", status = "ready_for_decontx",
    script = "seurat-qc-metrics.R",
    inputs = list(input_rds = input_file, input_md5 = file_md5(input_file)),
    outputs = list(qc_metrics_rds = output_file, qc_metrics_md5 = file_md5(output_file)),
    metrics = list(
      n_cells = ncol(object), n_genes = nrow(object),
      nCount_RNA_recalculated = TRUE, nFeature_RNA_recalculated = TRUE,
      scDblFinder = object@misc$scDblFinder
    ),
    packages = packages
  )
  message("Created: ", output_file)
}, error = function(e) {
  record_failure(state_file, project, "qc_metrics", NULL, "seurat-qc-metrics.R", e,
                 inputs = list(input_rds = input_file), packages = packages)
  stop(e)
})
