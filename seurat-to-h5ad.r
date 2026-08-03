#!/usr/bin/env Rscript
# 输入：
#   1. 已完成 NormalizeData 的 Seurat v4/v5 RDS
# 输出：
#   1. 默认在输入 RDS 同目录生成的同名 .h5ad 文件，或由 --output 指定的文件
# 示例命令：Rscript seurat-to-h5ad.r --input data/seurat.rds --output output/seurat.h5ad --assay RNA
# 注意事项：需由 reticulate 绑定一个安装 anndata 的 Python；导出的 adata.X 是 data layer
#           （归一化且通常对数转换后的矩阵），不是原始 counts。
#
# 通用脚本：Seurat RDS -> h5ad（reticulate + anndata 直连，不依赖 sceasy）
#
# 特点：
#   - Seurat v4 (legacy Assay) 与 v5 (Assay5, 多 layer) 均兼容；
#     Assay5 多 layer 先 JoinLayers 合并再导出
#   - adata.X = data layer（归一化且通常经对数转换的数据，infercnvpy/scanpy 可直接使用）
#   - 默认导出 DefaultAssay（通常 RNA），--assay 可指定
#   - obs = meta.data（因子列转字符）；var = assay 的 feature 元数据
#   - Python 绑定顺序：RETICULATE_PYTHON -> CONDA_PREFIX -> conda env "seurat"，
#     避免 reticulate 自动探测到错误的解释器
#
# 用法：
#   Rscript seurat-to-h5ad.r --input <seurat.rds> [--output <out.h5ad>] [--assay <name>]
#
# 依赖：
#   R: Seurat, reticulate, Matrix
#   Python (conda env seurat): anndata

suppressPackageStartupMessages({
  library(Seurat)
  library(reticulate)
  library(Matrix)
})

# ---------- 参数解析 ----------

parse_args <- function(argv) {
  args <- list(input = NULL, output = NULL, assay = NULL)
  i <- 1
  while (i <= length(argv)) {
    a <- argv[i]
    if (a %in% c("--input", "--output", "--assay")) {
      if (i + 1 > length(argv)) stop("Missing value for ", a)
      args[[sub("^--", "", a)]] <- argv[i + 1]
      i <- i + 2
    } else {
      stop("Unknown argument: ", a)
    }
  }
  args
}

# ---------- Python 绑定 ----------

# 绑定流程：先探测候选 Python（RETICULATE_PYTHON -> CONDA_PREFIX/bin/python ->
# conda env "seurat"）中哪个装有 anndata，确定后才调用一次 use_python ——
# reticulate 一旦初始化某个 python 便无法在同一会话切换，所以探测用 system2 直查，
# 不经过 reticulate。conda env 定位不依赖 PATH 里的 conda 命令，
# 改用 ~/.conda/environments.txt 与常见 envs 目录。

py_has_anndata <- function(py) {
  suppressWarnings(
    system2(py, c("-c", shQuote("import anndata")), stdout = FALSE, stderr = FALSE) == 0
  )
}

find_env_python <- function(env = "seurat") {
  hits <- character(0)
# 1) conda environments.txt（conda 官方记录所有 env 路径，格式为裸目录）
  envs_txt <- file.path(path.expand("~"), ".conda", "environments.txt")
  if (file.exists(envs_txt)) {
    hits <- readLines(envs_txt, warn = FALSE)
    hits <- hits[endsWith(hits, paste0("/", env))]
  }
  # 2) 常见 envs 目录
  envs_dirs <- c(
    file.path(path.expand("~"), ".conda", "envs"),
    file.path(path.expand("~"), "miniconda3", "envs"),
    file.path(path.expand("~"), "anaconda3", "envs"),
    file.path(path.expand("~"), "mambaforge", "envs"),
    "/opt/conda/envs"
  )
  hits <- c(hits, file.path(envs_dirs, env, "bin", "python"))
  # 统一补全为 bin/python 路径
  hits <- ifelse(
    endsWith(hits, "bin/python"),
    hits,
    file.path(hits, "bin", "python")
  )
  hits <- unique(hits[nzchar(hits) & file.exists(hits)])
  if (length(hits) == 0) NULL else hits[1]
}

bind_python <- function(env = "seurat") {
  py_env <- Sys.getenv("RETICULATE_PYTHON")
  candidates <- c(
    if (nzchar(py_env)) py_env else character(0),
    file.path(Sys.getenv("CONDA_PREFIX"), "bin", "python"),
    find_env_python(env)
  )
  candidates <- unique(candidates[nzchar(candidates) & file.exists(candidates)])
  for (py in candidates) {
    if (py_has_anndata(py)) {
      message("Python: ", py)
      reticulate::use_python(py, required = TRUE)
      return(invisible(TRUE))
    }
  }
  stop(
    "No python with anndata found (tried: ",
    if (length(candidates) == 0) "(none located)" else paste(candidates, collapse = "; "),
    "). Install anndata (pip install anndata) in conda env '", env,
    "' or set RETICULATE_PYTHON to a python with anndata."
  )
}

# ---------- 主流程 ----------

args <- parse_args(commandArgs(trailingOnly = TRUE))
input_file <- args$input
if (is.null(input_file) || !nzchar(input_file)) {
  stop("Usage: Rscript seurat-to-h5ad.r --input <seurat.rds> [--output <out.h5ad>] [--assay <name>]")
}
if (!file.exists(input_file)) {
  stop("Input file not found: ", input_file)
}

output_file <- args$output
if (is.null(output_file) || !nzchar(output_file)) {
  output_file <- paste0(tools::file_path_sans_ext(input_file), ".h5ad")
}
dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)

message("Input : ", input_file)
message("Output: ", output_file)

message("1. Read Seurat object ...")
seurat <- readRDS(input_file)
if (!inherits(seurat, "Seurat")) {
  stop("Input object is not a Seurat object: ", class(seurat)[1])
}

# 确定要导出的 assay（默认 DefaultAssay，可指定）
if (!is.null(args$assay)) {
  if (!args$assay %in% names(seurat@assays)) {
    stop(
      "Assay '", args$assay, "' not found. Available assays: ",
      paste(names(seurat@assays), collapse = ", ")
    )
  }
  DefaultAssay(seurat) <- args$assay
}
assay_name <- DefaultAssay(seurat)
message("Assay: ", assay_name)
assay_obj <- seurat[[assay_name]]

# Assay5 多 layer 合并（merge/双样本产物会有 data.1, data.2 ... split layers）
if (inherits(assay_obj, "Assay5")) {
  layer_names <- Layers(assay_obj)
  if (length(layer_names) > 1) {
    message("  Join layers for assay '", assay_name, "': ",
            paste(layer_names, collapse = ", "))
    seurat <- JoinLayers(seurat, assay = assay_name)
    assay_obj <- seurat[[assay_name]]
  }
}

# data 矩阵（基因×细胞），转置为细胞×基因（anndata 标准方向）
message("2. Extract data matrix (layer = data, genes x cells) ...")
data_mat <- GetAssayData(assay_obj, layer = "data")
if (is.null(data_mat) || length(data_mat) == 0) {
  stop("Assay '", assay_name, "' has no 'data' layer (run NormalizeData first).")
}
message("   Dim: ", nrow(data_mat), " x ", ncol(data_mat),
        " | sparse: ", is(data_mat, "dgCMatrix") || is(data_mat, "dgTMatrix"))
data_mat_t <- Matrix::t(data_mat)  # dgCMatrix, cells x genes

# obs：细胞级 metadata
message("3. Build obs / var ...")
obs_df <- seurat@meta.data
# 因子列转字符，避免 reticulate 转换问题
for (col in colnames(obs_df)) {
  if (is.factor(obs_df[[col]])) {
    obs_df[[col]] <- as.character(obs_df[[col]])
  }
}

# var：基因级 metadata（Assay5 用 [[]]，Assay 用 @meta.features）
if (inherits(assay_obj, "Assay5")) {
  var_df <- assay_obj[[]]
} else {
  var_df <- assay_obj@meta.features
}
if (is.null(var_df) || nrow(var_df) == 0) {
  var_df <- data.frame(
    feature_name = rownames(assay_obj),
    row.names = rownames(assay_obj),
    stringsAsFactors = FALSE
  )
}
if (anyDuplicated(rownames(var_df))) {
  stop(
    "Duplicated feature names in var (", anyDuplicated(rownames(var_df)), " duplicates); ",
    "anndata requires a unique var index."
  )
}

# 写 h5ad（gzip 压缩，最大兼容性）
message("4. Write h5ad ...")
bind_python()
anndata <- import("anndata", convert = FALSE)
adata <- anndata$AnnData(
  X   = data_mat_t,
  obs = obs_df,
  var = var_df
)
adata$write_h5ad(output_file, compression = "gzip")

message("Done: ", output_file)
