#!/usr/bin/env Rscript
# 输入：
#   1. 已注释 Seurat RDS：metadata 含由 --sample-col、--group-col、--celltype-col 指定的列；
#      每个样本必须仅属于一个分组，且三列均无 NA 或空字符串。
# 输出：
#   1. output/propeller/ 下以输入 RDS 文件名为前缀的细胞类型比例差异检验 CSV。
#   2. 同前缀的 sample×celltype 比例 CSV、每组样本数 CSV 和 sessionInfo TXT。
# 示例命令：
#   Rscript propeller.r data/annotated.rds --sample-col orig.ident --group-col condition \\
#     --celltype-col celltype
# 注意事项：检验单位是生物学样本而不是细胞；每个比较组至少需要 2 个独立生物学重复。
#           不应将技术重复或通过 bootstrap 构造的伪重复当作真实重复；多次采样请另建模型。

usage <- paste(
  "Usage:",
  "  Rscript propeller.r <input.rds> --sample-col COLUMN --group-col COLUMN",
  "    --celltype-col COLUMN [--transform logit|asin] [--output-dir DIR]",
  sep = "\n"
)

parse_args <- function(args) {
  opts <- list(
    input = NULL,
    sample_col = NULL,
    group_col = NULL,
    celltype_col = NULL,
    transform = "logit",
    output_dir = file.path("output", "propeller")
  )
  positional <- character(0)
  i <- 1L
  while (i <= length(args)) {
    arg <- args[[i]]
    if (arg %in% c("-h", "--help")) {
      cat(usage, "\n")
      quit(save = "no", status = 0)
    }
    if (!startsWith(arg, "--")) {
      positional <- c(positional, arg)
      i <- i + 1L
      next
    }
    key <- gsub("-", "_", sub("^--", "", arg))
    if (!key %in% names(opts) || key == "input") {
      stop("未知参数: ", arg, "\n", usage, call. = FALSE)
    }
    if (i == length(args)) stop("参数缺少取值: ", arg, "\n", usage, call. = FALSE)
    opts[[key]] <- args[[i + 1L]]
    i <- i + 2L
  }
  if (length(positional) != 1L) stop("必须提供且仅提供一个输入 RDS。\n", usage, call. = FALSE)
  opts$input <- positional[[1]]
  required <- c("sample_col", "group_col", "celltype_col")
  if (anyNA(unlist(opts[required])) || any(!nzchar(unlist(opts[required])))) {
    stop("--sample-col、--group-col 和 --celltype-col 均为必填参数。\n", usage, call. = FALSE)
  }
  if (!opts$transform %in% c("logit", "asin")) {
    stop("--transform 只能为 logit 或 asin。", call. = FALSE)
  }
  opts
}

require_packages <- function(packages) {
  missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing) > 0) {
    stop(
      "缺少 R 包: ", paste(missing, collapse = ", "),
      ". 请按 Bioconductor 官方方式安装 speckle 及其依赖后重试。",
      call. = FALSE
    )
  }
}

opts <- parse_args(commandArgs(trailingOnly = TRUE))
require_packages(c("Seurat", "speckle"))
if (!file.exists(opts$input)) stop("输入 RDS 不存在: ", opts$input, call. = FALSE)
if (!nzchar(opts$output_dir)) stop("--output-dir 不能为空。", call. = FALSE)

seurat_obj <- readRDS(opts$input)
if (!inherits(seurat_obj, "Seurat")) stop("输入 RDS 不是 Seurat 对象。", call. = FALSE)
meta <- seurat_obj[[]]
required_cols <- c(opts$sample_col, opts$group_col, opts$celltype_col)
missing_cols <- setdiff(required_cols, colnames(meta))
if (length(missing_cols) > 0L) {
  stop("metadata 缺少列: ", paste(missing_cols, collapse = ", "), call. = FALSE)
}

cell_info <- data.frame(
  sample = as.character(meta[[opts$sample_col]]),
  group = as.character(meta[[opts$group_col]]),
  celltype = as.character(meta[[opts$celltype_col]]),
  stringsAsFactors = FALSE
)
if (anyNA(cell_info) || any(!nzchar(as.matrix(cell_info)))) {
  stop("样本、分组和细胞类型列均不能含 NA 或空字符串。", call. = FALSE)
}
sample_group_n <- aggregate(group ~ sample, data = cell_info, FUN = function(x) length(unique(x)))
if (any(sample_group_n$group != 1L)) {
  bad_samples <- sample_group_n$sample[sample_group_n$group != 1L]
  stop("以下样本对应多个分组: ", paste(bad_samples, collapse = ", "), call. = FALSE)
}
sample_groups <- unique(cell_info[c("sample", "group")])
group_n <- as.data.frame(table(sample_groups$group), stringsAsFactors = FALSE)
colnames(group_n) <- c("group", "n_samples")
if (nrow(group_n) < 2L) stop("至少需要两个分组。", call. = FALSE)
if (any(group_n$n_samples < 2L)) {
  stop("每个分组至少需要 2 个独立生物学样本；当前样本数：",
       paste(paste0(group_n$group, "=", group_n$n_samples), collapse = ", "), call. = FALSE)
}

result <- speckle::propeller(
  clusters = cell_info$celltype,
  sample = cell_info$sample,
  group = cell_info$group,
  transform = opts$transform
)
result <- data.frame(celltype = rownames(result), result, row.names = NULL, check.names = FALSE)

count_table <- xtabs(~ celltype + sample, data = cell_info)
prop_table <- prop.table(count_table, margin = 2L)
prop_result <- data.frame(
  celltype = rownames(prop_table),
  as.data.frame.matrix(prop_table, check.names = FALSE),
  row.names = NULL,
  check.names = FALSE
)

input_tag <- tools::file_path_sans_ext(basename(opts$input))
dir.create(opts$output_dir, recursive = TRUE, showWarnings = FALSE)
output_prefix <- file.path(opts$output_dir, paste0(input_tag, "_propeller"))
utils::write.csv(result, paste0(output_prefix, "_results.csv"), row.names = FALSE)
utils::write.csv(prop_result, paste0(output_prefix, "_sample_proportions.csv"), row.names = FALSE)
utils::write.csv(group_n, paste0(output_prefix, "_group_sample_counts.csv"), row.names = FALSE)
writeLines(capture.output(sessionInfo()), paste0(output_prefix, "_sessionInfo.txt"))

message("propeller 完成：", output_prefix)
