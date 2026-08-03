#!/usr/bin/env Rscript
## 输入：
##   1. 已注释 Seurat RDS（含 celltype 列和可用于 DimPlot 的 UMAP；可用 --celltype-col 指定列名）
##   2. 亚群映射 CSV：第一列 celltype（原亚群名，须唯一且存在于对象），第二列 display_name
##      （含 marker 的展示名，如 AR+ / CD4_Tn_CCR7 风格）；行序即集群编号顺序
##      （第 1 行 = 编号 1，依此类推）；可选第三列 color（十六进制颜色，不填则用内置默认调色板）
## 输出：
##   1. output/plot-lineage-umap/ 下以输入 RDS 文件名为前缀的两页 PDF（第 1 页 UMAP 编号图，第 2 页编号与展示名图例）
##   2. output/plot-lineage-umap/ 下以输入 RDS 文件名为前缀的第 1 页对应 PNG
## 示例命令：Rscript plot-lineage-umap.r --input data/annotated.rds --mapping mapping.csv
##   Rscript plot-lineage-umap.r --input data/annotated.rds --mapping mapping.csv \
##         --celltype-col celltype --title "PCa subtypes" --output-dir output/custom
## 注意事项：对象中未出现在映射表中的 celltype 会以灰色 "Other" 显示并提示；
##           系统需提供 Ghostscript 命令 gs 用于合并 PDF。

library(Seurat)
library(patchwork)
library(ggplot2)

## ----------------------------- 参数解析 ---------------------------------
parse_args <- function(argv) {
  opts <- list(
    input = NULL, mapping = NULL, celltype_col = "celltype",
    title = NULL, output_dir = file.path("output", "plot-lineage-umap")
  )
  i <- 1L
  while (i <= length(argv)) {
    arg <- argv[i]
    if (arg %in% c("-h", "--help")) {
      cat(
        "Usage: Rscript plot-lineage-umap.r --input <annotated.rds> --mapping <mapping.csv>",
        " [--celltype-col COL] [--title TITLE] [--output-dir DIR]",
        sep = "\n"
      )
      quit(save = "no", status = 0)
    }
    if (!startsWith(arg, "--")) {
      stop("Unexpected argument: ", arg, call. = FALSE)
    }
    key <- substring(arg, 3L)
    if (!key %in% names(opts)) {
      stop("Unknown argument: --", key, call. = FALSE)
    }
    if (i == length(argv) || startsWith(argv[i + 1L], "--")) {
      stop("Missing value for argument: --", key, call. = FALSE)
    }
    opts[[key]] <- argv[i + 1L]
    i <- i + 2L
  }
  opts
}

opts <- parse_args(commandArgs(trailingOnly = TRUE))
if (is.null(opts$input) || is.null(opts$mapping)) {
  stop("--input 与 --mapping 为必需参数。\n",
       "Usage: Rscript plot-lineage-umap.r --input <annotated.rds> --mapping <mapping.csv>",
       call. = FALSE)
}
if (!file.exists(opts$input)) {
  stop("Input RDS does not exist: ", opts$input, call. = FALSE)
}
if (!file.exists(opts$mapping)) {
  stop("Mapping CSV does not exist: ", opts$mapping, call. = FALSE)
}

## ----------------------------- 读取映射表 --------------------------------
mapping <- read.csv(opts$mapping, check.names = FALSE, stringsAsFactors = FALSE)
if (ncol(mapping) < 2) {
  stop("Mapping CSV 至少需要两列: celltype, display_name", call. = FALSE)
}
celltype_map <- trimws(as.character(mapping[[1]]))
display_map  <- trimws(as.character(mapping[[2]]))
if (anyNA(celltype_map) || any(!nzchar(celltype_map))) {
  stop("Mapping CSV 第一列(celltype)不能含 NA 或空值", call. = FALSE)
}
if (anyDuplicated(celltype_map)) {
  stop("Mapping CSV 第一列(celltype)必须唯一: ",
       paste(unique(celltype_map[duplicated(celltype_map)]), collapse = ", "), call. = FALSE)
}
if (anyNA(display_map) || any(!nzchar(display_map))) {
  stop("Mapping CSV 第二列(display_name)不能含 NA 或空值", call. = FALSE)
}

## 行序 = 集群编号顺序（1..N）
cluster_num <- setNames(seq_along(celltype_map), celltype_map)
new_names   <- setNames(display_map, celltype_map)

## 颜色：优先用第三列 color，缺失时用内置默认调色板按行数循环
n_clusters <- length(celltype_map)
default_colors <- c(
  "#FDBF6F", "#6A3D9A", "#FF7F00", "#B2DF8A", "#FFFF99", "#CAB2D6",
  "#FB9A99", "#A6CEE3", "#E31A1C", "#33A02C", "#FC8D62", "#8DA0CB",
  "#E78AC3", "#66C2A5", "#A6761D", "#1F78B4", "#E7298A", "#66A61E",
  "#7570B3", "#D95F02", "#1B9E77", "#E6AB02"
)
if (ncol(mapping) >= 3 && all(nzchar(trimws(as.character(mapping[[3]]))))) {
  colors_used <- trimws(as.character(mapping[[3]]))
  if (anyNA(colors_used) || any(!nzchar(colors_used))) {
    stop("Mapping CSV 第三列(color)不能含 NA 或空值", call. = FALSE)
  }
} else {
  colors_used <- rep(default_colors, length.out = n_clusters)
}
names(colors_used) <- celltype_map

## ----------------------------- 读取对象并加标签 ----------------------------
message("Reading: ", opts$input)
obj <- readRDS(opts$input)
if (!inherits(obj, "Seurat")) {
  stop("Input RDS is not a Seurat object.", call. = FALSE)
}
if (!opts$celltype_col %in% colnames(obj@meta.data)) {
  stop("meta.data 中未找到 celltype 列: ", opts$celltype_col, call. = FALSE)
}
ct_values <- as.character(obj[[opts$celltype_col]][, 1])
if (anyNA(ct_values) || any(!nzchar(ct_values))) {
  stop(opts$celltype_col, " 列包含 NA 或空值", call. = FALSE)
}

missing_map <- setdiff(celltype_map, unique(ct_values))
if (length(missing_map) > 0) {
  stop("映射表中的 celltype 在对象中不存在: ",
       paste(missing_map, collapse = ", "), call. = FALSE)
}
unmapped_obs <- setdiff(unique(ct_values), celltype_map)
if (length(unmapped_obs) > 0) {
  message("对象中存在但未在映射表中的 celltype（将以灰色 Other 显示）: ",
          paste(unmapped_obs, collapse = ", "))
}

## 集群编号（字符型，1..N）；未映射的细胞标为 "Other"
obj$cluster_num <- as.character(unname(cluster_num[ct_values]))
obj$cluster_num[is.na(obj$cluster_num)] <- "Other"
## 因子水平 = 编号 1..N + Other（保证 DimPlot 的 cols 一一对应）
level_order <- c(as.character(seq_along(celltype_map)), "Other")
obj$cluster_num <- factor(obj$cluster_num, levels = level_order)

## 颜色按因子水平（编号 1..N + Other）命名，DimPlot 才能正确匹配
plot_colors <- setNames(c(unname(colors_used), "grey75"), level_order)

## ----------------------------- 绘图 --------------------------------------
title_text <- if (is.null(opts$title)) "" else opts$title

p_umap <- DimPlot(obj, group.by = "cluster_num", cols = plot_colors,
                  label = TRUE, label.size = 6, repel = TRUE) +
  ggtitle(title_text) +
  NoLegend() +
  theme(
    plot.title   = element_text(hjust = 0.5, size = 14),
    aspect.ratio = 1,
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5)
  )

## 图例页：编号 + 展示名 + 颜色（按编号 1..N 顺序）
df_legend <- data.frame(
  num   = seq_along(celltype_map),
  label = paste0(seq_along(celltype_map), "  ", unname(new_names)),
  color = unname(colors_used),
  y     = rev(seq_along(celltype_map)), ## 编号 1 在顶部
  stringsAsFactors = FALSE
)
p_legend <- ggplot(df_legend, aes(x = 1, y = y)) +
  geom_point(color = df_legend$color, size = 4) +
  geom_text(aes(label = label), hjust = 0, nudge_x = 0.15, size = 4.5) +
  labs(title = title_text) +
  xlim(0.5, 3.5) +
  theme_void() +
  theme(
    plot.title  = element_text(hjust = 0, size = 12, face = "bold",
                               margin = margin(b = 2, t = 6)),
    plot.margin = margin(4, 4, 4, 4)
  )

## ----------------------------- 保存输出 ----------------------------------
dir.create(opts$output_dir, showWarnings = FALSE, recursive = TRUE)
input_tag <- tools::file_path_sans_ext(basename(opts$input))
output_prefix <- file.path(opts$output_dir, paste0(input_tag, "_lineage_umap"))

tmp1 <- tempfile(fileext = ".pdf")
tmp2 <- tempfile(fileext = ".pdf")
outfile_pdf <- paste0(output_prefix, ".pdf")

ggsave(tmp1, p_umap, width = 210, height = 210, units = "mm")
ggsave(tmp2, p_legend, width = 210, height = 297, units = "mm")
ggsave(paste0(output_prefix, ".png"), p_umap, width = 8, height = 8, dpi = 300, bg = "white")

## 合并两页 PDF（需要系统 Ghostscript gs）
gs_cmd <- paste("gs -q -dNOPAUSE -dBATCH -sDEVICE=pdfwrite",
  paste0("-sOutputFile=", outfile_pdf),
  tmp1, tmp2)
if (Sys.which("gs") == "") {
  warning("未找到 Ghostscript (gs)，跳过两页 PDF 合并；单页 PDF 仍在临时目录: ",
          tmp1, ", ", tmp2)
} else {
  system(gs_cmd, intern = TRUE)
  message("Saved: ", outfile_pdf)
}

unlink(c(tmp1, tmp2))
message("Saved: ", paste0(output_prefix, ".png"))
