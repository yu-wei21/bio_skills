## 输入：Seurat RDS（RNA data/data.*层，feature名为gene symbol）、无缺失的细胞类型metadata列名、含gene_name/gene_biotype列的注释CSV。
## 输出：output/presto/下以输入RDS文件名开头的完整DEG/ORA XLSX、背景CSV及逐细胞类型GSEA CSV。
## 示例：Rscript presto.r cancer.rds celltype gtf.csv

## 脚本目的：对 Seurat 对象中指定细胞类型列进行 presto 差异分析
##           匹配 protein-coding gene，输出差异基因列表
##           同时准备 ORA (top100 logFC) 和 GSEA (全基因 ranked list) 的输入文件
##
## 后续分析：
##   ORA:  Rscript ORA.r <输出_xlsx> SYMBOL <输出_universe.csv> --gmt-dir=<GMT目录>
##   GSEA: Rscript GSEA.r <输出_GSEA_csv> <GMT目录>
##
## 注意事项：
##   - presto v1.0.0 的 Seurat 方法使用 GetAssayData(slot=) 已废弃于 Seurat v5，
##     因此这里手动提取矩阵，调用 wilcoxauc.default 绕过兼容性问题。
##   - 全基因矩阵可能占用较多内存；脚本在释放大对象后显式 gc()。

suppressPackageStartupMessages({
  library(Seurat)
  library(presto)
  library(tidyverse)
  library(openxlsx)
  library(data.table)
})

## ---- 参数解析 ----
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3) {
  stop("Usage: Rscript presto.r <input.rds> <celltype_col> <gtf.csv>")
}
rds_file <- args[1]
celltype_col <- args[2]
gtf_file <- args[3]

## ---- ORA marker 预设 ----
marker_padj_max <- 0.05
marker_logfc_min <- 0.25
marker_auc_min <- 0.60
marker_pct_in_min <- 10
marker_n_max <- 100

## Excel sheet 名最多 31 字符，且不允许 \ / : * ? [ ]
## 清洗和截断后仍需保证名称唯一（Excel sheet 名大小写不敏感）
make_excel_sheet_names <- function(x) {
  base_names <- as.character(x)
  for (invalid_char in c("\\", "/", ":", "*", "?", "[", "]")) {
    base_names <- gsub(invalid_char, "_", base_names, fixed = TRUE)
  }
  base_names <- trimws(base_names)
  base_names[base_names == ""] <- "Sheet"

  result <- character(length(base_names))
  for (i in seq_along(base_names)) {
    suffix_id <- 0
    repeat {
      suffix <- if (suffix_id == 0) "" else paste0("_", suffix_id)
      candidate <- paste0(
        substr(base_names[i], 1, 31 - nchar(suffix)),
        suffix
      )
      previous <- if (i == 1) character() else result[seq_len(i - 1)]
      if (!tolower(candidate) %in% tolower(previous)) {
        result[i] <- candidate
        break
      }
      suffix_id <- suffix_id + 1
    }
  }
  result
}

make_unique_file_stems <- function(x) {
  base_names <- gsub("[/\\\\]", "_", as.character(x))
  base_names <- gsub("[<>:\"|?*]", "_", base_names)
  base_names <- gsub("[[:cntrl:]]", "_", base_names)
  base_names <- gsub("[. ]+$", "_", base_names)
  base_names[!nzchar(base_names) | base_names %in% c(".", "..")] <- "celltype"

  result <- character(length(base_names))
  for (i in seq_along(base_names)) {
    suffix_id <- 0
    repeat {
      suffix <- if (suffix_id == 0) "" else paste0("_", suffix_id)
      candidate <- paste0(base_names[i], suffix)
      previous <- if (i == 1) character() else result[seq_len(i - 1)]
      if (!tolower(candidate) %in% tolower(previous)) {
        result[i] <- candidate
        break
      }
      suffix_id <- suffix_id + 1
    }
  }
  result
}

## ---- 根据输入文件名推导输出前缀 ----
input_stem <- tools::file_path_sans_ext(basename(rds_file))
celltype_col_tag <- make_unique_file_stems(celltype_col)[[1]]
output_dir <- file.path("output", "presto")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
output_prefix <- file.path(
  output_dir,
  paste0(input_stem, "_", celltype_col_tag, "_presto")
)

message("Input RDS  : ", rds_file)
message("Input GTF  : ", gtf_file)
message("Celltype col: ", celltype_col)
message("Output prefix: ", output_prefix)

## ---- 读入 Seurat 对象 ----
message("Loading Seurat object ...")
obj <- readRDS(rds_file)
message("  Cells: ", ncol(obj), ", Features: ", nrow(obj))

if (!inherits(obj, "Seurat")) {
  stop("输入RDS不是Seurat对象")
}
if (!"RNA" %in% Assays(obj)) {
  stop("Seurat对象中缺少RNA assay")
}
if (!celltype_col %in% colnames(obj@meta.data)) {
  stop("meta.data 中未找到 ", celltype_col, " 列，请检查输入文件")
}

celltype_values <- as.character(obj[[celltype_col]][, 1])
if (anyNA(celltype_values) || any(!nzchar(celltype_values))) {
  stop("细胞类型metadata列包含NA或空值: ", celltype_col)
}
ct_levels <- unique(celltype_values)
message("  Cell types (", length(ct_levels), "): ", paste(sort(ct_levels), collapse = ", "))

## ---- 读入 GTF，筛选 protein-coding genes ----
message("Loading GTF reference ...")
gtf <- fread(gtf_file, data.table = FALSE)
required_gtf_cols <- c("gene_name", "gene_biotype")
missing_gtf_cols <- setdiff(required_gtf_cols, colnames(gtf))
if (length(missing_gtf_cols) > 0) {
  stop("注释CSV缺少必需列: ", paste(missing_gtf_cols, collapse = ", "))
}
pc_genes <- unique(trimws(as.character(
  gtf$gene_name[gtf$gene_biotype == "protein_coding"]
)))
pc_genes <- pc_genes[!is.na(pc_genes) & nzchar(pc_genes)]
if (length(pc_genes) == 0) {
  stop("注释CSV中未找到protein-coding gene_name")
}
message("  Protein-coding genes: ", length(pc_genes))

## ---- 提取表达矩阵和分组 (Seurat v5 兼容方式) ----
## 多个 RNA data layers 仅在内存中的对象副本上 join，不覆盖输入 RDS。
## 直接提取矩阵后调用 wilcoxauc.default，避免 presto v1.0.0 的 Seurat 接口兼容问题。
message("Extracting expression matrix (RNA data layer) ...")
## 注意：Layers() 的 search 参数在 v4 Assay 类下为精确匹配（非正则），
##       因此这里自行 grep，兼容 Assay 与 Assay5 两种格式。
data_layers <- grep("^data($|\\.)", Layers(obj, assay = "RNA"), value = TRUE)
if (length(data_layers) == 0) {
  stop("RNA assay中未找到data或data.* layer")
}
if (length(data_layers) > 1) {
  message("  Joining ", length(data_layers), " RNA data layers in memory ...")
  obj <- JoinLayers(obj, assay = "RNA", layers = "data", new = "data")
  data_layer <- "data"
} else {
  data_layer <- data_layers[[1]]
}
X_matrix <- LayerData(obj, assay = "RNA", layer = data_layer)
y <- FetchData(obj, celltype_col, cells = colnames(X_matrix)) %>%
  unlist() %>%
  as.character()
message("  Matrix: ", nrow(X_matrix), " genes x ", ncol(X_matrix), " cells")
rm(obj, data_layers, data_layer)
invisible(gc())

## ---- 运行 presto wilcoxauc ----
## 每个细胞类型 vs 其余所有细胞
message("Running presto wilcoxauc (one-vs-rest) ...")
degs <- wilcoxauc(X_matrix, y)
message("  Total rows: ", nrow(degs))
rm(y)

## ---- 匹配 protein-coding genes ----
message("Filtering to protein-coding genes ...")
degs_pc <- degs[degs$feature %in% pc_genes, ]
pc_universe <- base::intersect(rownames(X_matrix), pc_genes)
message("  Protein-coding rows: ", nrow(degs_pc))
message("  Tested protein-coding genes: ", length(pc_universe))
if (length(pc_universe) == 0) {
  stop("RNA feature与注释CSV的protein-coding gene_name没有交集，请确认使用SYMBOL")
}
rm(degs, X_matrix)
invisible(gc())

## ---- 按细胞类型拆分 ----
degs_list <- split(degs_pc, degs_pc$group)
if (length(degs_list) == 0) {
  stop("没有可输出的protein-coding差异分析结果")
}
rm(degs_pc)
invisible(gc())
sheet_names <- setNames(
  make_excel_sheet_names(names(degs_list)),
  names(degs_list)
)
gsea_file_stems <- setNames(
  make_unique_file_stems(names(degs_list)),
  names(degs_list)
)

for (ct in names(degs_list)) {
  n_deg <- nrow(degs_list[[ct]])
  n_sig <- sum(degs_list[[ct]]$padj < 0.05, na.rm = TRUE)
  message(sprintf("  %s: %d genes, %d padj<0.05", ct, n_deg, n_sig))
}

## ---- 保存完整 DEG 结果（每个细胞类型一个 sheet）----
message("Saving full DEG results to Excel ...")
full_wb <- createWorkbook()
for (ct in names(degs_list)) {
  df <- degs_list[[ct]] %>%
    arrange(desc(logFC), desc(auc), padj, feature)
  addWorksheet(full_wb, sheet_names[[ct]])
  writeData(full_wb, sheet_names[[ct]], df)
  freezePane(full_wb, sheet_names[[ct]], firstRow = TRUE)
}
full_xlsx <- paste0(output_prefix, "_full_DEG.xlsx")
saveWorkbook(full_wb, full_xlsx, overwrite = TRUE)
rm(full_wb, df)
invisible(gc())
message("Saved: ", full_xlsx)

## ---- 保存 ORA 背景基因 ----
universe_file <- paste0(output_prefix, "_ORA_universe.csv")
fwrite(data.frame(gene = pc_universe), universe_file)
message("Saved: ", universe_file)
rm(gtf, pc_genes, pc_universe)
invisible(gc())

## ---- 准备 ORA 输入：筛选可靠正向 marker 后按 logFC 取最多 100 个 ----
message(
  "Preparing ORA input (padj < ", marker_padj_max,
  ", logFC >= ", marker_logfc_min,
  ", auc >= ", marker_auc_min,
  ", pct_in >= ", marker_pct_in_min,
  "; top ", marker_n_max, " by logFC) ..."
)
top100_list <- lapply(degs_list, function(df) {
  df %>%
    filter(
      padj < marker_padj_max,
      logFC >= marker_logfc_min,
      auc >= marker_auc_min,
      pct_in >= marker_pct_in_min
    ) %>%
    arrange(desc(logFC), desc(auc), padj, feature) %>%
    slice_head(n = marker_n_max) %>%
    pull(feature)
})

## 检查基因数（ORA.r 要求每组合 >= 3 个基因才能做分析）
n_genes_per_ct <- sapply(top100_list, length)
message(
  "  Selected genes per cell type: ",
  paste(names(n_genes_per_ct), n_genes_per_ct, sep = "=", collapse = ", ")
)

wb <- createWorkbook()
for (ct in names(top100_list)) {
  genes <- top100_list[[ct]]
  addWorksheet(wb, sheet_names[[ct]])
  df_out <- data.frame(genes, stringsAsFactors = FALSE)
  colnames(df_out) <- ct
  writeData(wb, sheet_names[[ct]], df_out)
}
ora_xlsx <- paste0(output_prefix, "_ORA_top100.xlsx")
saveWorkbook(wb, ora_xlsx, overwrite = TRUE)
message("Saved: ", ora_xlsx)
rm(wb, top100_list)
invisible(gc())

## ---- 准备 GSEA 输入：每个细胞类型全基因按 logFC 降序排列 ----
## gsea.r 期望的输入格式：CSV，第一列 gene，第二列连续排序值（已降序）
message("Preparing GSEA inputs (full gene list ranked by logFC) ...")
gsea_files <- character(length(degs_list))
names(gsea_files) <- names(degs_list)

for (ct in names(degs_list)) {
  df <- degs_list[[ct]]
  gsea_df <- df[, c("feature", "logFC")]
  colnames(gsea_df) <- c("gene", "logFC")
  gsea_df <- gsea_df[
    !is.na(gsea_df$gene) &
      !duplicated(gsea_df$gene) &
      is.finite(gsea_df$logFC),
  ]
  ## 按 logFC 降序排列（gsea.r 内会再排一次，这里先排好）
  gsea_df <- gsea_df[order(gsea_df$logFC, decreasing = TRUE), ]
  ct_safe <- gsea_file_stems[[ct]]
  gsea_files[[ct]] <- paste0(output_prefix, "_GSEA_", ct_safe, ".csv")
  fwrite(gsea_df, gsea_files[[ct]])
  message("  ", basename(gsea_files[[ct]]), ": ", nrow(gsea_df), " genes")
}
message("Saved GSEA inputs under: ", output_dir)

message("===== Done =====")
message("Full DEG:    ", full_xlsx)
message("ORA input:   ", ora_xlsx)
message("ORA universe:", universe_file)
message("GSEA inputs: ", output_dir)
message("")
message("Next steps:")
message(
  "  ORA:  Rscript ORA.r ",
  shQuote(ora_xlsx), " SYMBOL ", shQuote(universe_file),
  " --gmt-dir=<GMT目录>"
)
for (ct in sort(names(degs_list))) {
  ct_safe <- gsea_file_stems[[ct]]
  message(
    "  GSEA: Rscript GSEA.r ",
    shQuote(gsea_files[[ct]]),
    " <GMT目录>"
  )
}
