## 输入：带表头的CSV/TSV/TXT基因排序表（前两列为与GMT一致的gene symbol和有限数值排序值）及含标准GMT文件的基因集目录。
## 输出：output/GSEA/下以排序表文件名开头的GSEA结果XLSX及图形PDF/PNG。
## 示例：Rscript GSEA.r data/gene_list.csv data/gmt

## 脚本目的：对基因表达数据进行GSEA富集分析
## GMT基因集使用gene symbol。

suppressPackageStartupMessages({
    library(clusterProfiler)
    library(openxlsx)
    library(cowplot)
    library(enrichplot)
    library(parallel)
    library(data.table)
    library(tidyverse)
})

## read files
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) {
    stop("Usage: Rscript GSEA.r <ranked_table> <gmt_dir>")
}
data_file <- args[1]
gmt_dir <- args[2]
if (!file.exists(data_file)) {
    stop("Gene ranking file does not exist: ", data_file)
}
output_dir <- file.path("output", "GSEA")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
input_tag <- tools::file_path_sans_ext(basename(data_file))
output_prefix <- file.path(output_dir, input_tag)
data <- fread(data_file, data.table = FALSE)

## prepare data
## 两列，第一列为基因名，第二列为排序值（大到小）
if (ncol(data) < 2) {
    stop("Gene ranking file must contain at least two columns: gene and ranking value.")
}
gene_ids <- trimws(as.character(data[[1]]))
rank_values <- suppressWarnings(as.numeric(data[[2]]))
if (anyNA(gene_ids) || any(!nzchar(gene_ids))) {
    stop("The first column contains NA or empty gene IDs.")
}
if (anyNA(rank_values) || any(!is.finite(rank_values))) {
    stop("The second column must contain only finite numeric ranking values.")
}
gene_list <- rank_values
names(gene_list) <- gene_ids
gene_list <- gene_list[!duplicated(names(gene_list))]
gene_list <- sort(gene_list, decreasing = TRUE)
if (length(gene_list) < 10) {
    stop("At least 10 unique ranked genes are required for GSEA.")
}

## read gmt file
if (!dir.exists(gmt_dir)) {
    stop("GMT directory does not exist: ", gmt_dir)
}
gmt_files <- list.files(gmt_dir, pattern = "\\.gmt$", full.names = TRUE)
if (length(gmt_files) == 0) {
    stop("No GMT files found in: ", gmt_dir)
}
gmt_list <- lapply(gmt_files, function(x) {
    gmt <- read.gmt(x)
    return(gmt)
})
names(gmt_list) <- gsub("\\.gmt$", "", basename(gmt_files))
gmt_genes <- unique(unlist(lapply(gmt_list, function(x) as.character(x$gene))))
n_overlap <- sum(names(gene_list) %in% gmt_genes)
message("Ranked genes overlapping GMT genes: ", n_overlap, " / ", length(gene_list))
if (n_overlap == 0) {
    stop("No ranked gene IDs overlap the GMT files; check that both use the same SYMBOL system.")
}

## GSEA analysis
set.seed(666)
gsea_results <- mclapply(gmt_list, function(gmt) {
    gsea_result <- GSEA(gene_list, TERM2GENE = gmt, pvalueCutoff = 0.05)
    return(gsea_result)
}, mc.cores = 4)
names(gsea_results) <- names(gmt_list)

## save results
gsea_tables <- lapply(gsea_results, as.data.frame)
write.xlsx(gsea_tables, file = paste0(output_prefix, "_gsea_results.xlsx"), overwrite = TRUE)

## 合并结果
gsea_results_nonempty <- Filter(function(x) nrow(as.data.frame(x)) > 0, gsea_results)
gsea_results_df <- do.call(rbind, lapply(names(gsea_results_nonempty), function(name) {
    df <- as.data.frame(gsea_results_nonempty[[name]])
    df$Category <- name
    return(df)
}))

## 可视化
if (is.null(gsea_results_df) || nrow(gsea_results_df) == 0) {
    message("No significant GSEA results; writing an empty-result plot.")
    p <- ggplot() +
        annotate("text", x = 0, y = 0, label = "No significant GSEA results") +
        xlim(-1, 1) +
        ylim(-1, 1) +
        theme_void()
} else {
    ## X轴标题改为Title Case
    gsea_results_df$Description <- gsub("_", " ", gsea_results_df$Description)
    gsea_results_df$Description <- str_to_title(gsea_results_df$Description)
    gsea_results_df$Description <- str_trunc(gsea_results_df$Description, width = 60, ellipsis = "...")

    p <- gsea_results_df %>%
        group_by(Category) %>%
        slice_max(order_by = abs(NES), n = 10) %>%
        ungroup() %>%
        ggplot(aes(x = reorder(Description, NES), y = NES, fill = Category)) +
        geom_bar(stat = "identity") +
        facet_grid(Category ~ ., scales = "free", space = "free") +
        coord_flip() +
        theme_cowplot() +
        labs(x = "Pathway", y = "Normalized Enrichment Score (NES)", title = "GSEA Enrichment Analysis") +
        theme(
            legend.position = "none",
            plot.margin = margin(5, 5, 5, 10, "mm")
        )
}
ggsave(paste0(output_prefix, "_gsea_plot.pdf"), plot = p, width = 10, height = 20)
ggsave(paste0(output_prefix, "_gsea_plot.png"), plot = p, width = 10, height = 20)
