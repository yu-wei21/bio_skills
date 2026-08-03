## 输入：
##   1. gene×sample 非负整数 count 矩阵（首列为唯一 SYMBOL）
##   2. 分组表（前两列为文本型唯一样本 ID 和分组，且恰好两组）
##   3. 分子组名称
##   4. 分母组名称
##   5. 人类 GTF 注释 CSV（含 gene_name、gene_biotype 列）
##   6. 下游分析使用的 GMT 目录
## 样本筛选：取 count 矩阵列与分组表样本的交集（两侧独有样本均忽略并提示）。
##   即分组表可包含 count 矩阵中不存在的样本（如 sc2pseudobulk 因该细胞类型无细胞
##   而剔除的样本），count 矩阵也可包含分组表中无临床信息的样本；分析只用交集。
## 输出：
##   1. output/DESeq2/ 下以 count 矩阵文件名和比较名称为前缀的 DESeq2 结果
##   2. PCA 与 VST 图形 PDF/PNG
##   3. ORA 背景与基因集输入文件
##   4. GSEA 排序输入文件
## 示例命令：
##   Rscript DESeq2.r counts.csv group.csv treated control genes.csv data/gmt

suppressPackageStartupMessages({
  library(DESeq2)
  library(ggplot2)
  library(cowplot)
  library(readr)
  library(openxlsx)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 6) {
  stop(paste0(
    "Usage: Rscript DESeq2.r <count_file> <group_file> ",
    "<numerator_group> <denominator_group> <gtf_file> <gmt_dir>"
  ))
}

count_file <- args[1]
group_file <- args[2]
numerator_group <- args[3]
denominator_group <- args[4]
gtf_file <- args[5]
gmt_dir <- args[6]
output_dir <- file.path("output", "DESeq2")

ora_padj_max <- 0.05
ora_log2fc_min <- 0
ora_n_max <- 100


dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
count_tag <- tools::file_path_sans_ext(basename(count_file))
contrast_name <- paste0(numerator_group, "_vs_", denominator_group)
contrast_tag <- gsub("[\\\\/:*?\"<>|]", "_", contrast_name)
contrast_tag <- gsub("[[:space:]]+", "_", contrast_tag)
suffix <- paste0(count_tag, "_", contrast_tag)

read_table_auto <- function(file, all_character = FALSE) {
  ext <- tolower(tools::file_ext(file))
  if (ext %in% c("xlsx", "xls")) {
    x <- openxlsx::read.xlsx(file, check.names = FALSE)
    if (all_character) {
      x[] <- lapply(x, as.character)
    }
    return(x)
  }
  if (ext %in% c("csv")) {
    return(read.csv(
      file,
      check.names = FALSE,
      stringsAsFactors = FALSE,
      colClasses = if (all_character) "character" else NA
    ))
  }
  if (ext %in% c("tsv", "txt")) {
    return(read.delim(
      file,
      check.names = FALSE,
      stringsAsFactors = FALSE,
      colClasses = if (all_character) "character" else NA
    ))
  }
  read.delim(
    file,
    check.names = FALSE,
    stringsAsFactors = FALSE,
    colClasses = if (all_character) "character" else NA
  )
}

message("[1/4] 读取输入文件...")
if (!file.exists(count_file)) {
  stop("count 矩阵文件不存在: ", count_file)
}
if (!file.exists(group_file)) {
  stop("分组文件不存在: ", group_file)
}
if (!file.exists(gtf_file)) {
  stop("GTF 注释文件不存在: ", gtf_file)
}
if (!dir.exists(gmt_dir)) {
  stop("GMT 目录不存在: ", gmt_dir)
}
if (!nzchar(numerator_group) || !nzchar(denominator_group)) {
  stop("分子组和分母组名称不能为空")
}
if (identical(numerator_group, denominator_group)) {
  stop("分子组和分母组不能相同: ", numerator_group)
}

count_df <- read_table_auto(count_file)
sample <- read_table_auto(group_file, all_character = TRUE)
gtf <- read_table_auto(gtf_file)

if (ncol(sample) < 2) {
  stop("分组文件至少需要2列: 样本ID列 + 分组列")
}
if (ncol(count_df) < 2) {
  stop("count矩阵至少需要2列: 基因ID列 + 样本列")
}
required_gtf_cols <- c("gene_name", "gene_biotype")
missing_gtf_cols <- setdiff(required_gtf_cols, colnames(gtf))
if (length(missing_gtf_cols) > 0) {
  stop("GTF 注释缺少必需列: ", paste(missing_gtf_cols, collapse = ", "))
}

sample <- sample[, 1:2, drop = FALSE]
colnames(sample) <- c("sample_id", "group")
sample$sample_id <- as.character(sample$sample_id)
sample$group <- as.character(sample$group)
invalid_sample <- is.na(sample$sample_id) | !nzchar(sample$sample_id) |
  trimws(sample$sample_id) != sample$sample_id
invalid_group <- is.na(sample$group) | !nzchar(sample$group) |
  trimws(sample$group) != sample$group
if (any(invalid_sample) || any(invalid_group)) {
  stop("分组表的样本ID和分组不得为NA、空值或带首尾空格")
}

if (anyDuplicated(sample$sample_id)) {
  stop("分组文件存在重复样本ID: ", paste(unique(sample$sample_id[duplicated(sample$sample_id)]), collapse = ", "))
}

observed_groups <- unique(sample$group)
if (length(observed_groups) != 2) {
  stop("当前脚本要求分组文件恰好包含2个组，当前为: ", paste(observed_groups, collapse = ", "))
}
if (!setequal(observed_groups, c(numerator_group, denominator_group))) {
  stop(
    "指定的分子/分母组与分组文件不一致。分组文件: ",
    paste(observed_groups, collapse = ", "),
    "; 指定比较: ", numerator_group, " vs ", denominator_group
  )
}

gene_ids <- as.character(count_df[[1]])
if (any(is.na(gene_ids) | gene_ids == "")) {
  stop("count 矩阵的基因 ID 存在 NA 或空值")
}
if (anyDuplicated(gene_ids)) {
  stop("count 矩阵存在重复基因 SYMBOL: ", paste(head(unique(gene_ids[duplicated(gene_ids)]), 10), collapse = ", "))
}
ensembl_ids <- gene_ids[grepl("^ENSG[0-9]+(?:\\.[0-9]+)?$", gene_ids, perl = TRUE)]
if (length(ensembl_ids) > 0) {
  stop(
    "count 矩阵必须使用标准 gene SYMBOL，不接受 Ensembl ID。示例: ",
    paste(head(ensembl_ids, 10), collapse = ", ")
  )
}
entrez_ids <- gene_ids[grepl("^[0-9]+$", gene_ids)]
if (length(entrez_ids) > 0) {
  stop(
    "count 矩阵必须使用标准 gene SYMBOL，不接受纯数字 ENTREZID。示例: ",
    paste(head(entrez_ids, 10), collapse = ", ")
  )
}

gtf$gene_name <- as.character(gtf$gene_name)
gtf$gene_biotype <- as.character(gtf$gene_biotype)
pc_genes <- unique(gtf$gene_name[
  !is.na(gtf$gene_name) &
    gtf$gene_name != "" &
    gtf$gene_biotype == "protein_coding"
])
if (length(pc_genes) == 0) {
  stop("GTF 注释中未找到 gene_biotype == protein_coding 的 gene_name")
}
all_gtf_symbols <- unique(gtf$gene_name[!is.na(gtf$gene_name) & gtf$gene_name != ""])
unmatched_gtf <- setdiff(gene_ids, all_gtf_symbols)
n_pc_input <- sum(gene_ids %in% pc_genes)
message("GTF protein-coding SYMBOL 数: ", length(pc_genes))
message("count 矩阵中 protein-coding SYMBOL 数: ", n_pc_input, " / ", length(gene_ids))
if (n_pc_input == 0) {
  stop("count 矩阵中没有任何基因能够匹配 GTF protein-coding SYMBOL")
}
if (length(unmatched_gtf) > 0) {
  message(
    "警告: count 矩阵中有 ", length(unmatched_gtf),
    " 个 SYMBOL 无法匹配当前 GTF，将保留在完整 DESeq2 结果中但不进入富集输出。示例: ",
    paste(head(unmatched_gtf, 10), collapse = ", ")
  )
}

rownames(count_df) <- gene_ids
count_df <- count_df[, -1, drop = FALSE]
if (anyDuplicated(colnames(count_df))) {
  stop(
    "count矩阵存在重复样本列名: ",
    paste(unique(colnames(count_df)[duplicated(colnames(count_df))]), collapse = ", ")
  )
}
## 样本筛选：取交集（分组表可含被 sc2pseudobulk 剔除的样本；count 矩阵可含无临床信息的样本）
group_only <- setdiff(sample$sample_id, colnames(count_df))
count_only <- setdiff(colnames(count_df), sample$sample_id)
if (length(group_only) > 0) {
  message(
    "忽略分组表独有（count矩阵中不存在，可能因该细胞类型无细胞被 sc2pseudobulk 剔除）样本: ",
    paste(group_only, collapse = ", ")
  )
}
if (length(count_only) > 0) {
  message(
    "忽略count矩阵独有（分组表无临床分组信息）样本: ",
    paste(count_only, collapse = ", ")
  )
}
shared_samples <- intersect(colnames(count_df), sample$sample_id)
if (length(shared_samples) == 0) {
  stop("count矩阵与分组表没有共同样本，无法进行分析。")
}

count_df <- count_df[, shared_samples, drop = FALSE]
sample <- sample[match(shared_samples, sample$sample_id), , drop = FALSE]
rownames(sample) <- sample$sample_id
count_df <- as.matrix(count_df)
mode(count_df) <- "numeric"
if (anyNA(count_df) || any(!is.finite(count_df))) {
  stop("count矩阵包含NA或非有限数值，可能存在非数值单元格")
}
if (any(count_df < 0)) {
  stop("count矩阵包含负值")
}
if (any(abs(count_df - round(count_df)) > sqrt(.Machine$double.eps))) {
  stop("count矩阵必须为整数原始counts")
}
zero_samples <- colnames(count_df)[colSums(count_df) == 0]
if (length(zero_samples) > 0) {
  stop(
    "以下样本总counts为0，无法估计DESeq2 size factor: ",
    paste(zero_samples, collapse = ", "),
    "。请核查该样本在对应细胞大类中是否没有细胞。"
  )
}
if (!identical(rownames(sample), colnames(count_df))) {
  stop("内部错误: colData行名与count矩阵列名未严格对齐")
}
sample$group <- factor(sample$group, levels = c(denominator_group, numerator_group))
group_counts <- table(sample$group)
if (any(group_counts < 2)) {
  stop(
    "每组至少需要2个生物学重复；当前为: ",
    paste(names(group_counts), as.integer(group_counts), sep = "=", collapse = ", ")
  )
}

message("纳入样本数: ", nrow(sample))
message("分组计数: ", paste(names(table(sample$group)), as.integer(table(sample$group)), sep = "=", collapse = ", "))
message("样本顺序: ", paste(sample$sample_id, collapse = ", "))
message("差异比较: ", numerator_group, " vs ", denominator_group)

message("[2/4] 构建DESeq2对象并过滤低表达基因...")
dds <- DESeqDataSetFromMatrix(
  countData = count_df,
  colData = sample,
  design = ~group
)
n_genes_before_filter <- nrow(dds)
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep, ]
message(
  "低表达过滤: ", n_genes_before_filter,
  " -> ", nrow(dds), " genes (rowSums(counts) >= 10)"
)
if (nrow(dds) == 0) {
  stop("低表达过滤后没有剩余基因，无法继续分析")
}

message("[3/4] PCA分析并输出图和坐标...")
if (nrow(dds) < 1000) {
  vsd <- varianceStabilizingTransformation(dds, blind = TRUE)
} else {
  ## vst() 需要足够多的 mean normalized count > 5 的基因，否则报错
  ## （如稀有细胞类型的低表达过滤后基因集）；回退到 VST（官方建议，
  ## 仅影响 PCA/VST 输出，不影响 DESeq2 差异检验结果）。
  vsd <- tryCatch(
    vst(dds, blind = TRUE),
    error = function(e) {
      message(
        "vst() failed (", conditionMessage(e),
        "); falling back to varianceStabilizingTransformation."
      )
      varianceStabilizingTransformation(dds, blind = TRUE)
    }
  )
}
vsd_data <- assay(vsd)
pca_data <- plotPCA(vsd, intgroup = "group", returnData = TRUE)
pca_data$sample_id <- rownames(pca_data)
percent_var <- round(100 * attr(pca_data, "percentVar"))

p <- ggplot(pca_data, aes(PC1, PC2, color = group, label = sample_id)) +
  geom_point(size = 3) +
  geom_text(aes(label = sample_id),
    nudge_y = 2, nudge_x = 1
  ) +
  xlab(paste0("PC1: ", percent_var[1], "% variance")) +
  ylab(paste0("PC2: ", percent_var[2], "% variance")) +
  coord_fixed() +
  theme_cowplot()

pca_plot_pdf_file <- file.path(output_dir, paste0(suffix, "_PCA_plot.pdf"))
pca_plot_png_file <- file.path(output_dir, paste0(suffix, "_PCA_plot.png"))
pca_coord_file <- file.path(output_dir, paste0(suffix, "_PCA_coordinates.csv"))
deg_file <- file.path(output_dir, paste0(suffix, "_DESeq2_results.csv"))
vsd_file <- file.path(output_dir, paste0(suffix, "_DESeq2_vsd.csv"))
ora_file <- file.path(output_dir, paste0(suffix, "_DESeq2_ORA_top100.xlsx"))
ora_universe_file <- file.path(output_dir, paste0(suffix, "_DESeq2_ORA_universe.csv"))
gsea_file <- file.path(output_dir, paste0(suffix, "_DESeq2_GSEA_ranked.csv"))

ggsave(pca_plot_pdf_file, p, width = 6, height = 6)
ggsave(pca_plot_png_file, p, width = 6, height = 6)
write.csv(pca_data, pca_coord_file, row.names = FALSE)
write.csv(vsd_data, vsd_file, row.names = TRUE)

message("[4/4] DESeq2差异分析并导出结果...")
dds <- DESeq(dds)
res <- results(dds, contrast = c("group", numerator_group, denominator_group))
message("contrast is: ", numerator_group, " vs ", denominator_group)

res_df <- as.data.frame(res)
res_df$gene <- rownames(res_df)
res_df <- res_df[, c("gene", setdiff(colnames(res_df), "gene"))]

write.csv(res_df, deg_file, row.names = FALSE)

## ORA 输入：显著上调基因，按 log2 fold change 排序。
ora_df <- res_df[
  !is.na(res_df$gene) &
    res_df$gene != "" &
    res_df$gene %in% pc_genes &
    is.finite(res_df$pvalue) &
    is.finite(res_df$log2FoldChange) &
    is.finite(res_df$padj) &
    res_df$padj < ora_padj_max &
    res_df$log2FoldChange > ora_log2fc_min,
  c("gene", "log2FoldChange", "padj"),
  drop = FALSE
]
ora_df <- ora_df[order(-ora_df$log2FoldChange, ora_df$padj, ora_df$gene), , drop = FALSE]
ora_df <- ora_df[!duplicated(ora_df$gene), , drop = FALSE]
ora_df <- head(ora_df, ora_n_max)

ora_sheet <- gsub("[\\\\/:*?\\[\\]]", "_", contrast_name)
ora_sheet <- substr(ora_sheet, 1, 31)
ora_input <- data.frame(ora_df$gene, stringsAsFactors = FALSE)
colnames(ora_input) <- contrast_name

ora_wb <- createWorkbook()
addWorksheet(ora_wb, ora_sheet)
writeData(ora_wb, ora_sheet, ora_input)
saveWorkbook(ora_wb, ora_file, overwrite = TRUE)

## ORA 背景：已检验、具有有限原始 p 值的蛋白编码 SYMBOL。
ora_universe <- unique(as.character(
  res_df$gene[
    !is.na(res_df$gene) &
      res_df$gene != "" &
      res_df$gene %in% pc_genes &
      is.finite(res_df$pvalue)
  ]
))
write.csv(data.frame(gene = ora_universe), ora_universe_file, row.names = FALSE)

## GSEA 输入：已检验的蛋白编码 SYMBOL，按未收缩的 log2 fold change 排序。
gsea_df <- res_df[
  !is.na(res_df$gene) &
    res_df$gene != "" &
    res_df$gene %in% pc_genes &
    is.finite(res_df$pvalue) &
    is.finite(res_df$log2FoldChange),
  c("gene", "log2FoldChange"),
  drop = FALSE
]
gsea_df <- gsea_df[!duplicated(gsea_df$gene), , drop = FALSE]
gsea_df <- gsea_df[order(gsea_df$log2FoldChange, decreasing = TRUE), , drop = FALSE]
write.csv(gsea_df, gsea_file, row.names = FALSE)

message(
  "ORA genes selected: ", nrow(ora_df),
  " (padj < ", ora_padj_max,
  ", log2FoldChange > ", ora_log2fc_min,
  "; top ", ora_n_max, ")"
)
if (nrow(ora_df) == 0) {
  message("警告: 没有符合 ORA 筛选条件的基因，已输出只含列标题的 Excel。")
} else if (nrow(ora_df) < 3) {
  message("警告: ORA 候选基因少于3个，ORA.r 将跳过该基因集。")
}
message("ORA universe genes: ", length(ora_universe))
message("GSEA ranked genes: ", nrow(gsea_df))

message("完成。输出文件:")
message(" - ", pca_plot_pdf_file)
message(" - ", pca_plot_png_file)
message(" - ", pca_coord_file)
message(" - ", vsd_file)
message(" - ", deg_file)
message(" - ", ora_file)
message(" - ", ora_universe_file)
message(" - ", gsea_file)
message("")
message("Next steps:")
message(
  "  ORA:  Rscript ORA.r ",
  shQuote(ora_file), " SYMBOL ", shQuote(ora_universe_file),
  " --gmt-dir=", shQuote(gmt_dir)
)
message(
  "  GSEA: Rscript GSEA.r ",
  shQuote(gsea_file), " ", shQuote(gmt_dir)
)
