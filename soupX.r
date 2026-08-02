## 输入：raw及filtered标准10X目录（matrix、features/genes、barcodes；基因顺序一致、filtered barcode为raw子集，并含HBB/IGKC/PTPRC/KLK3/EPCAM）；可选样本名。
## 输出：output/soupX/下以样本名开头的污染率CSV、校正对比PDF/PNG及校正后Seurat RDS。
## 示例：Rscript soupX.r data/sample01/outs/raw_feature_bc_matrix data/sample01/outs/filtered_feature_bc_matrix sample01

library(Seurat)
library(SoupX)
library(tidyverse)


## args 参数
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2 || length(args) > 3) {
    stop(
        "Usage: Rscript soupX.r <raw_10x_dir> <filtered_10x_dir> ",
        "[sample_name]"
    )
}

raw_matrix_path <- args[1]
filter_matrix_path <- args[2]
samples <- if (length(args) >= 3) {
    args[3]
} else {
    basename(dirname(dirname(normalizePath(filter_matrix_path, mustWork = TRUE))))
}
output_dir <- file.path("output", "soupX")

if (!dir.exists(raw_matrix_path)) {
    stop("Raw 10X directory does not exist: ", raw_matrix_path)
}
if (!dir.exists(filter_matrix_path)) {
    stop("Filtered 10X directory does not exist: ", filter_matrix_path)
}
if (!nzchar(samples)) {
    stop("sample_name must not be empty.")
}
cat("使用SoupX处理样本:", samples, "\n")

## 步骤1.1：读取原始矩阵和过滤后矩阵
raw_data <- Read10X(data.dir = raw_matrix_path)
filt_data <- Read10X(data.dir = filter_matrix_path)

## 步骤1.2：创建SoupChannel对象（核心对象）
sc <- SoupChannel(raw_data, filt_data)

## 步骤1.3：预处理过滤后数据以获取聚类信息（SoupX需要聚类辅助估计污染）
seob_temp <- CreateSeuratObject(counts = filt_data)
seob_temp <- NormalizeData(seob_temp)
seob_temp <- FindVariableFeatures(seob_temp, nfeatures = 2000)
seob_temp <- ScaleData(seob_temp)
seob_temp <- RunPCA(seob_temp)
seob_temp <- FindNeighbors(seob_temp, dims = 1:20)
seob_temp <- FindClusters(seob_temp, resolution = 0.8)
seob_temp <- RunUMAP(seob_temp, dims = 1:20)

## 步骤1.4：将聚类信息传入SoupChannel
clusters <- seob_temp@meta.data$seurat_clusters
names(clusters) <- colnames(seob_temp)
sc <- setClusters(sc, clusters)
sc <- setDR(sc, Embeddings(seob_temp, reduction = "umap"))

## 步骤1.5：自动估计污染比例
sc <- autoEstCont(sc)   
global_rho <- sc$fit$rhoEst # SoupX 1.6.2: 全局污染率标量
cat("估计全局污染率rho:", round(global_rho, 4), "\n")

## 保存污染率估计结果
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
if (!dir.exists(output_dir)) {
    stop("Failed to create output directory: ", output_dir)
}
contamination_csv <- file.path(output_dir, paste0(samples, "_contamination_estimate.csv"))
plot_pdf <- file.path(output_dir, paste0(samples, "_soupx_correction_plots.pdf"))
plot_png <- file.path(output_dir, paste0(samples, "_soupx_correction_plots.png"))
corrected_rds <- file.path(output_dir, paste0(samples, "_soupx_corrected.rds"))
contamination_result <- data.frame(
    sample = samples,
    global_rho = round(global_rho, 4),
    n_cells = ncol(sc$toc),
    n_genes = nrow(sc$toc),
    interpretation = dplyr::case_when(
        global_rho < 0.02 ~ "可忽略(<2%)",
        global_rho < 0.05 ~ "轻度(2-5%)",
        global_rho < 0.10 ~ "建议校正(5-10%)",
        TRUE ~ "强烈建议校正(>10%)"
    )
)
write.csv(contamination_result,
    file = contamination_csv,
    row.names = FALSE
)
cat("污染率估计结果已保存\n")

## 步骤1.6：校正表达矩阵（去除背景污染）
adj_data <- adjustCounts(sc, roundToInt = TRUE) # 输出校正后的count矩阵

## 步骤1.7：用校正后的矩阵创建Seurat对象
seob <- CreateSeuratObject(
    counts = adj_data
)

## 步骤1.8：对比校正前后表达矩阵的差异
title_theme <- theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"))

p1 <- plotChangeMap(sc, adj_data, "HBB") + ggtitle("HBB") + title_theme
p2 <- plotChangeMap(sc, adj_data, "IGKC") + ggtitle("IGKC") + title_theme
p3 <- plotChangeMap(sc, adj_data, "PTPRC") + ggtitle("PTPRC") + title_theme
p4 <- plotChangeMap(sc, adj_data, "KLK3") + ggtitle("KLK3") + title_theme
p5 <- plotChangeMap(sc, adj_data, "EPCAM") + ggtitle("EPCAM") + title_theme
p <- cowplot::plot_grid(p1, p2, p3, p4, p5, ncol = 3)
ggsave(filename = plot_pdf, plot = p, width = 12, height = 8)
ggsave(filename = plot_png, plot = p, width = 12, height = 8)
cat("校正前后表达差异图已保存为:", plot_pdf, "和", plot_png, "\n")

## 步骤1.9：保存校正后的Seurat对象
saveRDS(seob, file = corrected_rds)
cat("SoupX处理完成，校正后的对象已保存为:", corrected_rds, "\n")
