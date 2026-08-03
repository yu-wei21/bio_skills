## 输入：
##   1. GSVA 评分 CSV（首列为唯一通路名，其余列为样本数值评分）
##   2. 临床 CSV（IDbn 唯一且完整匹配样本列；FDG_groupe 取 0/2 且两组均存在）
## 输出：
##   1. output/GSVA_2_limma/ 下以 GSVA 评分文件名为前缀的 limma 差异结果 CSV
## 示例命令：Rscript GSVA_2_limma.r output/GSVA/input_expression_TPM_pathway_GSVA_result.csv data/clinical.csv
## 注意事项：项目特异脚本；仅比较 FDG_groupe 为 0/2 的样本（映射为 Low/High）。

library(limma)
library(tidyverse)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) {
    stop("Usage: Rscript GSVA_2_limma.r <gsva_result.csv> <clinical.csv>")
}
gsva_file <- args[1]
clinical_file <- args[2]
output_dir <- file.path("output", "GSVA_2_limma")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
input_tag <- tools::file_path_sans_ext(basename(gsva_file))
output_file <- file.path(output_dir, paste0(input_tag, "_limma_result.csv"))

## 读取GSVA结果
## check.names=FALSE 保留原始样本 ID，确保与 clinical CSV 的 IDbn 精确匹配
gsva_mat <- read.csv(gsva_file, row.names = 1, check.names = FALSE)
gsva_mat %>% head()

## 读取临床信息
clinical <- read.csv(clinical_file)
clinical %>% head()

## 提取FDG_groupe中0,2样本
clinical <- clinical %>%
    filter(FDG_groupe %in% c("0", "2")) %>%
    mutate(FDG_groupe = case_when(
        FDG_groupe == "0" ~ "Low",
        FDG_groupe == "2" ~ "High"
    ))
clinical %>% head()

## 提取GSVA评分中的0,2样本
gsva_mat <- gsva_mat %>%
    dplyr::select(any_of(clinical$IDbn))
gsva_mat %>% head()

## 构建设计矩阵
factor <- factor(clinical$FDG_groupe)
design <- model.matrix(~ 0 + factor)
colnames(design) <- levels(factor)
rownames(design) <- clinical$IDbn
design

## 构建差异分析
# 定义对比组合：High组 vs Low组的差异
compare <- makeContrasts(High - Low, levels = design)
# 线性模型拟合：使用GSVA评分矩阵和设计矩阵进行拟合
fit <- lmFit(gsva_mat, design)
# 应用对比矩阵：计算High组与Low组之间的差异
fit <- contrasts.fit(fit, compare)
# 经验贝叶斯统计：计算统计量和p值
fit <- eBayes(fit)
# 提取差异分析结果：获取所有通路的统计结果
Diff <- topTable(fit, coef = 1, number = Inf)

## 导出结果
write.csv(Diff, file = output_file, row.names = TRUE)
