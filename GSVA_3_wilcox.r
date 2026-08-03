## 输入：GSVA评分CSV（首列唯一且可作R列名的通路名，其余列为样本数值评分）和临床CSV（IDbn唯一且完整匹配样本列；FDG_groupe取0/2且两组均存在）。
## 输出：output/GSVA_3_wilcox/下以GSVA评分文件名开头的Wilcoxon差异结果CSV。
## 示例：Rscript GSVA_3_wilcox.r output/GSVA/input_expression_TPM_pathway_GSVA_result.csv data/clinical.csv

## 脚本目的： 使用Wilcoxon秩和检验进行差异分析
library(tidyverse)
library(rstatix)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) {
    stop("Usage: Rscript GSVA_3_wilcox.r <gsva_result.csv> <clinical.csv>")
}
gsva_file <- args[1]
clinical_file <- args[2]
output_dir <- file.path("output", "GSVA_3_wilcox")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
input_tag <- tools::file_path_sans_ext(basename(gsva_file))
output_file <- file.path(output_dir, paste0(input_tag, "_wilcox_result.csv"))

## 读取GSVA结果
## check.names=FALSE 保留原始样本 ID，确保与 clinical CSV 的 IDbn 精确匹配
gsva_mat <- read.csv(gsva_file, row.names = 1, check.names = FALSE)
gsva_mat <- gsva_mat %>%
    t() %>%
    as.data.frame() %>%
    rownames_to_column(var = "IDbn")
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
data <- clinical %>%
    dplyr::select(any_of(c("IDbn", "FDG_groupe"))) %>%
    left_join(gsva_mat, by = "IDbn")
data %>% head()

## 逐个通路进行差异分析
colnames(data) <- gsub("\\/", "", colnames(data))
colnames(data) <- gsub("\\,", "", colnames(data))
colnames(data) <- gsub("\\(", "", colnames(data))
colnames(data) <- gsub("\\)", "", colnames(data))
var <- colnames(data)[-c(1, 2)]
var %>% head()

res <- lapply(var, function(x) {
    print(x)
    res <- data %>%
        wilcox_test(as.formula(paste(x, "~ FDG_groupe")), detailed = TRUE)
    res <- res %>%
        mutate(pathway = x)
    return(res)
})
res <- do.call(rbind, res)
res %>% head()

## 导出结果
write.csv(res, file = output_file, row.names = FALSE)
