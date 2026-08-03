## 输入：
##   1. gene×sample 非负数值表达矩阵 CSV（首列为唯一 gene symbol、与 GMT 一致且无 NA）
##   2. 表达类型：TPM、RPKM、CPM 或 count
##   3. 标准 GMT 文件
## 输出：
##   1. output/GSVA/ 下以表达矩阵文件名为前缀的 GSVA 评分 CSV
## 示例命令：Rscript GSVA.r input_expression.csv TPM pathway.gmt

## GSVA分析脚本

library(GSVA)
library(tidyverse)
library(clusterProfiler)

## 参数解析
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3) {
    stop("Usage: Rscript GSVA.r <expression.csv> <TPM|RPKM|CPM|count> <pathway.gmt>")
}
input_rna <- args[1]
input_style <- args[2]
input_gmt <- args[3]
output_dir <- file.path("output", "GSVA")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
input_tag <- tools::file_path_sans_ext(basename(input_rna))
gmt_tag <- tools::file_path_sans_ext(basename(input_gmt))
output_gsva <- file.path(
    output_dir,
    paste0(input_tag, "_", input_style, "_", gmt_tag, "_GSVA_result.csv")
)

## 读取表达矩阵文件
## check.names=FALSE 保留原始样本 ID（如含 - 的样本名，避免被改写成 QXH.G 导致下游匹配失败）
matrix <- read.csv(input_rna, row.names = 1, check.names = FALSE) %>% as.matrix()
matrix %>% head()

## 读取 GMT 文件
pathway <- read.gmt(input_gmt)
pathway %>% head()
pathway <- split(pathway$gene, pathway$term)

## GSVA分析
## GSVA >= 2.0 使用 gsvaParam() 参数对象（旧签名 gsva(expr=, gset.idx.list=) 已 defunct）
if (input_style == "TPM" | input_style == "RPKM" | input_style == "CPM") {
    matrix <- log2(matrix + 1)
    matrix %>% head()
    gsva_mat <- gsva(
        gsvaParam(
            exprData = matrix,
            geneSets = pathway,
            kcdf = "Gaussian" # "Gaussian" for logCPM,logRPKM,logTPM, "Poisson" for counts
        )
    )
} else if (input_style == "count") {
    gsva_mat <- gsva(
        gsvaParam(
            exprData = matrix,
            geneSets = pathway,
            kcdf = "Poisson" # "Gaussian" for logCPM,logRPKM,logTPM, "Poisson" for counts
        )
    )
}
gsva_mat %>% head()

## 导出结果
write.csv(gsva_mat, file = output_gsva, row.names = TRUE)

## 脚本结束
cat("GSVA analysis completed successfully.\n")
