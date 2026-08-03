## 输入：
##   1. GSEA 结果 XLSX：每个 sheet 须含区分大小写的 ID、Description、NES、p.adjust 列，后两列为数值
## 输出：
##   1. output/plot-GSEA-dotplot/ 下以输入 XLSX 文件名为前缀的配对点图 PDF
##   2. output/plot-GSEA-dotplot/ 下以输入 XLSX 文件名为前缀的配对点图 PNG
## 示例命令：Rscript plot-GSEA-dotplot.r results/gsea_results.xlsx
## 注意事项：项目主题图；仅绘制脚本内置 FDG 相关通路，输入 ID 必须与内置列表完全匹配。

library(tidyverse)
library(cowplot)
library(openxlsx)

## 读取输入文件
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1) {
    stop("Usage: Rscript plot-GSEA-dotplot.r <gsea_results.xlsx>")
}
input_file <- args[1]
if (!file.exists(input_file)) {
    stop("Input XLSX does not exist: ", input_file)
}
output_dir <- file.path("output", "plot-GSEA-dotplot")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
input_tag <- tools::file_path_sans_ext(basename(input_file))
output_prefix <- file.path(output_dir, input_tag)
# 将各 sheet 读入为 data.frame 列表
sheets <- openxlsx::getSheetNames(input_file)
data_list <- lapply(sheets, function(s) read.xlsx(input_file, sheet = s))
data <- do.call(rbind, data_list) %>%
    as.data.frame()

## 筛选感兴趣的通路
pathway <- c(
    ## OXPHOS 抑制（FDG 核心机制）：NES 为负，在 FDG 阴性细胞中富集
    "HALLMARK_OXIDATIVE_PHOSPHORYLATION", # NES=-1.94, padj=2.9e-07
    "Oxidative phosphorylation", # NES=-2.00, padj=1.1e-05
    "Cellular Respiration (GO:0045333)", # NES=-2.05, padj=3.8e-04
    "KEGG_OXIDATIVE_PHOSPHORYLATION", # NES=-1.56, padj=1.5e-02
    "GOBP_OXIDATIVE_PHOSPHORYLATION", # NES=-1.51, padj=4.2e-02
    ## 缺氧 / HIF-1α 激活（增强 GLUT1/LDHA，稳定糖酵解）
    "HALLMARK_HYPOXIA", # NES=+1.43, padj=2.1e-03
    "GROSS_HYPOXIA_VIA_ELK3_DN", # NES=+1.69, padj=2.6e-07
    ## TNF-α / NF-κB 炎症信号（糖酵解的转录驱动因素）
    "HALLMARK_TNFA_SIGNALING_VIA_NFKB", # NES=+1.71, padj=5.0e-09
    "HALLMARK_INFLAMMATORY_RESPONSE", # NES=+1.30, padj=4.7e-02
    ## RAS / TGF-β / EMT（高合成、高侵袭且葡萄糖需求高的状态）
    "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION", # NES=+1.43, padj=2.1e-03
    ## 增殖（对葡萄糖来源碳的生物量需求增加）
    "HALLMARK_MITOTIC_SPINDLE", # NES=+1.34, padj=3.3e-02
    "GOBP_GLIAL_CELL_PROLIFERATION",
    "GOBP_EPITHELIAL_CELL_PROLIFERATION",
    "PULVER_FOREY_CELLCYCLE_ENRICHED_TFS_G2_M"
)
data.pathway <- data %>% filter(ID %in% pathway)


## 绘制图形
## NES 排序
data.pathway <- data.pathway %>%
    mutate(Description = str_to_title(Description)) %>%
    mutate(Description = str_replace_all(Description, "_", " ")) %>%
    mutate(
        group = ifelse(NES > 0, "up", "down"),
        Description = factor(Description, levels = Description[order(NES)])
    )
p <- data.pathway %>%
    ggplot(aes(x = NES, y = Description)) +
    geom_point(aes(color = -log10(p.adjust), shape = group), size = 3) +
    scale_color_gradientn(colors = c("green4", "yellowgreen", "yellow", "orange", "red")) +
    scale_shape_manual(values = c("up" = 16, "down" = 17)) +
    scale_x_continuous(expand = expansion(mult = 0.2)) +
    labs(x = "NES", y = "", color = "-log10(p.adjust)", shape = "Group") +
    theme_bw() +
    theme(
        panel.grid = element_blank(),
        axis.text.y = element_text(size = 8),
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 8)
    )
ggsave(paste0(output_prefix, "_gsea_h_plot_pair.pdf"), p, width = 4.5, height = 2.5)
ggsave(paste0(output_prefix, "_gsea_h_plot_pair.png"), p, width = 4.5, height = 2.5)
