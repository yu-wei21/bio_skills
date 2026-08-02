## 输入：GSEA结果XLSX；每个sheet须含区分大小写的ID、Description、NES、p.adjust列，后两列为数值，ID须匹配脚本内置pathway列表。
## 输出：output/GSEA.dotplot.pair/下以输入XLSX文件名开头的配对点图PDF/PNG。
## 示例：Rscript GSEA.dotplot.pair.r results/gsea_results.xlsx

library(tidyverse)
library(cowplot)
library(openxlsx)

## read file
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1) {
    stop("Usage: Rscript GSEA.dotplot.pair.r <gsea_results.xlsx>")
}
input_file <- args[1]
if (!file.exists(input_file)) {
    stop("Input XLSX does not exist: ", input_file)
}
output_dir <- file.path("output", "GSEA.dotplot.pair")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
input_tag <- tools::file_path_sans_ext(basename(input_file))
output_prefix <- file.path(output_dir, input_tag)
# read input_file as a list of data.frames from the sheets
sheets <- openxlsx::getSheetNames(input_file)
data_list <- lapply(sheets, function(s) read.xlsx(input_file, sheet = s))
data <- do.call(rbind, data_list) %>%
    as.data.frame()

## 筛选感兴趣的pathway
pathway <- c(
    ## OXPHOS suppression (core FDG mechanism) — NES negative, enriched in FDG-negative cells
    "HALLMARK_OXIDATIVE_PHOSPHORYLATION", # NES=-1.94, padj=2.9e-07
    "Oxidative phosphorylation", # NES=-2.00, padj=1.1e-05
    "Cellular Respiration (GO:0045333)", # NES=-2.05, padj=3.8e-04
    "KEGG_OXIDATIVE_PHOSPHORYLATION", # NES=-1.56, padj=1.5e-02
    "GOBP_OXIDATIVE_PHOSPHORYLATION", # NES=-1.51, padj=4.2e-02
    ## Hypoxia / HIF-1α activation (amplifies GLUT1/LDHA, stabilises glycolysis)
    "HALLMARK_HYPOXIA", # NES=+1.43, padj=2.1e-03
    "GROSS_HYPOXIA_VIA_ELK3_DN", # NES=+1.69, padj=2.6e-07
    ## TNF-α / NF-κB inflammatory signalling (transcriptional driver of glycolysis)
    "HALLMARK_TNFA_SIGNALING_VIA_NFKB", # NES=+1.71, padj=5.0e-09
    "HALLMARK_INFLAMMATORY_RESPONSE", # NES=+1.30, padj=4.7e-02
    ## RAS / TGF-β / EMT (anabolic, invasive state with high glucose demand)
    "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION", # NES=+1.43, padj=2.1e-03
    ## Proliferation (elevated biomass demand for glucose-derived carbons)
    "HALLMARK_MITOTIC_SPINDLE", # NES=+1.34, padj=3.3e-02
    "GOBP_GLIAL_CELL_PROLIFERATION",
    "GOBP_EPITHELIAL_CELL_PROLIFERATION",
    "PULVER_FOREY_CELLCYCLE_ENRICHED_TFS_G2_M"
)
data.pathway <- data %>% filter(ID %in% pathway)


## 绘制plot图
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
