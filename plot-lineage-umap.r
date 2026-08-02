#!/usr/bin/env Rscript
## 输入：5 个已注释的谱系 Seurat RDS；路径、原始 celltype 名称、展示名称与颜色均在下方
##       rds_files、cluster_num、new_names、lineage_colors 中硬编码，运行前必须按项目替换。
## 输出：./output/lineage_umap.pdf（两页：各谱系 UMAP 与编号图例）。
## 示例：Rscript plot-lineage-umap.r
## 注意事项：项目特异脚本；每个对象须有 celltype 列和可用于 DimPlot 的 UMAP。另需系统提供
##           Ghostscript 命令 gs 用于合并 PDF。

library(Seurat)
library(patchwork)
library(ggplot2)

## input rds files
rds_files <- c(
  Tcell       = "./output/Tcell.annotation.rds",
  B           = "./output/B.annotation.rds",
  Myeloid     = "./data/myeloid.annotation.rds",
  Fibroblast  = "./data/fibroblast.annotation.rds",
  Endothelial = "./output/endothelial.annotation.rds"
)

## titles
titles <- c(
  Tcell       = "T cells",
  B           = "B cells",
  Myeloid     = "Myeloid cells",
  Fibroblast  = "Fibroblasts",
  Endothelial = "Endothelial cells"
)

## ---- cluster mapping (old celltype -> cluster number) ----
## ordered by cluster number ascending
cluster_num <- list(
  Tcell = c(
    "CCR7_CD4_Tn"        = 1L,
    "ANXA1_CD4_Tcm"      = 2L,
    "PDCD1_CD4_T"        = 3L,
    "FOXP3_CD4_Treg"     = 4L,
    "CCR7_CD8_Tn"        = 5L,
    "ZNF683_CD8_Trm_Tex" = 6L,
    "GZMK_CD8_Tem"       = 7L,
    "MKI67_CD8_T"        = 8L,
    "gdT"                = 9L,
    "NK"                 = 10L
  ),
  B = c(
    "Bnaive"     = 1L,
    "Bmem"       = 2L,
    "Bgc"        = 3L,
    "B_HSP"      = 4L,
    "Bgc_MKI67"  = 5L,
    "B_ISG15"    = 6L
  ),
  Myeloid = c(
    "cDC1"           = 1L,
    "cDC2"           = 2L,
    "Mature_DC"      = 3L,
    "CD14_mono"      = 4L,
    "CD16_mono"      = 5L,
    "IL1B_Macro"     = 6L,
    "FOLR2_Macro"    = 7L,
    "MKI67_myeloid"  = 8L,
    "PR_high_myeloid" = 9L
  ),
  Fibroblast = c(
    "FOS_CAF"     = 1L,
    "POSTN_CAF"   = 2L,
    "PLA2G2A_CAF" = 3L,
    "MKI67_CAF"   = 4L,
    "SMC"         = 5L,
    "pericyte"    = 6L
  ),
  Endothelial = c(
    "vein"      = 1L,
    "artery"    = 2L,
    "HSPx"      = 3L,
    "tip"       = 4L,
    "lymphatic" = 5L
  )
)

## ---- new display names (old celltype -> new name) ----
new_names <- list(
  Tcell = c(
    "CCR7_CD4_Tn"        = "CD4_Tn_CCR7",
    "ANXA1_CD4_Tcm"      = "CD4_Tcm_ANXA1",
    "PDCD1_CD4_T"        = "CD4_Tex_PDCD1",
    "FOXP3_CD4_Treg"     = "CD4_Treg_FOXAP3",
    "CCR7_CD8_Tn"        = "CD8_Tn_CCR7",
    "ZNF683_CD8_Trm_Tex" = "CD8_Trm_ZNF683",
    "GZMK_CD8_Tem"       = "CD8_Tem_GZMK",
    "MKI67_CD8_T"        = "CD8_T_MKI67",
    "gdT"                = "Tgd",
    "NK"                 = "NK"
  ),
  B = c(
    "Bnaive"     = "Bn_TCL1A",
    "Bmem"       = "Bmem_CD27",
    "Bgc"        = "Bgc_BCL6",
    "B_HSP"      = "B_HSP",
    "Bgc_MKI67"  = "Bgc_MKI67",
    "B_ISG15"    = "B_ISG15"
  ),
  Myeloid = c(
    "cDC1"           = "cDC1",
    "cDC2"           = "cDC2",
    "Mature_DC"      = "cDC3",
    "CD14_mono"      = "Mo_CD14",
    "CD16_mono"      = "Mo_CD16",
    "IL1B_Macro"     = "Ma_IL1B",
    "FOLR2_Macro"    = "Ma_FOLR2",
    "MKI67_myeloid"  = "My_MKI67",
    "PR_high_myeloid" = "My_RP"
  ),
  Fibroblast = c(
    "FOS_CAF"     = "CAF_FOS",
    "POSTN_CAF"   = "CAF_POSTN",
    "PLA2G2A_CAF" = "CAF_PLA2G2A",
    "MKI67_CAF"   = "CAF_MKI67",
    "SMC"         = "SMC",
    "pericyte"    = "Pericyte"
  ),
  Endothelial = c(
    "vein"      = "En_vein",
    "artery"    = "En_artery",
    "HSPx"      = "En_HSP",
    "tip"       = "En_tip",
    "lymphatic" = "En_lym"
  )
)

## ---- colors (in cluster-number order 1, 2, 3, ...) ----
lineage_colors <- list(
  Tcell = c(
    "#FDBF6F", "#6A3D9A", "#FF7F00", "#B2DF8A",
    "#FFFF99", "#CAB2D6", "#FB9A99", "#A6CEE3",
    "#E31A1C", "#33A02C"
  ),
  B = c(
    "#FC8D62", "#6A3D9A", "#8DA0CB", "#E78AC3",
    "#66C2A5", "#A6761D"
  ),
  Myeloid = c(
    "#6A3D9A", "#1F78B4", "#CAB2D6", "#FF7F00",
    "#33A02C", "#FDBF6F", "#E31A1C", "#B2DF8A",
    "#A6CEE3"
  ),
  Fibroblast = c(
    "#E7298A", "#66A61E", "#7570B3", "#D95F02",
    "#66C2A5", "#E6AB02"
  ),
  Endothelial = c(
    "#1B9E77", "#D95F02", "#7570B3", "#E7298A",
    "#A6761D"
  )
)

## ---- read objects & add cluster labels ----
objs <- list()
for (name in names(rds_files)) {
  message("Reading: ", rds_files[name])
  obj <- readRDS(rds_files[name])

  map_cl <- cluster_num[[name]]
  map_nm <- new_names[[name]]

  ## lookup: old celltype -> cluster number (character for factor)
  obj$cluster_num <- as.character(unname(map_cl[obj$celltype]))
  obj$new_name    <- unname(map_nm[obj$celltype])

  ## factor levels = cluster numbers 1..N
  obj$cluster_num <- factor(obj$cluster_num,
    levels = as.character(sort(unique(map_cl))))

  objs[[name]] <- obj
}

## ---- generate DimPlots (no legend, labels show cluster number) ----
plots <- list()
for (name in names(objs)) {
  obj  <- objs[[name]]
  cols <- lineage_colors[[name]]

  p <- DimPlot(obj, group.by = "cluster_num", cols = cols,
               label = TRUE, label.size = 6, repel = TRUE) +
    ggtitle(titles[name]) +
    NoLegend() +
    theme(
      plot.title   = element_text(hjust = 0.5, size = 14),
      aspect.ratio = 1,
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5)
    )
  plots[[name]] <- p
}

## ---- legend panel (sorted by cluster number 1→N) ----
legend_plots <- list()
for (name in names(objs)) {
  map_cl <- cluster_num[[name]]
  map_nm <- new_names[[name]]
  cols   <- lineage_colors[[name]]

  ## ensure sorted by cluster number ascending
  srt <- order(map_cl)
  df <- data.frame(
    num   = map_cl[srt],
    label = paste0(map_cl[srt], "  ", map_nm[srt]),
    color = cols[srt],
    y     = seq_along(map_cl)
  )
  df$y <- max(df$y) - df$y + 1  ## reverse: cluster 1 on top

  p <- ggplot(df, aes(x = 1, y = y)) +
    geom_point(color = df$color, size = 3) +
    geom_text(aes(label = label), hjust = 0, nudge_x = 0.12, size = 4) +
    labs(title = titles[name]) +
    xlim(0.5, 3) +
    theme_void() +
    theme(
      plot.title  = element_text(hjust = 0, size = 10, face = "bold",
                                 margin = margin(b = 1, t = 4)),
      plot.margin = margin(2, 2, 2, 2)
    )
  ## legend for page 2: larger fonts, more spacing
  p2 <- ggplot(df, aes(x = 1, y = y)) +
    geom_point(color = df$color, size = 4) +
    geom_text(aes(label = label), hjust = 0, nudge_x = 0.15, size = 4.5) +
    labs(title = titles[name]) +
    xlim(0.5, 3.5) +
    theme_void() +
    theme(
      plot.title  = element_text(hjust = 0, size = 12, face = "bold",
                                 margin = margin(b = 2, t = 6)),
      plot.margin = margin(4, 4, 4, 4)
    )
  legend_plots[[name]] <- p2
}

## ---- Page 1: 2x3 UMAP (6th panel = blank) ----
all_panels <- list(
  plots[["Tcell"]], plots[["B"]], plots[["Myeloid"]],
  plots[["Fibroblast"]], plots[["Endothelial"]], plot_spacer()
)
page1 <- wrap_plots(all_panels, nrow = 2, ncol = 3)

## ---- Page 2: Legends in 2 columns ----
page2 <- wrap_plots(legend_plots, ncol = 2)

## ---- save two-page PDF via ghostscript merge ----
tmp1 <- tempfile(fileext = ".pdf")
tmp2 <- tempfile(fileext = ".pdf")
outfile <- "./output/lineage_umap.pdf"

ggsave(tmp1, page1, width = 210, height = 140, units = "mm")
ggsave(tmp2, page2, width = 210, height = 297, units = "mm")

## merge
gs_cmd <- paste("gs -q -dNOPAUSE -dBATCH -sDEVICE=pdfwrite",
  paste0("-sOutputFile=", outfile),
  tmp1, tmp2)
system(gs_cmd, intern = TRUE)

unlink(c(tmp1, tmp2))
message("Saved: ", outfile)
