## 输入：
##   1. 已归一化 Seurat RDS：RNA data/data.* 层的 feature 名为 gene symbol，metadata 含无缺失的 RNA_snn_res.0.5 聚类列
## 输出：
##   1. output/seurat.4.annotation/ 下以输入 RDS 文件名为前缀的 marker panel PDF/PNG
##   2. output/seurat.4.annotation/ 下以输入 RDS 文件名为前缀的 top-marker heatmap PDF/PNG
## 示例命令：Rscript seurat.4.annotation.r data/pbmc.rds

library(Seurat)
library(tidyverse)
library(viridis)

## 参数解析
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1) {
  stop("Usage: Rscript seurat.4.annotation.r <input.rds>")
}
file.sc <- args[1]
output_dir <- file.path("output", "seurat.4.annotation")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
input_tag <- tools::file_path_sans_ext(basename(file.sc))
file.out <- file.path(output_dir, input_tag)
cat("file is : ", file.sc, "\n")
## 按分辨率绘制 UMAP

## 读取文件
sc <- readRDS(file.sc)

## 步骤 1：固定当前聚类身份
## 后续所有 DotPlot / marker 统计 / DEG 计算都基于这一列 cluster 身份展开。
## 这里先把 cluster 转成按出现顺序排列的 factor，保证作图时行顺序稳定。
## -------------------------------------------------------------------------
sc$RNA_snn_res.0.5 <- fct_inseq(sc$RNA_snn_res.0.5)
Idents(sc) <- "RNA_snn_res.0.5"
cluster_col <- "RNA_snn_res.0.5"

## -------------------------------------------------------------------------
## 步骤 2：定义用于注释的 marker 面板
## 最外层 list 代表不同主题的图：
##   1) all_markers: 全局大类群快速浏览
##   2) *_marker: 某一谱系内部的亚群/状态 marker
## 内层 list 代表图中的列分组，这些名字会显示在每个分组上方作为标题。
## -------------------------------------------------------------------------
marker_panels <- list(
  all_markers = list(
    Epithelial = c("EPCAM", "KRT8", "KRT18"),
    Pan_immune = c("PTPRC"),
    T_cell = c("CD3D", "CD3E", "TRAC"),
    NK_cell = c("NKG7", "GNLY", "KLRD1"),
    B_cell = c("CD19", "MS4A1", "CD79A"),
    Plasma = c("MZB1", "JCHAIN", "SDC1"),
    Pan_myeloid = c("LYZ", "LST1", "TYROBP"),
    Monocyte = c("FCN1", "VCAN", "CD14"),
    macro = c("CD68", "FCGR1A", "CSF1R"),
    Mast = c("KIT", "TPSAB1", "CPA3"),
    pDC = c("IL3RA", "LILRA4", "GZMB"),
    Neutrophil = c("CSF3R", "FCGR3B", "CXCR2"),
    Fibroblast = c("COL1A1", "DCN", "LUM"),
    Endothelial = c("PECAM1", "VWF", "PLVAP"),
    Pericyte = c("RGS5", "PDGFRB", "CSPG4"),
    Smooth_muscle = c("ACTA2", "ACTG2", "MYH11", "CNN1"),
    Proliferation = c("MKI67", "TOP2A", "HMGB2"),
    Erythrocyte = c("HBA1", "HBA2", "HBB"),
    Platelet = c("PPBP", "PF4"),
    Immunoglobulin = c("IGKC", "IGLC2", "IGHG1", "IGHM")
  ),
  epi_marker = list(
    Pan_epi = c("EPCAM"),
    Luminal = c("KRT8", "KRT18", "KRT19"),
    Basal = c("KRT5", "KRT15", "TP63"),
    Club_like = c("SCGB1A1", "PIGR", "LTF", "CP", "MMP7"),
    Prostate_luminal = c("KLK3", "KLK2", "ACPP", "AR", "NKX3-1"),
    Tumor_like = c("AMACR", "PCA3", "FOLH1", "ERG", "ETV1", "ETV4", "ETV5", "FLI1"),
    NE_like = c("SYP", "CHGA", "INSM1", "NCAM1", "ASCL1"),
    Cycling = c("MKI67", "TOP2A", "HMGB2")
  ),
  T_subtype_marker = list(
    pan_T = c("CD3D", "CD3E", "TRAC"),
    CD4_CD8 = c("CD4", "CD8A", "CD8B"),
    gd_T = c("TRDC", "TRGC1", "TRGC2", "TRDV1", "TRDV2", "TRGV9"),
    MAIT = c("TRAV1-2", "SLC4A10", "KLRB1"),
    NK_cell = c("NKG7", "GNLY", "KLRD1", "NCAM1", "FCGR3A"),
    Treg = c("FOXP3", "IL2RA", "CTLA4", "TNFRSF18"),
    Tfh = c("BCL6", "ICOS", "CXCR5", "IL21"),
    Naive = c("CCR7", "SELL", "TCF7", "LEF1", "IL7R"),
    Trm = c("CD69", "ITGAE", "ZNF683", "ITGA1", "CXCR6"),
    Tem = c("GZMK", "CXCR3", "CCL5"),
    Tcm = c("AQP3", "ANXA1", "S100A4", "GPR183", "CD44"),
    cytotoxic = c("GZMA", "GZMH", "GZMB", "PRF1", "IFNG"),
    Exhaustion = c("HAVCR2", "LAG3", "TIGIT", "TOX", "PDCD1", "CXCL13"),
    Th1 = c("TBX21", "IFNG", "CCR5", "IL12RB2"),
    Th2 = c("GATA3", "IL4", "IL5", "CCR4", "PTGDR2"),
    Th17 = c("CCR6", "RORC", "IL23R", "IL17A"),
    Cycling = c("MKI67", "TOP2A", "HMGB2")
  ),
  B_subtype_marker = list(
    pan_B = c("CD19", "MS4A1", "CD79A", "CD79B"),
    Naive_B = c("TCL1A", "IGHD", "IGHM", "SELL", "CCR7"),
    switchedMemory_B = c("CD27", "TNFRSF13B", "IGHG1", "IGHG2", "IGHG3", "IGHG4", "IGHA1", "IGHA2", "IGHE"),
    GC_B = c("AICDA", "BCL6", "RGS13"),
    Plasma = c("PRDM1", "JCHAIN", "XBP1", "MZB1", "SDC1"),
    Breg = c("IL10", "CD1D", "CD5"),
    cycling = c("MKI67", "TOP2A", "HMGB2")
  ),
  mye_subtype_marker = list(
    pan_myeloid = c("LYZ", "LST1", "TYROBP", "ITGAM"),
    cDC1 = c("CLEC9A", "XCR1", "BATF3"),
    cDC2 = c("CD1C", "FCER1A", "CLEC10A"),
    Mature_DC = c("LAMP3", "FSCN1", "CCL19"),
    pDC = c("LILRA4", "IL3RA", "TCF4", "GZMB"),
    Mast = c("KIT", "TPSAB1", "TPSB2"),
    Neutrophil = c("CSF3R", "FCGR3B", "CXCR2"),
    CD14_mono = c("CD14", "FCN1", "VCAN", "S100A8"),
    CD16_mono = c("FCGR3A", "MS4A7", "CX3CR1"),
    macro = c("CD68", "FCGR1A", "CSF1R"),
    TREM2_APOE_macrophage = c(
      "TREM2", "APOE", "APOC1",
      "C1QA", "C1QB", "C1QC"
    ),
    SPP1_macrophage = c(
      "SPP1", "MMP9", "MMP12",
      "FN1", "VEGFA"
    ),
    IL1B_macrophage = c(
      "IL1B",
      "CXCL8", "PLAUR", "TREM1"
    ),
    FOLR2_resident_macrophage = c(
      "FOLR2", "MRC1", "CD163", "MSR1",
      "LYVE1", "SELENOP", "MAF", "MAFB",
      "MARCO", "STAB1"
    ),
    IFN_macrophage = c(
      "ISG15", "IFIT1", "IFIT2", "IFIT3",
      "IFI6", "IFI27", "MX1", "MX2", "CXCL10"
    ),
    MHCII_macrophage = c(
      "HLA-DPA1", "HLA-DPB1",
      "HLA-DQA1", "HLA-DQB1",
      "HLA-DRA", "HLA-DRB1",
      "CD74", "CIITA"
    )
  ),
  fib_subtype_marker = list(
    pan_fib = c("COL1A1", "LUM", "DCN", "PDGFRA"),
    resting = c("DPT", "PI16", "APOD"),
    myCAF = c("POSTN", "MMP11", "COL11A1", "LRRC15", "CTHRC1", "FAP", "ITGA11"),
    iCAF = c("CXCL12", "CXCL14", "IL6", "C3", "CFD", "PLA2G2A"),
    apCAF = c("CD74", "HLA-DRA", "HLA-DRB1", "HLA-DPA1", "HLA-DPB1", "HLA-DMA", "HLA-DMB", "CIITA", "IFI30"),
    ifnCAF = c(
      "CXCL9", "CXCL10", "CXCL11", "IDO1",
      "STAT1", "IRF1", "ISG15", "IFIT1", "MX1", "GBP1"
    ),
    cycling_CAF = c(
      "MKI67", "TOP2A", "STMN1"
    )
  ),
  ec_subtype_marker = list(
    pan_ec = c("PECAM1", "VWF"),
    Arteries = c("GJA5", "FBLN5", "GJA4"),
    Capillaries = c("CA4", "CD36", "RGCC"),
    Hypoxia = c("MT1X", "MT1E", "MT2A"),
    Lymphatics = c("PROX1", "LYVE1", "CCL21"),
    Tip_cell = c("COL4A1", "KDR", "ESM1"),
    Veins = c("ACKR1", "SELP", "CLU"),
    cycling = c(
      "MKI67", "TOP2A", "STMN1"
    )
  )
)

## -------------------------------------------------------------------------
## 步骤 3：准备气泡热图的公共参数
## 这里使用 ComplexHeatmap 自定义“圆点热图”，目的是实现：
##   1) cluster 在 y 轴
##   2) marker 在 x 轴
##   3) marker 按生物学类别分组，并在组间留白
##   4) 圆点颜色表示 scaled expression，大小表示表达细胞比例
## -------------------------------------------------------------------------
library(ComplexHeatmap)

# 表达量颜色映射：低表达为蓝色，中间为白色，高表达为红色。
col_fun <- circlize::colorRamp2(c(-0.8, 0, 2), c("#424da7", "#ffffff", "#dd2b19"))

# 固定 cluster 的显示顺序，后续矩阵重排和 top marker 排序都复用这一顺序。
cluster_levels <- levels(Idents(sc))

# 把嵌套 marker list 展平成两列表：
#   feature: 基因名
#   marker_group: 基因所属分组
# 这样后续既能保留基因顺序，也能直接用于列分组。
prepare_marker_panel <- function(marker_panel) {
  tibble(
    feature = unname(unlist(marker_panel, use.names = FALSE)),
    marker_group = rep(names(marker_panel), lengths(marker_panel))
  ) %>%
    dplyr::filter(!duplicated(feature)) %>%
    dplyr::mutate(marker_group = factor(marker_group, levels = names(marker_panel)))
}

# 图上分组标题把下划线替换为空格，读起来更自然。
format_group_title <- function(x) {
  gsub("_", " ", x)
}

# 自定义百分比图例：圆点越大表示该 cluster 内表达该基因的细胞比例越高。
lgd_list <- list(
  Legend(
    labels = c("0.0", "0.2", "0.4", "0.6", "0.8"),
    title = "Percentage",
    graphics = list(
      function(x, y, w, h) grid.circle(x = x, y = y, r = 0 * unit(3, "mm"), gp = gpar(fill = "black")),
      function(x, y, w, h) grid.circle(x = x, y = y, r = 0.2 * unit(3, "mm"), gp = gpar(fill = "black")),
      function(x, y, w, h) grid.circle(x = x, y = y, r = 0.4 * unit(3, "mm"), gp = gpar(fill = "black")),
      function(x, y, w, h) grid.circle(x = x, y = y, r = 0.6 * unit(3, "mm"), gp = gpar(fill = "black")),
      function(x, y, w, h) grid.circle(x = x, y = y, r = 0.8 * unit(3, "mm"), gp = gpar(fill = "black"))
    )
  )
)

# --------------------------------------------------------------------------
# 步骤 4：预计算一次 DotPlot 数据，再按 panel 切片复用
# 这样避免在 for 循环里对每个 panel 重复跑 DotPlot，提高大对象时的速度。
# --------------------------------------------------------------------------
marker_feature_order <- unique(unlist(lapply(marker_panels, unlist, use.names = FALSE), use.names = FALSE))
marker_feature_order <- marker_feature_order[marker_feature_order %in% rownames(sc)]
dotplot_df_all <- DotPlot(sc, features = marker_feature_order)$data

# --------------------------------------------------------------------------
# 步骤 5：逐个 marker panel 出图
# 每轮循环做的事情是：
#   1) 取出当前 panel 的基因和分组
#   2) 从预计算 DotPlot 结果里切出当前 panel 的数据
#   3) 整理成表达矩阵和比例矩阵
#   4) 用 ComplexHeatmap 画成分组气泡热图
# --------------------------------------------------------------------------
for (i in seq_along(marker_panels)) {
  marker_df <- prepare_marker_panel(marker_panels[[i]])
  marker_df <- marker_df %>%
    dplyr::filter(feature %in% rownames(sc)) %>%
    dplyr::mutate(marker_group = forcats::fct_drop(marker_group))

  # 当前 panel 在图上方显示的分组标题。
  group_titles <- levels(marker_df$marker_group) %>%
    format_group_title()

  # 如果这个 panel 的 marker 在对象里一个都不存在，就跳过该图。
  if (nrow(marker_df) == 0) {
    next
  }

  # 从一次性计算好的 DotPlot 结果中切出当前 panel，避免重复计算。
  # 同时再次显式固定 gene / cluster 顺序，保证最终矩阵顺序可控。
  df <- dotplot_df_all %>%
    dplyr::filter(features.plot %in% marker_df$feature) %>%
    dplyr::mutate(
      features.plot = factor(features.plot, levels = marker_df$feature),
      id = factor(id, levels = cluster_levels)
    )

  # 表达矩阵：
  # 行 = cluster，列 = gene，值 = avg.exp.scaled
  # 用于映射圆点颜色。
  exp_mat <- df %>%
    dplyr::select(-pct.exp, -avg.exp) %>%
    pivot_wider(names_from = id, values_from = avg.exp.scaled) %>%
    as.data.frame()
  row.names(exp_mat) <- exp_mat$features.plot
  exp_mat <- exp_mat[, -1] %>% as.matrix()
  exp_mat <- t(exp_mat) %>% as.matrix()
  exp_mat <- exp_mat[cluster_levels, marker_df$feature, drop = FALSE]

  # 百分比矩阵：
  # 行 = cluster，列 = gene，值 = pct.exp
  # 用于映射圆点大小。
  percent_mat <- df %>%
    dplyr::select(-avg.exp, -avg.exp.scaled) %>%
    pivot_wider(names_from = id, values_from = pct.exp) %>%
    as.data.frame()
  row.names(percent_mat) <- percent_mat$features.plot
  percent_mat <- percent_mat[, -1] %>% as.matrix()
  percent_mat <- t(percent_mat) %>% as.matrix()
  percent_mat <- percent_mat[cluster_levels, marker_df$feature, drop = FALSE]

  # 自定义绘图层：
  # 不画默认热图方格，改成在每个格子中心画一个圆。
  # 圆点半径按 pct.exp 缩放，填充颜色按 avg.exp.scaled 映射。
  layer_fun <- local({
    pm <- percent_mat
    em <- exp_mat
    function(j, ii, x, y, w, h, fill) {
      grid.rect(x = x, y = y, width = w, height = h, gp = gpar(col = NA, fill = NA))
      grid.circle(
        x = x, y = y,
        r = pindex(pm, ii, j) / 100 * unit(2, "mm"),
        gp = gpar(fill = col_fun(pindex(em, ii, j)), col = NA)
      )
    }
  })

  # 用 Heatmap 承担排版职责：
  #   1) 列按 marker_group 分片
  #   2) 分片之间留白
  #   3) 每个分片顶部显示标题
  # 最终视觉效果接近常见单细胞 marker dotplot。
  hp <- Heatmap(exp_mat,
    name = "Expression",
    heatmap_legend_param = list(title = "Expression"),
    col = col_fun,
    rect_gp = gpar(type = "none"),
    layer_fun = layer_fun,
    row_names_gp = gpar(fontsize = 11),
    row_names_side = "left",
    row_title = NULL,
    column_names_rot = 90,
    column_names_gp = gpar(fontsize = 10),
    border = "black",
    column_split = marker_df$marker_group,
    column_title = group_titles,
    column_title_gp = gpar(fontsize = 10, fontface = "bold"),
    column_title_rot = 45,
    column_gap = unit(3, "mm"),
    cluster_column_slices = FALSE,
    cluster_rows = FALSE,
    cluster_columns = FALSE
  )

  # 每个 panel 单独输出 pdf 和 png，文件名与 panel 名一致，便于后续挑选查看。
  pdf(paste0(file.out, ".", names(marker_panels)[i], ".pdf"), width = 18, height = 5)
  draw(hp, annotation_legend_list = lgd_list)
  dev.off()
  png(paste0(file.out, ".", names(marker_panels)[i], ".png"), width = 18, height = 5, units = "in", res = 300)
  draw(hp, annotation_legend_list = lgd_list)
  dev.off()
}

## -------------------------------------------------------------------------
## 步骤 6：计算每个 cluster 的 top DEG
## 前面的 marker 面板是“人工指定 marker”视角；
## 这里转为“数据驱动的 cluster 特异基因”视角，补充绘制一张 top DEG 热图。
## -------------------------------------------------------------------------
library(presto)

# 兼容 Seurat v5 的 Assay5 / layer 结构。
# 如果指定 layer 不存在但存在同名前缀的多个 layer，则先合并再计算。
run_presto_wilcox <- function(seu, group_by, assay = "RNA", layer = "data", groups_use = NULL) {
  assay_obj <- seu[[assay]]
  if (inherits(assay_obj, "Assay5")) {
    layer_names <- Layers(assay_obj)
    if (!(layer %in% layer_names)) {
      layer_matches <- grep(paste0("^", layer, "($|\\.)"), layer_names, value = TRUE)
      if (length(layer_matches) == 0) {
        stop("Layer '", layer, "' not found in assay '", assay, "'.")
      }
      if (length(layer_matches) > 1) {
        seu <- JoinLayers(seu, assay = assay, layers = layer_matches, new = layer)
      } else {
        layer <- layer_matches[[1]]
      }
    }
  }

  expr_mat <- GetAssayData(seu, assay = assay, layer = layer)
  groups <- FetchData(seu, vars = group_by)[[group_by]] %>% as.character()
  presto::wilcoxauc(expr_mat, y = groups, groups_use = groups_use)
}

# 对当前 cluster 身份做 Wilcoxon 差异分析，得到每个 cluster 的 marker 统计量。
markers <- run_presto_wilcox(sc, group_by = cluster_col, assay = "RNA", layer = "data")

# 每个 cluster 取 logFC 最高的前 15 个基因，并按 cluster 原始显示顺序重排。
markers.top15 <- markers %>%
  group_by(group) %>%
  slice_max(order_by = logFC, n = 15, with_ties = FALSE) %>%
  ungroup() %>%
  mutate(group = factor(group, levels = cluster_levels)) %>%
  arrange(group, desc(logFC))

# 为 DoHeatmap 预先 scale 这些 top markers，提升不同 cluster 间的对比度。
sc <- ScaleData(sc, features = markers.top15$feature)

# --------------------------------------------------------------------------
# 步骤 7：输出 top marker 热图
# 这一部分和前面的气泡热图互补：
#   前者强调“已知 marker 是否符合注释预期”
#   这里强调“每个 cluster 自动找出的 top 基因长什么样”
# --------------------------------------------------------------------------
p <- DoHeatmap(sc, features = markers.top15$feature, size = 3) +
  scale_fill_gradientn(colors = viridis(100)) +
  theme(
    axis.text.x = element_blank(),
    axis.text.y = element_text(size = 10),
    legend.position = "right",
    plot.title = element_text(size = 14, face = "bold")
  ) +
  labs(title = "Top DEGs Heatmap", x = "Genes", y = "Clusters")
ggsave(paste0(file.out, ".heatmap.top15.", cluster_col, ".pdf"), plot = p, width = 15, height = 25)
ggsave(paste0(file.out, ".heatmap.top15.", cluster_col, ".png"), plot = p, width = 15, height = 25)
