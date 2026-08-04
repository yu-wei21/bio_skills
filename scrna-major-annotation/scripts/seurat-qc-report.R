#!/usr/bin/env Rscript
#
# 功能：按已批准分辨率生成 QC 图、统计证据和单一 QC 审批 CSV。
# 输入：
#   1. --input：本轮 clustered Seurat RDS；必须含 UMAP、orig.ident、8 个连续 QC 指标（含 decontX_contamination）及 scDblFinder.class。
#   2. --resolution-review：本轮 resolution-review.csv；同目录必须有 resolution-candidates.csv，且人工已批准。
#   3. --round、--project、--output-dir：必须与输入对象和审批绑定一致。
# 输出：
#   1. roundN/qc.metric-umap.png/.pdf：8 个连续 QC 指标 UMAP 加 1 个 singlet/doublet 分类 UMAP。
#   2. roundN/qc.violin-by-cluster.png/.pdf、qc.violin-by-sample.png/.pdf：各含 4 个 Round 1 可审批过滤指标的 violin，及一个 singlet/doublet 比例图。
#   3. roundN/qc-summary.csv、qc-review.csv；Round 1 承载阈值、cluster、doublet 决策，Round 2/3 仅承载 cluster 决策。
# 示例命令：
#   Rscript scripts/seurat-qc-report.R --input output/round1/prostate.round1.clustered.rds --resolution-review output/round1/resolution-review.csv --project prostate --round 1 --output-dir output

script_arg <- grep("^--file=", commandArgs(), value = TRUE)
script_dir <- if (length(script_arg)) dirname(normalizePath(sub("^--file=", "", script_arg[[1]]))) else getwd()
source(file.path(script_dir, "common.R"))
source(file.path(script_dir, "workflow-state.R"))
source(file.path(script_dir, "review-utils.R"))

qc_summary <- function(meta, group_col, group_type) {
  groups <- sort(unique(as.character(meta[[group_col]])))
  do.call(rbind, lapply(groups, function(group) {
    idx <- as.character(meta[[group_col]]) == group
    do.call(rbind, lapply(WORKFLOW_QC_METRICS, function(metric) {
      x <- as.numeric(meta[[metric]][idx])
      data.frame(
        group_type = group_type, group = group, metric = metric, n_cells = sum(idx),
        min = min(x, na.rm = TRUE), q01 = stats::quantile(x, 0.01, na.rm = TRUE),
        q05 = stats::quantile(x, 0.05, na.rm = TRUE), median = stats::median(x, na.rm = TRUE),
        q95 = stats::quantile(x, 0.95, na.rm = TRUE), q99 = stats::quantile(x, 0.99, na.rm = TRUE),
        max = max(x, na.rm = TRUE), mad = stats::mad(x, na.rm = TRUE), stringsAsFactors = FALSE
      )
    }))
  }))
}

doublet_plot <- function(meta, group_col, x_title, text_size = 8) {
  tab <- as.data.frame(table(
    group = as.character(meta[[group_col]]),
    class = factor(meta$scDblFinder.class, levels = c("singlet", "doublet"))
  ), stringsAsFactors = FALSE)
  totals <- aggregate(Freq ~ group, tab, sum)
  names(totals)[[2]] <- "total"
  tab <- merge(tab, totals, by = "group", all.x = TRUE)
  doublet_pct <- aggregate(Freq ~ group, tab[tab$class == "doublet", , drop = FALSE], sum)
  names(doublet_pct)[[2]] <- "doublet_n"
  totals <- merge(totals, doublet_pct, by = "group", all.x = TRUE)
  totals$doublet_n[is.na(totals$doublet_n)] <- 0
  totals$label <- sprintf("%.1f%% | n=%d", 100 * totals$doublet_n / totals$total, totals$total)
  ggplot2::ggplot(tab, ggplot2::aes(group, Freq, fill = class)) +
    ggplot2::geom_col(position = "fill") +
    ggplot2::geom_text(data = totals, ggplot2::aes(group, 1.03, label = label), inherit.aes = FALSE, size = 2.5, angle = 45, hjust = 0) +
    ggplot2::scale_y_continuous(labels = function(x) paste0(round(100 * x), "%"), expand = c(0, 0)) +
    ggplot2::coord_cartesian(ylim = c(0, 1.18), clip = "off") +
    ggplot2::scale_fill_manual(values = c(singlet = "#4C78A8", doublet = "#E45756"), drop = FALSE) +
    ggplot2::labs(x = x_title, y = "Cell proportion", fill = NULL, title = "scDblFinder class composition") +
    ggplot2::theme_classic(base_size = 9) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, size = text_size))
}

make_review <- function(object, project, round, input_file, cluster_col) {
  clusters <- data.frame(
    decision_type = "cluster", scope = "global", profile = "",
    cluster_id = sort(unique(as.character(object[[cluster_col, drop = TRUE]]))), metric = "",
    stringsAsFactors = FALSE
  )
  if (round == 1L) {
    threshold <- expand.grid(
      decision_type = "threshold", scope = "global", profile = c("general", "neutrophil"),
      cluster_id = "", metric = WORKFLOW_HARD_QC_METRICS, stringsAsFactors = FALSE
    )
    threshold <- threshold[c("decision_type", "scope", "profile", "cluster_id", "metric")]
    samples <- data.frame(
      decision_type = "doublet", scope = sort(unique(as.character(object$orig.ident))),
      profile = "", cluster_id = "", metric = "", stringsAsFactors = FALSE
    )
    review <- rbind(threshold, clusters, samples)
  } else {
    review <- clusters
  }
  review$agent_action <- ""
  review$agent_min <- ""
  review$agent_max <- ""
  review$agent_reason <- ""
  review$agent_evidence <- ""
  review$human_action <- ""
  review$human_min <- ""
  review$human_max <- ""
  review$human_reason <- ""
  review$approved <- "0"
  review$reviewer <- ""
  review$reviewed_at <- ""
  add_review_binding(review, review_binding(project, round, input_file, object))
}

opts <- parse_cli_args()
input_file <- normalizePath(required_opt(opts, "input"), mustWork = TRUE)
resolution_review <- normalizePath(required_opt(opts, "resolution_review"), mustWork = TRUE)
project <- required_opt(opts, "project")
round <- parse_integer(required_opt(opts, "round"), "round", 1L, 3L)
paths <- project_paths(required_opt(opts, "output_dir"), project, round)
candidates_file <- file.path(dirname(resolution_review), "resolution-candidates.csv")
summary_file <- file.path(paths$round_dir, "qc-summary.csv")
packages <- c("Seurat", "SeuratObject", "ggplot2", "patchwork", "jsonlite")

tryCatch({
  require_packages(packages, minimum = list(Seurat = "5.0.0", SeuratObject = "5.0.0"))
  assert_clustered_stage_input(paths$state, project, round, input_file)
  object <- readRDS(input_file)
  assert_seurat_input(object)
  missing <- setdiff(c(WORKFLOW_QC_METRICS, "scDblFinder.class"), colnames(object[[]]))
  if (length(missing)) stop("Missing QC metadata: ", paste(missing, collapse = ", "))
  if (anyNA(object$decontX_contamination) || any(!is.finite(object$decontX_contamination)) ||
      any(object$decontX_contamination < 0 | object$decontX_contamination > 1)) {
    stop("decontX_contamination must contain finite values in [0, 1].")
  }
  selected <- validate_resolution_review(resolution_review, candidates_file, project, round, input_file, object)
  cluster_col <- selected$column
  meta <- object[[]]
  violin_metrics <- unique(WORKFLOW_HARD_QC_METRICS)

  feature_plots <- Seurat::FeaturePlot(object, features = WORKFLOW_QC_METRICS, reduction = "umap", order = TRUE, combine = FALSE)
  doublet_umap <- Seurat::DimPlot(
    object, reduction = "umap", group.by = "scDblFinder.class",
    order = "doublet", cols = c(singlet = "#C7C7C7", doublet = "#D73027")
  ) + ggplot2::ggtitle("scDblFinder class") + ggplot2::theme(legend.position = "bottom")
  umap_plots <- c(feature_plots, list(doublet_umap))
  umap_size <- grid_plot_dimensions(length(umap_plots), ncol = 2L, panel_width = 5.0, panel_height = 3.8)
  feature_grid <- patchwork::wrap_plots(umap_plots, ncol = umap_size$ncol)
  save_plot_pair(feature_grid, file.path(paths$round_dir, "qc.metric-umap"), umap_size$width, umap_size$height)

  n_clusters <- length(unique(as.character(meta[[cluster_col]])))
  cluster_text_size <- axis_text_size(n_clusters)
  cluster_vln <- lapply(violin_metrics, function(metric) {
    Seurat::VlnPlot(object, features = metric, group.by = cluster_col, pt.size = 0) +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, size = cluster_text_size), legend.position = "none")
  })
  cluster_vln[[length(cluster_vln) + 1L]] <- doublet_plot(meta, cluster_col, "Cluster", cluster_text_size)
  cluster_size <- discrete_grid_dimensions(length(cluster_vln), n_clusters, ncol = 2L)
  save_plot_pair(patchwork::wrap_plots(cluster_vln, ncol = cluster_size$ncol), file.path(paths$round_dir, "qc.violin-by-cluster"), cluster_size$width, cluster_size$height)

  n_samples <- length(unique(as.character(meta$orig.ident)))
  sample_text_size <- axis_text_size(n_samples)
  sample_vln <- lapply(violin_metrics, function(metric) {
    Seurat::VlnPlot(object, features = metric, group.by = "orig.ident", pt.size = 0) +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, size = sample_text_size), legend.position = "none")
  })
  sample_vln[[length(sample_vln) + 1L]] <- doublet_plot(meta, "orig.ident", "Sample", sample_text_size)
  sample_size <- discrete_grid_dimensions(length(sample_vln), n_samples, ncol = 2L)
  save_plot_pair(patchwork::wrap_plots(sample_vln, ncol = sample_size$ncol), file.path(paths$round_dir, "qc.violin-by-sample"), sample_size$width, sample_size$height)

  summary <- rbind(qc_summary(meta, cluster_col, "cluster"), qc_summary(meta, "orig.ident", "sample"))
  atomic_write_csv(summary, summary_file)
  review <- make_review(object, project, round, input_file, cluster_col)
  write_review_template(review, paths$qc_review)
  record_workflow_event(
    paths$state, project, stage = "qc_report", status = "waiting_qc_review", round = round,
    script = "seurat-qc-report.R",
    inputs = list(clustered_rds = input_file, clustered_md5 = file_md5(input_file), resolution_review = resolution_review, resolution_review_md5 = file_md5(resolution_review)),
    outputs = list(qc_summary = summary_file, qc_review = paths$qc_review),
    metrics = list(
      cluster_col = cluster_col, n_cells = ncol(object), review_level = if (round == 1L) "cell_and_cluster" else "cluster_only",
      plot_dimensions = list(umap = umap_size, cluster = cluster_size, sample = sample_size)
    ), packages = packages
  )
  message("Created QC review bundle in: ", paths$round_dir)
}, error = function(e) {
  record_failure(paths$state, project, "qc_report", round, "seurat-qc-report.R", e,
                 inputs = list(clustered_rds = input_file), packages = packages)
  stop(e)
})
