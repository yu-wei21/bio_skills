#!/usr/bin/env Rscript
#
# 功能：比较固定 0.1–1.5 分辨率的 UMAP、clustree、marker heatmap 和 cluster 稳定性证据。
# 输入：
#   1. --input：本轮 seurat-cluster.R 输出的 Seurat v5 RDS；必须有 UMAP、RNA data 和 15 个 RNA_snn_res.* 列。
#   2. --round：1、2 或 3；--project 与 --output-dir 必须与聚类步骤一致。
#   3. Marker 使用 RNA data layer 的 presto::wilcoxauc；筛选固定为 padj<=0.05、pct_in>=10、logFC>0.25、auc>0.5、每群 top15。
# 输出：
#   1. roundN/resolution.umap-grid.png/.pdf、resolution.clustree.png/.pdf、resolution.marker-heatmaps.pdf。
#   2. roundN/resolution-candidates.csv、resolution-markers.tsv.gz、resolution-review.csv。
#   3. workflow_state.json 更新为 waiting_resolution_review；既有 resolution-review.csv 默认覆盖。
# 示例命令：
#   Rscript scripts/seurat-find-resolution.R --input output/round1/prostate.round1.clustered.rds --project prostate --round 1 --output-dir output

script_arg <- grep("^--file=", commandArgs(), value = TRUE)
script_dir <- if (length(script_arg)) dirname(normalizePath(sub("^--file=", "", script_arg[[1]]))) else getwd()
source(file.path(script_dir, "common.R"))
source(file.path(script_dir, "workflow-state.R"))
source(file.path(script_dir, "review-utils.R"))

opts <- parse_cli_args()
input_file <- normalizePath(required_opt(opts, "input"), mustWork = TRUE)
project <- required_opt(opts, "project")
round <- parse_integer(required_opt(opts, "round"), "round", 1L, 3L)
paths <- project_paths(required_opt(opts, "output_dir"), project, round)
candidates_file <- file.path(paths$round_dir, "resolution-candidates.csv")
markers_file <- file.path(paths$round_dir, "resolution-markers.tsv.gz")
review_file <- paths$resolution_review
packages <- c("Seurat", "SeuratObject", "presto", "clustree", "ggraph", "ggplot2", "patchwork", "Matrix", "jsonlite")

tryCatch({
  require_packages(packages, minimum = list(Seurat = "5.0.0", SeuratObject = "5.0.0"))
  assert_clustered_stage_input(paths$state, project, round, input_file)
  object <- readRDS(input_file)
  assert_seurat_input(object)
  object <- join_rna_layers(object, require_data = TRUE)
  if (!"umap" %in% names(object@reductions)) stop("Input object lacks UMAP.")
  columns <- vapply(WORKFLOW_RESOLUTIONS, resolution_column, character(1))
  missing <- setdiff(columns, colnames(object[[]]))
  if (length(missing)) stop("Missing fixed resolution columns: ", paste(missing, collapse = ", "))

  # clustree 0.5.1 resolves these legacy guide names from the global environment
  # when used with ggraph 2.x and ggplot2 4.x.
  assign("guide_edge_colourbar", ggraph::guide_edge_colourbar, envir = .GlobalEnv)
  assign("guide_edge_coloursteps", ggraph::guide_edge_coloursteps, envir = .GlobalEnv)

  umaps <- lapply(seq_along(columns), function(i) {
    Seurat::DimPlot(object, reduction = "umap", group.by = columns[[i]], label = TRUE, repel = TRUE) +
      ggplot2::ggtitle(paste0("resolution ", format_resolution(WORKFLOW_RESOLUTIONS[[i]]))) +
      ggplot2::theme(legend.position = "none")
  })
  umap_size <- grid_plot_dimensions(length(umaps), ncol = 3L, panel_width = 4.4, panel_height = 3.9)
  umap_grid <- patchwork::wrap_plots(umaps, ncol = umap_size$ncol)
  save_plot_pair(umap_grid, file.path(paths$round_dir, "resolution.umap-grid"), umap_size$width, umap_size$height)

  tree <- clustree::clustree(object[[]], prefix = "RNA_snn_res.")
  max_clusters <- max(vapply(columns, function(column) length(unique(object[[column, drop = TRUE]])), integer(1)))
  tree_width <- clamp_value(9 + 0.35 * max_clusters, 10, 18)
  tree_height <- clamp_value(7 + 0.25 * length(columns), 9, 14)
  save_plot_pair(tree, file.path(paths$round_dir, "resolution.clustree"), tree_width, tree_height)

  marker_results <- vector("list", length(columns))
  candidate_rows <- vector("list", length(columns))
  heatmap_plots <- vector("list", length(columns))
  heatmap_rows <- integer(length(columns))
  heatmap_columns <- integer(length(columns))
  for (i in seq_along(columns)) {
    column <- columns[[i]]
    all_markers <- run_presto(object, column)
    top <- filter_presto_markers(all_markers, top_n = 15L)
    top$resolution <- rep(WORKFLOW_RESOLUTIONS[[i]], nrow(top))
    marker_results[[i]] <- top
    groups <- as.character(object[[column, drop = TRUE]])
    counts <- table(groups)
    candidate_rows[[i]] <- data.frame(
      resolution = WORKFLOW_RESOLUTIONS[[i]], cluster_col = column,
      n_clusters = length(counts), min_cluster_cells = min(counts),
      median_cluster_cells = stats::median(as.numeric(counts)),
      n_passing_markers = nrow(top),
      n_clusters_with_passing_markers = length(unique(top$group)),
      stringsAsFactors = FALSE
    )
    if (nrow(top)) {
      genes <- unique(top$feature)
      z <- row_zscore(cluster_average(object, genes, column))
      heatmap_plots[[i]] <- heatmap_plot(z, paste0("Resolution ", format_resolution(WORKFLOW_RESOLUTIONS[[i]])))
      heatmap_rows[[i]] <- nrow(z)
      heatmap_columns[[i]] <- ncol(z)
    } else {
      heatmap_plots[[i]] <- ggplot2::ggplot() + ggplot2::theme_void() + ggplot2::ggtitle(
        paste0("Resolution ", format_resolution(WORKFLOW_RESOLUTIONS[[i]]), ": no markers pass fixed filters")
      )
    }
  }
  heatmap_size <- heatmap_dimensions(max(heatmap_rows, 1L), max(heatmap_columns, 1L))
  pdf_file <- file.path(paths$round_dir, "resolution.marker-heatmaps.pdf")
  grDevices::pdf(pdf_file, width = heatmap_size$width, height = heatmap_size$height, onefile = TRUE)
  pdf_open <- TRUE
  on.exit(if (isTRUE(pdf_open)) grDevices::dev.off(), add = TRUE)
  for (plot in heatmap_plots) print(plot)
  grDevices::dev.off()
  pdf_open <- FALSE

  candidates <- do.call(rbind, candidate_rows)
  markers <- do.call(rbind, marker_results)
  atomic_write_csv(candidates, candidates_file)
  atomic_write_tsv_gz(markers, markers_file)
  binding <- review_binding(project, round, input_file, object)
  review <- data.frame(
    candidate_md5 = file_md5(candidates_file),
    agent_resolution = "", agent_reason = "", agent_evidence = "",
    human_resolution = "", approved = "0", reviewer = "", reviewed_at = "",
    stringsAsFactors = FALSE
  )
  review <- add_review_binding(review, binding)
  write_review_template(review, review_file)
  record_workflow_event(
    paths$state, project, stage = "resolution", status = "waiting_resolution_review", round = round,
    script = "seurat-find-resolution.R",
    inputs = list(clustered_rds = input_file, clustered_md5 = file_md5(input_file)),
    outputs = list(candidates = candidates_file, review = review_file, umap_grid = paste0(file.path(paths$round_dir, "resolution.umap-grid"), ".png"), clustree = paste0(file.path(paths$round_dir, "resolution.clustree"), ".png"), heatmaps = pdf_file),
    metrics = list(
      n_candidates = nrow(candidates),
      plot_dimensions = list(
        umap = umap_size,
        clustree = list(width = tree_width, height = tree_height),
        marker_heatmap = heatmap_size
      )
    ), packages = packages
  )
  message("Created resolution review bundle in: ", paths$round_dir)
}, error = function(e) {
  record_failure(paths$state, project, "resolution", round, "seurat-find-resolution.R", e,
                 inputs = list(clustered_rds = input_file), packages = packages)
  stop(e)
})
