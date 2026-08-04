#!/usr/bin/env Rscript
#
# 功能：生成大类细胞注释所需的 marker 统计、DotPlot、cluster-average heatmap 和审批模板。
# 输入：
#   1. --input：provisional 模式使用本轮 clustered RDS；final 模式使用本轮最终 filtered RDS，二者均须保留 RNA data 和批准的 cluster 列。
#   2. --mode：provisional 或 final；--resolution-review 为本轮已批准文件，同目录含 resolution-candidates.csv。
#   3. --project、--round、--output-dir；marker panel 固定为 references/marker-panels.tsv，物种为 human/GRCh38 gene symbols。
# 输出：
#   1. roundN/annotation.marker-dotplot.png/.pdf、annotation.top15-heatmap.png/.pdf、annotation-markers.tsv.gz。
#   2. provisional 输出 roundN/provisional-annotation-review.csv；final 输出 roundN/annotation-review.csv。
#   3. final 模式只准备正式审批，不写 celltype；既有 review CSV 默认覆盖。
# 示例命令：
#   Rscript scripts/seurat-annotation-report.R --mode provisional --input output/round1/prostate.round1.clustered.rds --resolution-review output/round1/resolution-review.csv --project prostate --round 1 --output-dir output
#   Rscript scripts/seurat-annotation-report.R --mode final --input output/round2/prostate.round2.filtered.rds --resolution-review output/round2/resolution-review.csv --project prostate --round 2 --output-dir output

script_arg <- grep("^--file=", commandArgs(), value = TRUE)
script_dir <- if (length(script_arg)) dirname(normalizePath(sub("^--file=", "", script_arg[[1]]))) else getwd()
skill_dir <- dirname(script_dir)
source(file.path(script_dir, "common.R"))
source(file.path(script_dir, "workflow-state.R"))
source(file.path(script_dir, "review-utils.R"))

opts <- parse_cli_args()
mode <- required_opt(opts, "mode")
if (!mode %in% c("provisional", "final")) stop("--mode must be provisional or final.")
input_file <- normalizePath(required_opt(opts, "input"), mustWork = TRUE)
resolution_review <- normalizePath(required_opt(opts, "resolution_review"), mustWork = TRUE)
project <- required_opt(opts, "project")
round <- parse_integer(required_opt(opts, "round"), "round", 1L, 3L)
paths <- project_paths(required_opt(opts, "output_dir"), project, round)
candidates_file <- file.path(dirname(resolution_review), "resolution-candidates.csv")
marker_panel_file <- file.path(skill_dir, "references", "marker-panels.tsv")
markers_file <- file.path(paths$round_dir, "annotation-markers.tsv.gz")
review_file <- if (mode == "final") paths$annotation_review else file.path(paths$round_dir, "provisional-annotation-review.csv")
packages <- c("Seurat", "SeuratObject", "presto", "ggplot2", "Matrix", "jsonlite")

tryCatch({
  require_packages(packages, minimum = list(Seurat = "5.0.0", SeuratObject = "5.0.0"))
  if (mode == "provisional") assert_clustered_stage_input(paths$state, project, round, input_file)
  object <- readRDS(input_file)
  assert_seurat_input(object)
  object <- join_rna_layers(object, require_data = TRUE)
  selected <- validate_resolution_for_current_object(resolution_review, candidates_file, project, round, input_file, object)
  cluster_col <- selected$column
  if (mode == "final") {
    state <- read_workflow_state(paths$state)
    event <- last_workflow_event(state, "filter_apply", round, "eligible_final_annotation")
    if (is.null(event)) stop("Final annotation requires an eligible_final_annotation filter result for this round.")
    if (!identical(normalizePath(input_file, mustWork = TRUE), normalizePath(event$outputs$filtered_rds, mustWork = TRUE)) ||
        !identical(file_md5(input_file), event$outputs$filtered_md5)) {
      stop("Final annotation input does not match the eligible filtered RDS and MD5 recorded in workflow state.")
    }
  }

  panel <- utils::read.delim(marker_panel_file, comment.char = "#", stringsAsFactors = FALSE, check.names = FALSE)
  require_review_columns(panel, c("panel", "marker_set", "gene", "role"), "Marker panel")
  expected_proliferation <- WORKFLOW_PROLIFERATION_GENES
  observed_proliferation <- panel$gene[panel$marker_set == "Proliferation"]
  if (!setequal(observed_proliferation, expected_proliferation)) stop("Proliferation marker panel does not match the fixed six-gene set.")
  present <- unique(panel$gene[panel$gene %in% rownames(object)])
  if (!length(present)) stop("No configured marker genes are present.")
  clusters <- sort(unique(as.character(object[[cluster_col, drop = TRUE]])))

  dot <- Seurat::DotPlot(object, features = present, group.by = cluster_col, assay = "RNA") +
    ggplot2::coord_flip() +
    ggplot2::labs(x = NULL, y = "Cluster", title = paste("Major-cell marker panel - round", round)) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(size = axis_text_size(length(clusters), base = 8, minimum = 5)),
      axis.text.y = ggplot2::element_text(size = axis_text_size(length(present), base = 7, minimum = 4.5))
    )
  dot_size <- dotplot_dimensions(length(present), length(clusters))
  save_plot_pair(dot, file.path(paths$round_dir, "annotation.marker-dotplot"), dot_size$width, dot_size$height)

  all_markers <- run_presto(object, cluster_col)
  top <- filter_presto_markers(all_markers, top_n = 15L)
  atomic_write_tsv_gz(top, markers_file)
  if (!nrow(top)) stop("No markers pass the fixed Presto filters; annotation cannot proceed.")
  z <- row_zscore(cluster_average(object, unique(top$feature), cluster_col))
  heat <- heatmap_plot(z, paste0("Round ", round, " cluster top-15 markers"))
  heat_size <- heatmap_dimensions(nrow(z), ncol(z))
  save_plot_pair(heat, file.path(paths$round_dir, "annotation.top15-heatmap"), heat_size$width, heat_size$height)

  review <- data.frame(
    cluster_col = cluster_col, cluster_id = clusters,
    agent_label = "", agent_reason = "", agent_evidence = "",
    human_label = "", approved = "0", reviewer = "", reviewed_at = "",
    stringsAsFactors = FALSE
  )
  review <- add_review_binding(review, review_binding(project, round, input_file, object))
  write_review_template(review, review_file)
  record_workflow_event(
    paths$state, project, stage = paste0("annotation_report_", mode),
    status = if (mode == "final") "waiting_annotation_review" else "provisional_annotation_ready",
    round = round, script = "seurat-annotation-report.R",
    inputs = list(input_rds = input_file, input_md5 = file_md5(input_file), resolution_review = resolution_review, resolution_review_md5 = file_md5(resolution_review)),
    outputs = list(markers = markers_file, review = review_file),
    metrics = list(
      cluster_col = cluster_col, n_clusters = length(clusters), n_cells = ncol(object),
      plot_dimensions = list(marker_dotplot = dot_size, top15_heatmap = heat_size)
    ), packages = packages
  )
  message("Created ", mode, " annotation bundle in: ", paths$round_dir)
}, error = function(e) {
  record_failure(paths$state, project, paste0("annotation_report_", mode), round, "seurat-annotation-report.R", e,
                 inputs = list(input_rds = input_file), packages = packages)
  stop(e)
})
