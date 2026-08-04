#!/usr/bin/env Rscript
#
# 功能：把最终轮人工批准的 cluster 标签写入固定 metadata 列 celltype，并完成工作流。
# 输入：
#   1. --input：最终轮 filtered Seurat RDS；必须是最近一次 eligible_final_annotation 的输出。
#   2. --resolution-review 与 --annotation-review：须绑定当前 round、对象、cluster 列与完整 cell set，且每行人工批准。
#   3. --project、--round、--output-dir；已有 celltype 默认在派生输出中直接覆盖，输入 RDS 不修改。
# 输出：
#   1. <output-dir>/<project>.final.annotated.rds：每个细胞恰有一个非空、非 Ambiguous 的 celltype。
#   2. <output-dir>/<project>.final.celltype-umap.png/.pdf 与 <project>.final.annotation-summary.csv。
#   3. <output-dir>/workflow_state.json 更新为 complete；自定义标签在 summary 中标记。
# 示例命令：
#   Rscript scripts/seurat-annotate-cells.R --input output/round2/prostate.round2.filtered.rds --resolution-review output/round2/resolution-review.csv --annotation-review output/round2/annotation-review.csv --project prostate --round 2 --output-dir output

script_arg <- grep("^--file=", commandArgs(), value = TRUE)
script_dir <- if (length(script_arg)) dirname(normalizePath(sub("^--file=", "", script_arg[[1]]))) else getwd()
source(file.path(script_dir, "common.R"))
source(file.path(script_dir, "workflow-state.R"))
source(file.path(script_dir, "review-utils.R"))

opts <- parse_cli_args()
input_file <- normalizePath(required_opt(opts, "input"), mustWork = TRUE)
resolution_review <- normalizePath(required_opt(opts, "resolution_review"), mustWork = TRUE)
annotation_review <- normalizePath(required_opt(opts, "annotation_review"), mustWork = TRUE)
project <- required_opt(opts, "project")
round <- parse_integer(required_opt(opts, "round"), "round", 1L, 3L)
paths <- project_paths(required_opt(opts, "output_dir"), project, round)
candidates_file <- file.path(dirname(resolution_review), "resolution-candidates.csv")
output_file <- file.path(paths$root, paste0(project, ".final.annotated.rds"))
plot_prefix <- file.path(paths$root, paste0(project, ".final.celltype-umap"))
summary_file <- file.path(paths$root, paste0(project, ".final.annotation-summary.csv"))
assert_output_not_input(input_file, output_file)
packages <- c("Seurat", "SeuratObject", "ggplot2", "jsonlite")

tryCatch({
  require_packages(packages, minimum = list(Seurat = "5.0.0", SeuratObject = "5.0.0"))
  state <- read_workflow_state(paths$state)
  if (is.null(state) || !identical(state$status %||% "", "waiting_annotation_review") ||
      as.integer(state$current_round %||% NA_integer_) != round) {
    stop("Workflow state is not waiting for final annotation review in this round.")
  }
  object <- readRDS(input_file)
  assert_seurat_input(object)
  input_cells <- colnames(object)
  selected <- validate_resolution_for_current_object(resolution_review, candidates_file, project, round, input_file, object)
  review <- validate_annotation_review(annotation_review, project, round, input_file, object, selected$column)
  labels <- setNames(trimws(review$human_label), review$cluster_id)
  clusters <- as.character(object[[selected$column, drop = TRUE]])
  celltype <- unname(labels[clusters])
  if (length(celltype) != ncol(object) || anyNA(celltype) || any(!nzchar(celltype)) || any(grepl("^Ambiguous", celltype, ignore.case = TRUE))) {
    stop("Final celltype assignment failed completeness checks.")
  }
  object$celltype <- celltype
  if (!identical(colnames(object), input_cells)) stop("Final annotation changed the cell set or order.")
  object <- join_rna_layers(object, require_data = TRUE)

  summary <- as.data.frame(table(celltype = object$celltype, orig.ident = object$orig.ident), stringsAsFactors = FALSE)
  summary$is_custom_label <- !summary$celltype %in% WORKFLOW_LABELS
  atomic_save_rds(object, output_file)
  atomic_write_csv(summary, summary_file)
  plot <- Seurat::DimPlot(object, reduction = "umap", group.by = "celltype", label = TRUE, repel = TRUE) +
    ggplot2::ggtitle(paste0(project, " final major-cell annotation"))
  umap_size <- umap_dimensions(length(unique(celltype)))
  save_plot_pair(plot, plot_prefix, umap_size$width, umap_size$height)
  custom <- sort(unique(summary$celltype[summary$is_custom_label]))
  record_workflow_event(
    paths$state, project, stage = "annotation_apply", status = "complete", round = round,
    script = "seurat-annotate-cells.R",
    inputs = list(input_rds = input_file, input_md5 = file_md5(input_file), annotation_review = annotation_review, annotation_review_md5 = file_md5(annotation_review)),
    outputs = list(final_rds = output_file, final_md5 = file_md5(output_file), celltype_umap_png = paste0(plot_prefix, ".png"), annotation_summary = summary_file),
    metrics = list(n_cells = ncol(object), n_labels = length(unique(celltype)), custom_labels = custom, plot_dimensions = umap_size, rna_layers_joined_before_save = TRUE), packages = packages
  )
  message("Created final annotated object: ", output_file)
  if (length(custom)) message("Custom labels without bundled downstream panels: ", paste(custom, collapse = ", "))
}, error = function(e) {
  record_failure(paths$state, project, "annotation_apply", round, "seurat-annotate-cells.R", e,
                 inputs = list(input_rds = input_file), packages = packages)
  stop(e)
})
