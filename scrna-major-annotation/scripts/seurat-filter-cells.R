#!/usr/bin/env Rscript
#
# 功能：预览或应用人工批准的 QC 决定；Round 1 支持 cell/cluster level，Round 2/3 仅支持 cluster level。
# 输入：
#   1. --mode：preview 或 apply；--input 为本轮 clustered Seurat RDS。
#   2. --resolution-review 与 --qc-review：均须绑定同一 project、round、输入 MD5 和 cell set；apply 要求每行人工明确批准。
#   3. Round 1 另须 --count-version-review：两个 count 分支都完成第一轮报告后，由人工明确选择 original 或 corrected；只有被选分支可继续。
#   4. --round、--project、--output-dir：Round 1 的阈值边界为 value<min 或 value>max；Round 2/3 禁止阈值和逐细胞 doublet 删除。
# 输出：
#   1. preview：roundN/filter-preview.csv；只计算逐规则和去重 union 影响，不写 RDS。
#   2. apply：roundN/<project>.roundN.filtered.rds、removed-cells.tsv.gz、filter-summary.csv，并更新 workflow_state.json。
#   3. Round 1 必须再聚类；Round 2/3 仅按删除 cluster 的 union 计算 5%/10% 触发，Round 3 不进入 Round 4。
# 示例命令：
#   Rscript scripts/seurat-filter-cells.R --mode preview --input output/count-versions/original/round1/prostate.round1.clustered.rds --resolution-review output/count-versions/original/round1/resolution-review.csv --qc-review output/count-versions/original/round1/qc-review.csv --count-version-review output/prostate.count-version-review.csv --project prostate --round 1 --output-dir output/count-versions/original
#   Rscript scripts/seurat-filter-cells.R --mode apply --input output/count-versions/original/round1/prostate.round1.clustered.rds --resolution-review output/count-versions/original/round1/resolution-review.csv --qc-review output/count-versions/original/round1/qc-review.csv --count-version-review output/prostate.count-version-review.csv --project prostate --round 1 --output-dir output/count-versions/original

script_arg <- grep("^--file=", commandArgs(), value = TRUE)
script_dir <- if (length(script_arg)) dirname(normalizePath(sub("^--file=", "", script_arg[[1]]))) else getwd()
source(file.path(script_dir, "common.R"))
source(file.path(script_dir, "workflow-state.R"))
source(file.path(script_dir, "review-utils.R"))

effective_text <- function(human, agent) {
  human <- trimws(as.character(human))
  ifelse(nzchar(human), human, trimws(as.character(agent)))
}

optional_bound <- function(x) {
  x <- trimws(tolower(as.character(x)))
  if (!nzchar(x) || identical(x, "none")) return(NA_real_)
  value <- suppressWarnings(as.numeric(x))
  if (is.na(value) || !is.finite(value)) stop("Threshold bound must be numeric, 'none', or blank before approval: ", x)
  value
}

evaluate_filter <- function(object, review, cluster_col, approved, round) {
  meta <- object[[]]
  n <- nrow(meta)
  cluster_rows <- review[review$decision_type == "cluster", , drop = FALSE]
  doublet_rows <- review[review$decision_type == "doublet", , drop = FALSE]
  threshold_rows <- review[review$decision_type == "threshold", , drop = FALSE]
  if (round >= 2L && (nrow(doublet_rows) || nrow(threshold_rows))) {
    stop("Round 2/3 permit cluster-level QC filtering only.")
  }
  cluster_rows$action <- if (approved) cluster_rows$human_action else effective_text(cluster_rows$human_action, cluster_rows$agent_action)
  doublet_rows$action <- if (approved) doublet_rows$human_action else effective_text(doublet_rows$human_action, doublet_rows$agent_action)
  if (nrow(threshold_rows)) {
    threshold_rows$min_value <- mapply(function(h, a) optional_bound(if (approved) h else effective_text(h, a)), threshold_rows$human_min, threshold_rows$agent_min)
    threshold_rows$max_value <- mapply(function(h, a) optional_bound(if (approved) h else effective_text(h, a)), threshold_rows$human_max, threshold_rows$agent_max)
  }

  clusters <- as.character(meta[[cluster_col]])
  samples <- as.character(meta$orig.ident)
  reasons <- vector("list", n)
  add_reason <- function(idx, reason) {
    if (any(idx)) for (j in which(idx)) reasons[[j]] <<- c(reasons[[j]], reason)
  }

  for (i in seq_len(nrow(cluster_rows))) {
    if (cluster_rows$action[[i]] == "drop") {
      category <- if (approved) cluster_rows$human_reason[[i]] else cluster_rows$agent_reason[[i]]
      add_reason(clusters == cluster_rows$cluster_id[[i]], paste0("cluster_drop:", category))
    }
  }
  for (i in seq_len(nrow(doublet_rows))) {
    if (doublet_rows$action[[i]] == "remove") {
      add_reason(samples == doublet_rows$scope[[i]] & meta$scDblFinder.class == "doublet",
                 paste0("scDblFinder.class:", doublet_rows$scope[[i]]))
    }
  }

  cluster_action <- setNames(cluster_rows$action, cluster_rows$cluster_id)
  profiles <- ifelse(cluster_action[clusters] == "neutrophil", "neutrophil", "general")
  combinations <- unique(data.frame(sample = samples, profile = profiles, stringsAsFactors = FALSE))
  for (k in seq_len(nrow(combinations))) {
    sample_id <- combinations$sample[[k]]
    profile <- combinations$profile[[k]]
    base_idx <- samples == sample_id & profiles == profile
    for (metric in WORKFLOW_HARD_QC_METRICS) {
      candidates <- threshold_rows[threshold_rows$profile == profile & threshold_rows$metric == metric, , drop = FALSE]
      sample_row <- candidates[candidates$scope == sample_id, , drop = FALSE]
      global_row <- candidates[candidates$scope == "global", , drop = FALSE]
      row <- if (nrow(sample_row)) sample_row[1, , drop = FALSE] else if (nrow(global_row)) global_row[1, , drop = FALSE] else NULL
      if (is.null(row)) next
      value <- as.numeric(meta[[metric]])
      if (!is.na(row$min_value[[1]])) {
        add_reason(base_idx & value < row$min_value[[1]], paste0("threshold:", row$scope[[1]], ":", profile, ":", metric, ":min"))
      }
      if (!is.na(row$max_value[[1]])) {
        add_reason(base_idx & value > row$max_value[[1]], paste0("threshold:", row$scope[[1]], ":", profile, ":", metric, ":max"))
      }
    }
  }
  matched <- vapply(reasons, function(x) paste(unique(x), collapse = ";"), character(1))
  list(remove = nzchar(matched), matched_rules = matched, cluster_action = cluster_action)
}

impact_table <- function(object, evaluation, round, decision_hash) {
  meta <- object[[]]
  removed <- evaluation$remove
  global_rate <- mean(removed)
  samples <- sort(unique(as.character(meta$orig.ident)))
  sample_rows <- do.call(rbind, lapply(samples, function(sample_id) {
    idx <- as.character(meta$orig.ident) == sample_id
    data.frame(row_type = "union", scope = sample_id, rule = "all_rules_union",
               n_input = sum(idx), n_remove = sum(removed[idx]), pct_remove = mean(removed[idx]),
               trigger_threshold = 0.10, triggered = mean(removed[idx]) >= 0.10, stringsAsFactors = FALSE)
  }))
  global <- data.frame(row_type = "union", scope = "global", rule = "all_rules_union",
                       n_input = nrow(meta), n_remove = sum(removed), pct_remove = global_rate,
                       trigger_threshold = 0.05, triggered = global_rate >= 0.05, stringsAsFactors = FALSE)
  rules <- sort(unique(unlist(strsplit(evaluation$matched_rules[nzchar(evaluation$matched_rules)], ";", fixed = TRUE))))
  rule_rows <- if (length(rules)) do.call(rbind, lapply(rules, function(rule) {
    hit <- vapply(strsplit(evaluation$matched_rules, ";", fixed = TRUE), function(x) rule %in% x, logical(1))
    data.frame(row_type = "rule", scope = "global", rule = rule,
               n_input = nrow(meta), n_remove = sum(hit), pct_remove = mean(hit),
               trigger_threshold = NA_real_, triggered = NA, stringsAsFactors = FALSE)
  })) else data.frame(row_type = character(), scope = character(), rule = character(), n_input = integer(), n_remove = integer(), pct_remove = numeric(), trigger_threshold = numeric(), triggered = logical())
  out <- rbind(global, sample_rows, rule_rows)
  trigger <- global$triggered || any(sample_rows$triggered)
  outcome <- if (round == 1L) "next_round_required" else if (round == 2L && trigger) "next_round_required" else if (round == 3L && trigger) "needs_manual_review" else "eligible_final_annotation"
  out$round <- round
  out$decision_hash <- decision_hash
  out$outcome <- outcome
  out
}

opts <- parse_cli_args()
mode <- required_opt(opts, "mode")
if (!mode %in% c("preview", "apply")) stop("--mode must be preview or apply.")
input_file <- normalizePath(required_opt(opts, "input"), mustWork = TRUE)
resolution_review <- normalizePath(required_opt(opts, "resolution_review"), mustWork = TRUE)
qc_review <- normalizePath(required_opt(opts, "qc_review"), mustWork = TRUE)
project <- required_opt(opts, "project")
round <- parse_integer(required_opt(opts, "round"), "round", 1L, 3L)
count_version_review <- if (round == 1L) {
  normalizePath(required_opt(opts, "count_version_review"), mustWork = TRUE)
} else {
  NULL
}
paths <- project_paths(required_opt(opts, "output_dir"), project, round)
assert_output_not_input(input_file, paths$filtered)
candidates_file <- file.path(dirname(resolution_review), "resolution-candidates.csv")
summary_file <- file.path(paths$round_dir, "filter-summary.csv")
removed_file <- file.path(paths$round_dir, "removed-cells.tsv.gz")
packages <- c("Seurat", "SeuratObject", "jsonlite")

tryCatch({
  require_packages(packages, minimum = list(Seurat = "5.0.0", SeuratObject = "5.0.0"))
  assert_clustered_stage_input(paths$state, project, round, input_file)
  assert_round_reports_ready(paths$state, project, round, input_file, resolution_review, qc_review)
  object <- readRDS(input_file)
  assert_seurat_input(object)
  if (round == 1L) validate_count_version_review(count_version_review, project, paths$state, object)
  selected <- validate_resolution_review(resolution_review, candidates_file, project, round, input_file, object)
  review <- validate_qc_review(qc_review, project, round, input_file, object, selected$column, require_approved = mode == "apply")
  decision_hash <- filter_decision_hash(resolution_review, qc_review, count_version_review)
  evaluation <- evaluate_filter(object, review, selected$column, approved = mode == "apply", round = round)
  impact <- impact_table(object, evaluation, round, decision_hash)

  if (mode == "preview") {
    atomic_write_csv(impact, paths$preview)
    record_workflow_event(
      paths$state, project, stage = "filter_preview", status = "waiting_qc_approval", round = round,
      script = "seurat-filter-cells.R", inputs = list(input_rds = input_file, qc_review = qc_review, count_version_review = count_version_review),
      outputs = list(filter_preview = paths$preview, decision_hash = decision_hash),
      metrics = list(n_input = ncol(object), n_remove = sum(evaluation$remove), outcome = impact$outcome[[1]]), packages = packages
    )
    message("Created preview only: ", paths$preview)
  } else {
    if (!file.exists(paths$preview)) stop("Run preview after the final review edit before apply.")
    prior <- utils::read.csv(paths$preview, stringsAsFactors = FALSE, check.names = FALSE)
    if (!length(unique(prior$decision_hash)) || !identical(as.character(unique(prior$decision_hash)), decision_hash)) {
      stop("filter-preview.csv does not match the current approved reviews; rerun preview.")
    }
    keep <- !evaluation$remove
    if (!any(keep)) stop("Approved decisions would remove every cell.")
    before_samples <- table(object$orig.ident)
    after_samples <- table(factor(as.character(object$orig.ident[keep]), levels = names(before_samples)))
    if (any(after_samples == 0L)) stop("Approved decisions would remove every cell from sample(s): ", paste(names(after_samples)[after_samples == 0L], collapse = ", "))
    clusters <- as.character(object[[selected$column, drop = TRUE]])
    protected <- names(evaluation$cluster_action)[evaluation$cluster_action %in% c("keep", "neutrophil")]
    emptied <- protected[vapply(protected, function(cluster) !any(keep & clusters == cluster), logical(1))]
    if (length(emptied)) stop("Approved decisions would empty kept cluster(s): ", paste(emptied, collapse = ", "))

    filtered <- subset(object, cells = colnames(object)[keep])
    check_cell_subset(object, filtered)
    filtered <- join_rna_layers(filtered, require_data = TRUE)
    filtered@misc$filter_provenance <- list(
      round = round, source_clustered_rds = input_file, source_md5 = file_md5(input_file),
      source_cell_set_md5 = cell_set_md5(colnames(object)), resolution_review = resolution_review,
      resolution_review_md5 = file_md5(resolution_review), qc_review = qc_review,
      qc_review_md5 = file_md5(qc_review), decision_hash = decision_hash,
      cluster_col = selected$column, filtered_after_clustering = any(evaluation$remove)
    )
    primary_reason <- vapply(strsplit(evaluation$matched_rules[evaluation$remove], ";", fixed = TRUE), function(rules) {
      categories <- character()
      if (any(startsWith(rules, "scDblFinder.class:"))) categories <- c(categories, "doublet")
      cluster_rules <- rules[startsWith(rules, "cluster_drop:")]
      if (length(cluster_rules)) categories <- c(categories, sub("^cluster_drop:", "", cluster_rules))
      if (any(startsWith(rules, "threshold:"))) categories <- c(categories, "qc_threshold")
      paste(unique(categories), collapse = ";")
    }, character(1))
    removed <- data.frame(
      barcode = colnames(object)[evaluation$remove],
      sample = as.character(object$orig.ident[evaluation$remove]),
      reason = primary_reason,
      matched_rules = evaluation$matched_rules[evaluation$remove],
      stringsAsFactors = FALSE
    )
    atomic_save_rds(filtered, paths$filtered)
    atomic_write_tsv_gz(removed, removed_file)
    atomic_write_csv(impact, summary_file)
    outcome <- impact$outcome[[1]]
    record_workflow_event(
      paths$state, project, stage = "filter_apply", status = outcome, round = round,
      script = "seurat-filter-cells.R",
      inputs = list(input_rds = input_file, input_md5 = file_md5(input_file), resolution_review = resolution_review, qc_review = qc_review, count_version_review = count_version_review, decision_hash = decision_hash),
      outputs = list(filtered_rds = paths$filtered, filtered_md5 = file_md5(paths$filtered), removed_cells = removed_file, filter_summary = summary_file),
      metrics = list(n_input = ncol(object), n_removed = nrow(removed), n_retained = ncol(filtered), global_remove_rate = mean(evaluation$remove), filter_level = if (round == 1L) "cell_and_cluster" else "cluster_only", outcome = outcome, rna_layers_joined_before_save = TRUE),
      packages = packages
    )
    message("Created: ", paths$filtered, " (", outcome, ")")
  }
}, error = function(e) {
  record_failure(paths$state, project, paste0("filter_", mode), round, "seurat-filter-cells.R", e,
                 inputs = list(input_rds = input_file), packages = packages)
  stop(e)
})
