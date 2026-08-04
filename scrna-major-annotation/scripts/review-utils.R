#!/usr/bin/env Rscript
# 内部辅助脚本：仅由本目录中的主脚本 source，不由用户直接运行。

read_review <- function(path) {
  if (!file.exists(path)) stop("Review CSV does not exist: ", path)
  utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE,
                  colClasses = "character", na.strings = character())
}

review_is_approved <- function(x) trimws(tolower(as.character(x))) %in% c("1", "true", "yes")

write_review_template <- function(x, path) {
  atomic_write_csv(x, path)
}

review_binding <- function(project, round, input_file, object) {
  data.frame(
    project = project,
    round = as.character(round),
    input_rds = normalizePath(input_file, mustWork = TRUE),
    input_md5 = file_md5(input_file),
    cell_set_md5 = cell_set_md5(colnames(object)),
    stringsAsFactors = FALSE
  )
}

add_review_binding <- function(x, binding) {
  cbind(binding[rep(1L, nrow(x)), , drop = FALSE], x, stringsAsFactors = FALSE)
}

require_review_columns <- function(x, required, name) {
  missing <- setdiff(required, colnames(x))
  if (length(missing)) stop(name, " lacks required column(s): ", paste(missing, collapse = ", "))
}

assert_review_binding <- function(review, project, round, input_file, object) {
  required <- c("project", "round", "input_rds", "input_md5", "cell_set_md5")
  require_review_columns(review, required, "Review")
  expected <- review_binding(project, round, input_file, object)
  for (column in required) {
    values <- unique(as.character(review[[column]]))
    if (length(values) != 1L || !identical(values, as.character(expected[[column]][[1]]))) {
      stop("Review binding mismatch for ", column, ". The review belongs to a different round or object.")
    }
  }
  invisible(TRUE)
}

validate_resolution_review <- function(path, candidates_path, project, round, input_file, object, require_approved = TRUE) {
  review <- read_review(path)
  require_review_columns(review, c(
    "candidate_md5", "agent_resolution", "agent_reason", "agent_evidence",
    "human_resolution", "approved", "reviewer", "reviewed_at"
  ), "Resolution review")
  if (nrow(review) != 1L) stop("Resolution review must contain exactly one row.")
  assert_review_binding(review, project, round, input_file, object)
  if (!identical(review$candidate_md5[[1]], file_md5(candidates_path))) {
    stop("Resolution review is not bound to the current candidate table.")
  }
  if (require_approved && !review_is_approved(review$approved[[1]])) stop("Resolution review is not approved.")
  if (require_approved && !nzchar(trimws(review$human_resolution[[1]]))) {
    stop("human_resolution must be filled explicitly; agent_resolution is never applied automatically.")
  }
  if (nzchar(trimws(review$agent_resolution[[1]]))) {
    agent_value <- suppressWarnings(as.numeric(review$agent_resolution[[1]]))
    if (is.na(agent_value)) stop("agent_resolution must be numeric when filled.")
  }
  value <- suppressWarnings(as.numeric(review$human_resolution[[1]]))
  if (is.na(value)) stop("human_resolution must be numeric.")
  candidates <- utils::read.csv(candidates_path, stringsAsFactors = FALSE, check.names = FALSE)
  if (!"resolution" %in% colnames(candidates)) stop("Candidate table lacks resolution.")
  allowed <- as.numeric(candidates$resolution)
  if (nzchar(trimws(review$agent_resolution[[1]]))) {
    agent_value <- as.numeric(review$agent_resolution[[1]])
    if (!any(abs(agent_value - allowed) < 1e-8)) stop("agent_resolution is not in this round's candidate set.")
  }
  if (!any(abs(value - allowed) < 1e-8)) stop("human_resolution is not in this round's candidate set.")
  column <- resolution_column(value)
  if (!column %in% colnames(object[[]])) stop("Approved cluster column is absent from the input object: ", column)
  list(value = value, column = column, review = review)
}

validate_resolution_for_current_object <- function(path, candidates_path, project, round, input_file, object) {
  provenance <- object@misc$filter_provenance %||% NULL
  if (is.null(provenance)) {
    return(validate_resolution_review(path, candidates_path, project, round, input_file, object))
  }
  source_file <- provenance$source_clustered_rds %||% ""
  if (!nzchar(source_file) || !file.exists(source_file) || !identical(file_md5(source_file), provenance$source_md5)) {
    stop("Filtered object cannot verify its source clustered RDS.")
  }
  if (!identical(normalizePath(path, mustWork = TRUE), normalizePath(provenance$resolution_review, mustWork = TRUE)) ||
      !identical(file_md5(path), provenance$resolution_review_md5)) {
    stop("Resolution review does not match the filter provenance.")
  }
  source <- readRDS(source_file)
  selected <- validate_resolution_review(path, candidates_path, project, round, source_file, source)
  check_cell_subset(source, object)
  if (!selected$column %in% colnames(object[[]])) stop("Filtered object lacks the approved cluster column.")
  selected
}

parse_human_bound <- function(x, field) {
  x <- trimws(tolower(as.character(x)))
  if (!nzchar(x)) stop(field, " must be filled explicitly with a number or 'none'.")
  if (identical(x, "none")) return(NA_real_)
  value <- suppressWarnings(as.numeric(x))
  if (is.na(value) || !is.finite(value)) stop(field, " must be a finite number or 'none'.")
  value
}

validate_qc_review <- function(path, project, round, input_file, object, cluster_col, require_approved = TRUE) {
  review <- read_review(path)
  require_review_columns(review, c(
    "decision_type", "scope", "profile", "cluster_id", "metric",
    "agent_action", "agent_min", "agent_max", "agent_reason", "agent_evidence",
    "human_action", "human_min", "human_max", "human_reason", "approved", "reviewer", "reviewed_at"
  ), "QC review")
  assert_review_binding(review, project, round, input_file, object)
  if (require_approved && any(!review_is_approved(review$approved))) stop("Every QC review row must be approved.")
  if (anyDuplicated(review[c("decision_type", "scope", "profile", "cluster_id", "metric")])) {
    stop("QC review contains duplicate decision keys.")
  }
  allowed_types <- c("threshold", "cluster", "doublet")
  if (any(!review$decision_type %in% allowed_types)) stop("Invalid QC decision_type.")
  if (anyNA(object$scDblFinder.class) || any(!object$scDblFinder.class %in% c("singlet", "doublet"))) {
    stop("scDblFinder.class must contain only singlet/doublet without missing values.")
  }

  threshold <- review[review$decision_type == "threshold", , drop = FALSE]
  if (round == 1L && !nrow(threshold)) stop("Round 1 QC review requires cell-level threshold rows.")
  if (round >= 2L && nrow(threshold)) stop("Round 2/3 permit cluster-level QC decisions only; threshold rows are forbidden.")
  if (round == 1L) {
    if (any(!threshold$metric %in% WORKFLOW_HARD_QC_METRICS)) {
      stop("Thresholds may use only: ", paste(WORKFLOW_HARD_QC_METRICS, collapse = ", "), ".")
    }
    if (any(!threshold$profile %in% c("general", "neutrophil"))) stop("Threshold profile must be general or neutrophil.")
    allowed_scopes <- c("global", sort(unique(as.character(object$orig.ident))))
    if (any(!threshold$scope %in% allowed_scopes)) stop("Threshold scope must be global or an exact orig.ident value.")
    required_global <- expand.grid(scope = "global", profile = c("general", "neutrophil"), metric = WORKFLOW_HARD_QC_METRICS, stringsAsFactors = FALSE)
    observed_global <- threshold[threshold$scope == "global", c("scope", "profile", "metric"), drop = FALSE]
    key <- function(x) apply(x, 1L, paste, collapse = "|")
    if (!setequal(key(required_global), key(observed_global))) stop("QC review must retain every global profile-by-metric threshold row.")
    if (require_approved) {
      mins <- vapply(seq_len(nrow(threshold)), function(i) parse_human_bound(threshold$human_min[[i]], "human_min"), numeric(1))
      maxs <- vapply(seq_len(nrow(threshold)), function(i) parse_human_bound(threshold$human_max[[i]], "human_max"), numeric(1))
      if (any(!is.na(mins) & !is.na(maxs) & mins > maxs)) stop("A human_min exceeds human_max.")
    }
  }

  cluster <- review[review$decision_type == "cluster", , drop = FALSE]
  if (any(nzchar(cluster$agent_action) & !cluster$agent_action %in% c("keep", "drop", "neutrophil", "review"))) {
    stop("Invalid cluster agent_action.")
  }
  expected_clusters <- sort(unique(as.character(object[[cluster_col, drop = TRUE]])))
  if (!identical(sort(cluster$cluster_id), expected_clusters)) stop("Cluster decisions must cover the current cluster set exactly.")
  if (require_approved && any(!cluster$human_action %in% c("keep", "drop", "neutrophil"))) {
    stop("Every cluster requires an explicit human_action: keep, drop, or neutrophil.")
  }
  dropped <- cluster$human_action == "drop"
  if (require_approved && any(dropped & !cluster$human_reason %in% c("low_quality", "doublet", "contaminant"))) {
    stop("Dropped clusters require human_reason low_quality, doublet, or contaminant.")
  }

  doublet <- review[review$decision_type == "doublet", , drop = FALSE]
  if (round >= 2L && nrow(doublet)) stop("Round 2/3 permit cluster-level QC decisions only; doublet removal rows are forbidden.")
  if (any(nzchar(doublet$agent_action) & !doublet$agent_action %in% c("keep", "remove", "review"))) {
    stop("Invalid doublet agent_action.")
  }
  expected_samples <- sort(unique(as.character(object$orig.ident)))
  if (round == 1L && !identical(sort(doublet$scope), expected_samples)) stop("Round 1 doublet decisions must cover orig.ident values exactly.")
  if (round == 1L && require_approved && any(!doublet$human_action %in% c("keep", "remove"))) {
    stop("Every sample requires an explicit doublet human_action: keep or remove.")
  }
  review
}

validate_annotation_review <- function(path, project, round, input_file, object, cluster_col, require_approved = TRUE) {
  review <- read_review(path)
  require_review_columns(review, c(
    "cluster_col", "cluster_id", "agent_label", "agent_reason", "agent_evidence",
    "human_label", "approved", "reviewer", "reviewed_at"
  ), "Annotation review")
  assert_review_binding(review, project, round, input_file, object)
  if (any(review$cluster_col != cluster_col)) stop("Annotation review cluster_col does not match the approved resolution.")
  expected <- sort(unique(as.character(object[[cluster_col, drop = TRUE]])))
  if (!identical(sort(review$cluster_id), expected) || anyDuplicated(review$cluster_id)) {
    stop("Annotation review must cover the current cluster set exactly once.")
  }
  if (require_approved && any(!review_is_approved(review$approved))) stop("Every annotation row must be approved.")
  labels <- trimws(review$human_label)
  agent_labels <- trimws(review$agent_label)
  if (any(nzchar(agent_labels) & !agent_labels %in% c(WORKFLOW_LABELS, "Ambiguous"))) {
    stop("Agent labels must use the 12 canonical labels or Ambiguous.")
  }
  if (require_approved && any(!nzchar(labels))) stop("Every cluster requires an explicit human_label; agent_label is never applied automatically.")
  if (require_approved && any(grepl("^Ambiguous", labels, ignore.case = TRUE))) stop("Ambiguous cannot be a final label.")
  if (require_approved && any(tolower(labels) %in% c("low_qc", "doublet", "contaminant"))) {
    stop("QC removal categories cannot be used as final celltype labels.")
  }
  review
}

validate_count_version_review <- function(path, project, state_file, object) {
  review <- read_review(path)
  required <- c(
    "project", "parent_qc_rds", "parent_qc_md5", "shared_cell_set_md5",
    "original_branch_dir", "corrected_branch_dir",
    "original_rds", "original_rds_md5", "corrected_rds", "corrected_rds_md5",
    "agent_version", "agent_reason", "agent_evidence",
    "human_version", "human_reason", "approved", "reviewer", "reviewed_at"
  )
  require_review_columns(review, required, "Count-version review")
  if (nrow(review) != 1L) stop("Count-version review must contain exactly one row.")
  if (!identical(review$project[[1]], project)) stop("Count-version review belongs to a different project.")
  for (column in c("parent_qc_rds", "original_branch_dir", "corrected_branch_dir", "original_rds", "corrected_rds")) {
    if (!file.exists(review[[column]][[1]]) && !dir.exists(review[[column]][[1]])) {
      stop("Count-version review path does not exist: ", review[[column]][[1]])
    }
  }
  if (!identical(file_md5(review$parent_qc_rds[[1]]), review$parent_qc_md5[[1]]) ||
      !identical(file_md5(review$original_rds[[1]]), review$original_rds_md5[[1]]) ||
      !identical(file_md5(review$corrected_rds[[1]]), review$corrected_rds_md5[[1]])) {
    stop("Count-version review MD5 binding no longer matches one or more input RDS files.")
  }
  agent_version <- trimws(review$agent_version[[1]])
  if (nzchar(agent_version) && !agent_version %in% c("original", "corrected")) {
    stop("agent_version must be original, corrected, or blank.")
  }
  human_version <- trimws(review$human_version[[1]])
  if (!review_is_approved(review$approved[[1]]) || !human_version %in% c("original", "corrected") ||
      !nzchar(trimws(review$human_reason[[1]])) || !nzchar(trimws(review$reviewer[[1]])) ||
      !nzchar(trimws(review$reviewed_at[[1]]))) {
    stop("Count-version review requires human_version, human_reason, reviewer, reviewed_at, and approved=1.")
  }

  branch_specs <- list(
    original = list(dir = review$original_branch_dir[[1]], rds = review$original_rds[[1]], md5 = review$original_rds_md5[[1]]),
    corrected = list(dir = review$corrected_branch_dir[[1]], rds = review$corrected_rds[[1]], md5 = review$corrected_rds_md5[[1]])
  )
  for (version in names(branch_specs)) {
    spec <- branch_specs[[version]]
    branch_state_file <- file.path(spec$dir, "workflow_state.json")
    branch_state <- read_workflow_state(branch_state_file)
    seed_event <- last_workflow_event(branch_state, "decontx", status = "ready_for_round1_clustering")
    clustering_event <- last_workflow_event(branch_state, "clustering", 1L, "waiting_resolution_review")
    qc_event <- last_workflow_event(branch_state, "qc_report", 1L)
    annotation_event <- last_workflow_event(branch_state, "annotation_report_provisional", 1L)
    if (is.null(seed_event) || is.null(clustering_event) || is.null(qc_event) || is.null(annotation_event)) {
      stop("Both count-version branches must complete Round 1 clustering, QC report, and provisional annotation before selection; incomplete branch: ", version, ".")
    }
    if (!identical(normalizePath(seed_event$outputs$decontx_rds, mustWork = TRUE), normalizePath(spec$rds, mustWork = TRUE)) ||
        !identical(seed_event$outputs$decontx_md5 %||% "", spec$md5) ||
        !identical(seed_event$outputs$cell_set_md5 %||% "", review$shared_cell_set_md5[[1]])) {
      stop("Count-version review does not match the DecontX seed object for branch: ", version, ".")
    }
    clustered_md5 <- clustering_event$outputs$clustered_md5 %||% ""
    if (!nzchar(clustered_md5) || !identical(qc_event$inputs$clustered_md5 %||% "", clustered_md5) ||
        !identical(annotation_event$inputs$input_md5 %||% "", clustered_md5)) {
      stop("Round 1 evidence is not bound to the same clustered object for branch: ", version, ".")
    }
  }

  current_state <- read_workflow_state(state_file)
  current_seed <- last_workflow_event(current_state, "decontx", status = "ready_for_round1_clustering")
  object_version <- object@misc$count_version$version %||% ""
  selected_spec <- branch_specs[[human_version]]
  if (is.null(current_seed) || !identical(object_version, human_version) ||
      !identical(normalizePath(current_seed$outputs$decontx_rds, mustWork = TRUE), normalizePath(selected_spec$rds, mustWork = TRUE))) {
    stop("The current branch/object is not the human-selected count version.")
  }
  invisible(human_version)
}

filter_decision_hash <- function(resolution_review, qc_review, count_version_review = NULL) {
  hashes <- c(file_md5(resolution_review), file_md5(qc_review))
  if (!is.null(count_version_review)) hashes <- c(hashes, file_md5(count_version_review))
  text_md5(paste(hashes, collapse = "\n"))
}
