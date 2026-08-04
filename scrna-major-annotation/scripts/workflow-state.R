#!/usr/bin/env Rscript
# 内部辅助脚本：仅由本目录中的主脚本 source，不由用户直接运行。

read_workflow_state <- function(path) {
  if (!file.exists(path)) return(NULL)
  jsonlite::fromJSON(path, simplifyVector = FALSE)
}

last_workflow_event <- function(state, stage, round = NULL, status = NULL) {
  if (is.null(state) || is.null(state$history)) return(NULL)
  matches <- vapply(state$history, function(event) {
    same <- identical(event$stage %||% "", stage)
    if (!is.null(round)) same <- same && as.integer(event$round %||% NA_integer_) == as.integer(round)
    if (!is.null(status)) same <- same && identical(event$status %||% "", status)
    same
  }, logical(1))
  if (!any(matches)) return(NULL)
  state$history[[tail(which(matches), 1L)]]
}

assert_decontx_input <- function(state_file, project, input_file) {
  state <- read_workflow_state(state_file)
  if (is.null(state) || !identical(state$project %||% "", project)) {
    stop("Missing workflow state for project ", project, ". Run seurat-qc-metrics.R first.")
  }
  event <- last_workflow_event(state, "qc_metrics", status = "ready_for_decontx")
  expected_path <- event$outputs$qc_metrics_rds %||% ""
  expected_md5 <- event$outputs$qc_metrics_md5 %||% ""
  if (is.null(event) || !nzchar(expected_path) || !file.exists(expected_path)) {
    stop("No successful qc_metrics output authorizes DecontX.")
  }
  if (!identical(normalizePath(input_file, mustWork = TRUE), normalizePath(expected_path, mustWork = TRUE)) ||
      !identical(file_md5(input_file), expected_md5)) {
    stop("DecontX input does not match the qc-metrics RDS and MD5 recorded in workflow_state.json.")
  }
  invisible(TRUE)
}

assert_round_clustering_input <- function(state_file, project, round, input_file) {
  state <- read_workflow_state(state_file)
  if (is.null(state) || !identical(state$project %||% "", project)) {
    stop("Missing workflow state for project ", project, ". Run the preceding scripts first.")
  }
  if (round == 1L) {
    event <- last_workflow_event(state, "decontx", status = "ready_for_round1_clustering")
    qc_event <- last_workflow_event(state, "qc_metrics", status = "ready_for_decontx")
    expected_path <- event$outputs$decontx_rds %||% ""
    expected_md5 <- event$outputs$decontx_md5 %||% ""
    if (!identical(state$current_stage %||% "", "decontx") ||
        !identical(state$status %||% "", "ready_for_round1_clustering") ||
        is.null(qc_event) ||
        !identical(event$inputs$qc_metrics_md5 %||% "", qc_event$outputs$qc_metrics_md5 %||% "")) {
      stop("Round 1 requires DecontX output derived from the latest successful qc-metrics RDS.")
    }
  } else {
    event <- last_workflow_event(state, "filter_apply", round - 1L, "next_round_required")
    expected_path <- event$outputs$filtered_rds %||% ""
    expected_md5 <- event$outputs$filtered_md5 %||% ""
  }
  if (is.null(event) || !nzchar(expected_path) || !file.exists(expected_path)) {
    stop("No valid preceding output authorizes clustering round ", round, ".")
  }
  if (!identical(normalizePath(input_file, mustWork = TRUE), normalizePath(expected_path, mustWork = TRUE)) ||
      !identical(file_md5(input_file), expected_md5)) {
    stop("Round ", round, " input does not match the authorized preceding RDS and MD5.")
  }
  invisible(TRUE)
}

assert_clustered_stage_input <- function(state_file, project, round, input_file) {
  state <- read_workflow_state(state_file)
  if (is.null(state) || !identical(state$project %||% "", project)) stop("Missing workflow state for this project.")
  event <- last_workflow_event(state, "clustering", round, "waiting_resolution_review")
  if (is.null(event)) stop("No successful clustering event exists for round ", round, ".")
  expected_path <- event$outputs$clustered_rds %||% ""
  expected_md5 <- event$outputs$clustered_md5 %||% ""
  if (!nzchar(expected_path) || !file.exists(expected_path) ||
      !identical(normalizePath(input_file, mustWork = TRUE), normalizePath(expected_path, mustWork = TRUE)) ||
      !identical(file_md5(input_file), expected_md5)) {
    stop("Input does not match the clustered RDS and MD5 recorded for round ", round, ".")
  }
  invisible(TRUE)
}

assert_round_reports_ready <- function(state_file, project, round, input_file,
                                       resolution_review, qc_review) {
  state <- read_workflow_state(state_file)
  if (is.null(state) || !identical(state$project %||% "", project)) stop("Missing workflow state for this project.")
  input_md5 <- file_md5(input_file)
  resolution_md5 <- file_md5(resolution_review)
  qc_event <- last_workflow_event(state, "qc_report", round)
  annotation_event <- last_workflow_event(state, "annotation_report_provisional", round)
  if (is.null(qc_event) || is.null(annotation_event)) {
    stop("Round ", round, " requires both QC report and provisional annotation report before filtering.")
  }
  recorded_qc_review <- qc_event$outputs$qc_review %||% ""
  qc_ok <- identical(qc_event$inputs$clustered_md5 %||% "", input_md5) &&
    identical(qc_event$inputs$resolution_review_md5 %||% "", resolution_md5) &&
    nzchar(recorded_qc_review) && file.exists(recorded_qc_review) &&
    identical(normalizePath(recorded_qc_review, mustWork = TRUE), normalizePath(qc_review, mustWork = TRUE))
  annotation_ok <- identical(annotation_event$inputs$input_md5 %||% "", input_md5) &&
    identical(annotation_event$inputs$resolution_review_md5 %||% "", resolution_md5)
  if (!qc_ok || !annotation_ok) {
    stop("Round ", round, " report provenance does not match the current clustered RDS or resolution review.")
  }
  invisible(TRUE)
}

write_workflow_state <- function(state, path) {
  path <- absolute_path(path, must_work = FALSE)
  tmp <- tempfile(pattern = ".workflow-state-", tmpdir = dirname(path), fileext = ".json")
  on.exit(if (file.exists(tmp)) unlink(tmp), add = TRUE)
  jsonlite::write_json(state, tmp, auto_unbox = TRUE, pretty = TRUE, null = "null", na = "null")
  replace_file(tmp, path)
}

record_workflow_event <- function(state_file, project, stage, status, round = NULL,
                                  script, inputs = list(), outputs = list(), metrics = list(),
                                  packages = character(), error = NULL) {
  require_packages("jsonlite")
  state <- read_workflow_state(state_file)
  now <- format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")
  if (is.null(state)) {
    state <- list(
      schema_version = 1L,
      project = project,
      species = "human",
      genome = "GRCh38",
      seed = WORKFLOW_SEED,
      max_rounds = 3L,
      created_at = now,
      history = list()
    )
  }
  if (!identical(as.character(state$project), as.character(project))) {
    stop("workflow_state.json belongs to project '", state$project, "', not '", project, "'.")
  }
  event <- list(
    timestamp = now,
    script = script,
    stage = stage,
    round = if (is.null(round)) NULL else as.integer(round),
    status = status,
    inputs = inputs,
    outputs = outputs,
    metrics = metrics,
    package_versions = package_versions(packages),
    error = error
  )
  state$current_stage <- stage
  state$current_round <- if (is.null(round)) state$current_round %||% NULL else as.integer(round)
  state$status <- status
  state$updated_at <- now
  state$history[[length(state$history) + 1L]] <- event
  write_workflow_state(state, state_file)
  invisible(state)
}

record_failure <- function(state_file, project, stage, round, script, condition, inputs = list(), packages = character()) {
  try(record_workflow_event(
    state_file = state_file, project = project, stage = stage, status = "failed",
    round = round, script = script, inputs = inputs, packages = packages,
    error = conditionMessage(condition)
  ), silent = TRUE)
  invisible(NULL)
}
