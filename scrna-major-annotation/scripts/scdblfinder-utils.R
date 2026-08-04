#!/usr/bin/env Rscript
# 内部辅助脚本：仅由 seurat-qc-metrics.R source，不由用户直接运行。

run_scdblfinder <- function(object) {
  require_packages(c("scDblFinder", "SingleCellExperiment", "BiocParallel", "parallelly", "S4Vectors", "SummarizedExperiment"))
  counts <- rna_layer(object, "counts")
  samples <- as.character(object$orig.ident)
  sce <- SingleCellExperiment::SingleCellExperiment(
    assays = list(counts = counts),
    colData = S4Vectors::DataFrame(orig.ident = samples, row.names = colnames(object))
  )
  workers <- worker_count()
  backend_requested <- "SerialParam"
  if (length(unique(samples)) == 1L || workers == 1L) {
    param <- BiocParallel::SerialParam(RNGseed = WORKFLOW_SEED, progressbar = FALSE)
  } else if (.Platform$OS.type == "windows") {
    backend_requested <- "SnowParam"
    param <- try(BiocParallel::SnowParam(workers = workers, RNGseed = WORKFLOW_SEED, progressbar = FALSE), silent = TRUE)
  } else {
    backend_requested <- "MulticoreParam"
    param <- try(BiocParallel::MulticoreParam(workers = workers, RNGseed = WORKFLOW_SEED, progressbar = FALSE), silent = TRUE)
  }
  if (inherits(param, "try-error")) {
    message("BiocParallel backend ", backend_requested,
            " could not be initialized; using reproducible SerialParam.")
    param <- BiocParallel::SerialParam(RNGseed = WORKFLOW_SEED, progressbar = FALSE)
  }
  set.seed(WORKFLOW_SEED)
  run_model <- function(bp) scDblFinder::scDblFinder(
    sce, samples = "orig.ident", clusters = FALSE,
    BPPARAM = bp, verbose = FALSE
  )
  fit <- try(run_model(param), silent = TRUE)
  if (inherits(fit, "try-error") && !inherits(param, "SerialParam")) {
    message("BiocParallel backend ", backend_requested,
            " could not start; retrying scDblFinder with reproducible SerialParam.")
    param <- BiocParallel::SerialParam(RNGseed = WORKFLOW_SEED, progressbar = FALSE)
    fit <- run_model(param)
  } else if (inherits(fit, "try-error")) {
    stop(as.character(fit))
  }
  sce <- fit
  object$scDblFinder.score <- as.numeric(SummarizedExperiment::colData(sce)$scDblFinder.score)
  object$scDblFinder.class <- as.character(SummarizedExperiment::colData(sce)$scDblFinder.class)
  summary <- do.call(rbind, lapply(sort(unique(samples)), function(sample_id) {
    idx <- samples == sample_id
    data.frame(
      orig.ident = sample_id,
      n_cells = sum(idx),
      n_doublet = sum(object$scDblFinder.class[idx] == "doublet", na.rm = TRUE),
      doublet_rate = mean(object$scDblFinder.class[idx] == "doublet", na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  }))
  object@misc$scDblFinder <- list(
    package_version = as.character(utils::packageVersion("scDblFinder")),
    seed = WORKFLOW_SEED,
    samples_column = "orig.ident",
    mode = "random",
    automatic_dbr = TRUE,
    parallel_backend_requested = backend_requested,
    parallel_backend_used = class(param)[[1]],
    workers_requested = workers,
    workers_used = BiocParallel::bpnworkers(param),
    per_sample = summary
  )
  object
}
