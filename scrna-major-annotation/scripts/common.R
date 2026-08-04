#!/usr/bin/env Rscript
# 内部辅助脚本：仅由本目录中的主脚本 source，不由用户直接运行。

`%||%` <- function(x, y) if (is.null(x) || length(x) == 0) y else x

WORKFLOW_SEED <- 666L
WORKFLOW_RESOLUTIONS <- seq(0.1, 1.5, by = 0.1)
WORKFLOW_DIMS <- 1:30
WORKFLOW_NFEATURES <- 2000L
WORKFLOW_LABELS <- c(
  "Epithelial", "T_NK_cell", "B_cell", "Plasma", "Myeloid", "pDC",
  "Mast", "Neutrophil", "Fibroblast", "Endothelial", "Mural", "Proliferation"
)
WORKFLOW_QC_METRICS <- c(
  "nFeature_RNA", "nCount_RNA", "percent_mito", "percent_ribo",
  "percent_hb", "S.Score", "G2M.Score", "decontX_contamination"
)
WORKFLOW_HARD_QC_METRICS <- c("nFeature_RNA", "nCount_RNA", "percent_mito", "decontX_contamination")
WORKFLOW_PROLIFERATION_GENES <- c("MKI67", "TOP2A", "UBE2C", "CENPF", "TYMS", "STMN1")

parse_cli_args <- function(args = commandArgs(trailingOnly = TRUE)) {
  out <- list()
  i <- 1L
  while (i <= length(args)) {
    token <- args[[i]]
    if (!startsWith(token, "--")) stop("Unexpected positional argument: ", token)
    key <- gsub("-", "_", substring(token, 3L), fixed = TRUE)
    if (!nzchar(key)) stop("Empty command-line option.")
    if (i == length(args) || startsWith(args[[i + 1L]], "--")) {
      out[[key]] <- TRUE
      i <- i + 1L
    } else {
      out[[key]] <- args[[i + 1L]]
      i <- i + 2L
    }
  }
  out
}

get_opt <- function(opts, key, default = NULL) opts[[gsub("-", "_", key, fixed = TRUE)]] %||% default

required_opt <- function(opts, key) {
  value <- get_opt(opts, key)
  if (is.null(value) || length(value) != 1L || !nzchar(as.character(value))) {
    stop("Missing required option --", gsub("_", "-", key, fixed = TRUE), ".")
  }
  as.character(value)
}

parse_integer <- function(x, name, min_value = NULL, max_value = NULL) {
  value <- suppressWarnings(as.integer(x))
  if (length(value) != 1L || is.na(value) || (!is.null(min_value) && value < min_value) ||
      (!is.null(max_value) && value > max_value)) {
    stop("Invalid ", name, ": ", x)
  }
  value
}

ensure_dir <- function(path) {
  if (!dir.exists(path) && !dir.create(path, recursive = TRUE, showWarnings = FALSE)) {
    stop("Cannot create directory: ", path)
  }
  normalizePath(path, mustWork = TRUE)
}

absolute_path <- function(path, must_work = FALSE) {
  path <- path.expand(path)
  if (must_work && !file.exists(path)) stop("Path does not exist: ", path)
  if (must_work) return(normalizePath(path, mustWork = TRUE))
  parent <- ensure_dir(dirname(path))
  file.path(parent, basename(path))
}

assert_output_not_input <- function(input, output) {
  input <- normalizePath(input, mustWork = TRUE)
  output <- absolute_path(output, must_work = FALSE)
  if (identical(input, output)) stop("Input and output paths must differ: ", input)
  invisible(TRUE)
}

replace_file <- function(tmp, target) {
  target <- absolute_path(target, must_work = FALSE)
  if (file.rename(tmp, target)) return(invisible(target))
  if (!file.copy(tmp, target, overwrite = TRUE)) stop("Cannot write output: ", target)
  unlink(tmp)
  invisible(target)
}

atomic_save_rds <- function(object, path, compress = TRUE) {
  path <- absolute_path(path, must_work = FALSE)
  tmp <- tempfile(pattern = paste0(".", basename(path), "."), tmpdir = dirname(path), fileext = ".tmp")
  on.exit(if (file.exists(tmp)) unlink(tmp), add = TRUE)
  saveRDS(object, tmp, compress = compress)
  replace_file(tmp, path)
}

atomic_write_csv <- function(x, path) {
  path <- absolute_path(path, must_work = FALSE)
  tmp <- tempfile(pattern = paste0(".", basename(path), "."), tmpdir = dirname(path), fileext = ".tmp")
  on.exit(if (file.exists(tmp)) unlink(tmp), add = TRUE)
  utils::write.csv(x, tmp, row.names = FALSE, na = "")
  replace_file(tmp, path)
}

atomic_write_tsv_gz <- function(x, path) {
  path <- absolute_path(path, must_work = FALSE)
  tmp <- tempfile(pattern = paste0(".", basename(path), "."), tmpdir = dirname(path), fileext = ".gz")
  on.exit(if (file.exists(tmp)) unlink(tmp), add = TRUE)
  con <- gzfile(tmp, open = "wt")
  on.exit(try(close(con), silent = TRUE), add = TRUE)
  utils::write.table(x, con, sep = "\t", quote = FALSE, row.names = FALSE, na = "")
  close(con)
  replace_file(tmp, path)
}

package_install_hint <- function(pkg) {
  bioc <- c("scDblFinder", "SingleCellExperiment", "BiocParallel", "decontX", "celda")
  github <- c(presto = "immunogenomics/presto")
  if (pkg %in% bioc) return(paste0('BiocManager::install("', pkg, '")'))
  if (pkg %in% names(github)) return(paste0('remotes::install_github("', github[[pkg]], '")'))
  paste0('install.packages("', pkg, '")')
}

require_packages <- function(packages, minimum = list()) {
  missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) {
    hints <- vapply(missing, package_install_hint, character(1))
    stop("Missing required package(s): ", paste(missing, collapse = ", "),
         "\nInstall with:\n", paste0("  ", hints, collapse = "\n"))
  }
  for (pkg in names(minimum)) {
    if (utils::packageVersion(pkg) < package_version(minimum[[pkg]])) {
      stop(pkg, " >= ", minimum[[pkg]], " is required; installed: ", utils::packageVersion(pkg), ".")
    }
  }
  invisible(TRUE)
}

package_versions <- function(packages) {
  out <- vapply(packages, function(pkg) {
    if (requireNamespace(pkg, quietly = TRUE)) as.character(utils::packageVersion(pkg)) else NA_character_
  }, character(1))
  as.list(out)
}

file_md5 <- function(path) {
  path <- normalizePath(path, mustWork = TRUE)
  unname(as.character(tools::md5sum(path)))
}

text_md5 <- function(x) {
  tmp <- tempfile(fileext = ".txt")
  on.exit(unlink(tmp), add = TRUE)
  writeLines(enc2utf8(x), tmp, useBytes = TRUE)
  unname(as.character(tools::md5sum(tmp)))
}

cell_set_md5 <- function(cells) text_md5(paste(sort(as.character(cells)), collapse = "\n"))

format_resolution <- function(x) formatC(as.numeric(x), format = "f", digits = 1)
resolution_column <- function(x) paste0("RNA_snn_res.", as.character(as.numeric(x)))

project_paths <- function(output_dir, project, round = NULL) {
  root <- ensure_dir(output_dir)
  if (is.null(round)) return(list(root = root, state = file.path(root, "workflow_state.json")))
  round_dir <- ensure_dir(file.path(root, paste0("round", round)))
  list(
    root = root,
    round_dir = round_dir,
    state = file.path(root, "workflow_state.json"),
    clustered = file.path(round_dir, paste0(project, ".round", round, ".clustered.rds")),
    resolution_prefix = file.path(round_dir, "resolution"),
    resolution_review = file.path(round_dir, "resolution-review.csv"),
    qc_review = file.path(round_dir, "qc-review.csv"),
    annotation_review = file.path(round_dir, "annotation-review.csv"),
    filtered = file.path(round_dir, paste0(project, ".round", round, ".filtered.rds")),
    preview = file.path(round_dir, "filter-preview.csv")
  )
}

assert_human_gene_symbols <- function(genes) {
  genes <- as.character(genes)
  if (!length(genes) || anyNA(genes) || any(!nzchar(genes)) || anyDuplicated(genes)) {
    stop("RNA feature names must be non-empty, unique human gene symbols.")
  }
  ensembl_fraction <- mean(grepl("^ENSG[0-9]+", genes))
  if (ensembl_fraction > 0.05) {
    stop("Input appears to use Ensembl IDs (", round(100 * ensembl_fraction, 1),
         "%). This workflow requires human gene symbols.")
  }
  if (!any(grepl("^MT-", genes))) {
    stop("No MT-* genes were found. Confirm GRCh38 human gene-symbol feature names.")
  }
  invisible(TRUE)
}

assert_seurat_input <- function(object, require_counts = TRUE, require_orig_ident = TRUE) {
  if (!inherits(object, "Seurat")) stop("Input RDS is not a Seurat object.")
  if (!"RNA" %in% names(object@assays)) stop("Input Seurat object lacks an RNA assay.")
  if (!ncol(object) || !nrow(object)) stop("Input Seurat object is empty.")
  if (anyDuplicated(colnames(object))) stop("Cell barcodes are not globally unique.")
  assert_human_gene_symbols(rownames(object[["RNA"]]))
  meta <- object[[]]
  if (require_orig_ident) {
    if (!"orig.ident" %in% colnames(meta) || anyNA(meta$orig.ident) || any(!nzchar(as.character(meta$orig.ident)))) {
      stop("Metadata must contain a non-empty orig.ident value for every cell.")
    }
  }
  if (require_counts) {
    layers <- SeuratObject::Layers(object[["RNA"]], search = "^counts")
    if (!length(layers)) stop("RNA assay has no counts layer.")
  }
  invisible(TRUE)
}

assert_joined_rna_layers <- function(object, require_data = FALSE) {
  expected <- c("counts", if (isTRUE(require_data)) "data")
  for (layer in expected) {
    found <- SeuratObject::Layers(object[["RNA"]], search = paste0("^", layer, "$"))
    if (length(found) != 1L) {
      stop("Expected exactly one joined RNA ", layer, " layer; found: ", paste(found, collapse = ", "), ".")
    }
  }
  invisible(TRUE)
}

assert_split_rna_layers <- function(object, sample_ids, require_data = FALSE) {
  sample_ids <- unique(as.character(sample_ids))
  if (length(sample_ids) < 2L) stop("Split-layer validation requires at least two samples.")
  expected <- c("counts", if (isTRUE(require_data)) "data")
  for (layer in expected) {
    found <- SeuratObject::Layers(object[["RNA"]], search = paste0("^", layer, "\\."))
    if (length(found) != length(sample_ids)) {
      stop(
        "Expected one RNA ", layer, " layer per orig.ident after split (", length(sample_ids),
        "); found ", length(found), ": ", paste(found, collapse = ", "), "."
      )
    }
  }
  invisible(TRUE)
}

join_rna_layers <- function(object, require_data = FALSE) {
  object <- SeuratObject::JoinLayers(object, assay = "RNA")
  assert_joined_rna_layers(object, require_data = require_data)
  object
}

rna_layer <- function(object, layer = c("counts", "data")) {
  layer <- match.arg(layer)
  layers <- SeuratObject::Layers(object[["RNA"]], search = paste0("^", layer, "$"))
  if (length(layers) != 1L) {
    stop("Expected one joined RNA ", layer, " layer; found: ", paste(layers, collapse = ", "),
         ". Run JoinLayers first.")
  }
  SeuratObject::LayerData(object[["RNA"]], layer = layers[[1]])
}

assert_positive_qc_counts <- function(object) {
  meta <- object[[]]
  count_layers <- SeuratObject::Layers(object[["RNA"]], search = "^counts$")
  if (length(count_layers) == 1L) {
    counts <- SeuratObject::LayerData(object[["RNA"]], layer = count_layers[[1]])
    n_count <- Matrix::colSums(counts)
    n_feature <- Matrix::colSums(counts > 0)
  } else {
    n_count <- meta$nCount_RNA
    n_feature <- meta$nFeature_RNA
  }
  bad <- which(n_count <= 0 | n_feature <= 0 | is.na(n_count) | is.na(n_feature))
  if (length(bad)) {
    stop(length(bad), " cells have zero/invalid nCount_RNA or nFeature_RNA. Resolve these cells explicitly; they are not silently removed.")
  }
  invisible(TRUE)
}

cell_cycle_genes <- function() {
  genes <- Seurat::cc.genes.updated.2019
  unique(c(genes$s.genes, genes$g2m.genes))
}

worker_count <- function() max(1L, min(4L, as.integer(parallelly::availableCores())))

future_parallel_plan <- function(workers = worker_count()) {
  if (.Platform$OS.type == "unix" && Sys.info()[["sysname"]] != "Windows") {
    future::plan(future::multicore, workers = workers)
  } else {
    future::plan(future::multisession, workers = workers)
  }
}

future_max_size_for <- function(object) {
  max(512 * 1024^2, as.numeric(object.size(object)) * 4)
}

clean_for_recluster <- function(object) {
  meta_cols <- colnames(object[[]])
  stale <- grep("(_snn_res\\.|^seurat_clusters$)", meta_cols, value = TRUE)
  for (column in stale) object[[column]] <- NULL
  object@reductions <- list()
  object@graphs <- list()
  object@neighbors <- list()
  if ("scale.data" %in% SeuratObject::Layers(object[["RNA"]])) {
    SeuratObject::LayerData(object[["RNA"]], layer = "scale.data") <- NULL
  }
  object
}

run_presto <- function(object, group_by) {
  require_packages(c("presto", "Matrix"))
  if (!group_by %in% colnames(object[[]])) stop("Missing grouping column: ", group_by)
  expression <- rna_layer(object, "data")
  groups <- factor(as.character(object[[group_by, drop = TRUE]]))
  if (nlevels(groups) < 2L) {
    return(data.frame(
      feature = character(), group = character(), auc = numeric(), pval = numeric(),
      padj = numeric(), pct_in = numeric(), pct_out = numeric(), avgExpr = numeric(),
      logFC = numeric(), stringsAsFactors = FALSE
    ))
  }
  result <- presto::wilcoxauc(X = expression, y = groups, verbose = FALSE)
  result$group <- as.character(result$group)
  result
}

filter_presto_markers <- function(markers, top_n = 15L) {
  required <- c("feature", "group", "auc", "padj", "pct_in", "logFC")
  missing <- setdiff(required, colnames(markers))
  if (length(missing)) stop("Presto result lacks: ", paste(missing, collapse = ", "))
  keep <- !is.na(markers$padj) & markers$padj <= 0.05 &
    !is.na(markers$pct_in) & markers$pct_in >= 10 &
    !is.na(markers$logFC) & markers$logFC > 0.25 &
    !is.na(markers$auc) & markers$auc > 0.5
  markers <- markers[keep, , drop = FALSE]
  if (!nrow(markers)) return(markers)
  split_markers <- split(markers, markers$group)
  out <- lapply(split_markers, function(x) {
    x <- x[order(-x$auc, x$padj, -x$logFC, x$feature), , drop = FALSE]
    utils::head(x, top_n)
  })
  do.call(rbind, out)
}

cluster_average <- function(object, features, group_by) {
  expression <- rna_layer(object, "data")
  features <- unique(intersect(features, rownames(expression)))
  if (!length(features)) stop("None of the requested genes are present in the RNA data layer.")
  groups <- as.character(object[[group_by, drop = TRUE]])
  levels <- sort(unique(groups))
  matrix <- vapply(levels, function(level) {
    Matrix::rowMeans(expression[features, groups == level, drop = FALSE])
  }, numeric(length(features)))
  if (is.null(dim(matrix))) matrix <- matrix(matrix, ncol = 1L)
  rownames(matrix) <- features
  colnames(matrix) <- levels
  matrix
}

row_zscore <- function(x) {
  means <- rowMeans(x)
  sds <- apply(x, 1L, stats::sd)
  sds[is.na(sds) | sds == 0] <- 1
  sweep(sweep(x, 1L, means, "-"), 1L, sds, "/")
}

heatmap_plot <- function(z, title = NULL) {
  df <- as.data.frame(as.table(z), stringsAsFactors = FALSE)
  colnames(df) <- c("gene", "cluster", "z")
  df$gene <- factor(df$gene, levels = rev(rownames(z)))
  df$cluster <- factor(df$cluster, levels = colnames(z))
  ggplot2::ggplot(df, ggplot2::aes(cluster, gene, fill = z)) +
    ggplot2::geom_tile() +
    ggplot2::scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0) +
    ggplot2::labs(x = "Cluster", y = NULL, fill = "row z", title = title) +
    ggplot2::theme_minimal(base_size = 9) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
      axis.text.y = ggplot2::element_text(size = axis_text_size(nrow(z), base = 7.5, minimum = 4.0)),
      panel.grid = ggplot2::element_blank()
    )
}

clamp_value <- function(x, lower, upper) max(lower, min(upper, as.numeric(x)))

axis_text_size <- function(n_categories, base = 8, minimum = 4.5) {
  clamp_value(base - 0.12 * max(0, n_categories - 10), minimum, base)
}

grid_plot_dimensions <- function(n_panels, ncol = 2L, panel_width = 5.0,
                                 panel_height = 3.8, max_width = 18, max_height = 20) {
  ncol <- min(as.integer(ncol), as.integer(n_panels))
  nrow <- ceiling(n_panels / ncol)
  list(
    width = clamp_value(ncol * panel_width, 7, max_width),
    height = clamp_value(nrow * panel_height, 6, max_height),
    ncol = ncol,
    nrow = nrow
  )
}

discrete_grid_dimensions <- function(n_panels, n_categories, ncol = 2L) {
  per_panel_width <- clamp_value(4.8 + 0.18 * max(0, n_categories - 8), 4.8, 9)
  grid_plot_dimensions(n_panels, ncol, panel_width = per_panel_width, panel_height = 3.9)
}

heatmap_dimensions <- function(n_rows, n_columns) {
  list(
    width = clamp_value(5.5 + 0.55 * n_columns, 8, 18),
    height = clamp_value(4.5 + 0.10 * n_rows, 7, 20)
  )
}

dotplot_dimensions <- function(n_genes, n_clusters) {
  list(
    width = clamp_value(6 + 0.50 * n_clusters, 9, 18),
    height = clamp_value(5 + 0.14 * n_genes, 8, 18)
  )
}

umap_dimensions <- function(n_labels) {
  side <- clamp_value(8.5 + 0.12 * max(0, n_labels - 10), 9, 14)
  list(width = side, height = clamp_value(side * 0.78, 7, 11))
}

save_plot_pair <- function(plot, prefix, width, height, dpi = 180) {
  width <- clamp_value(width, 4, 20)
  height <- clamp_value(height, 4, 20)
  ggplot2::ggsave(paste0(prefix, ".png"), plot = plot, width = width, height = height, dpi = dpi, limitsize = TRUE)
  ggplot2::ggsave(paste0(prefix, ".pdf"), plot = plot, width = width, height = height, limitsize = TRUE)
}

check_cell_subset <- function(before, after) {
  before_cells <- colnames(before)
  after_cells <- colnames(after)
  if (anyDuplicated(after_cells) || !all(after_cells %in% before_cells)) {
    stop("Filtered output must contain a unique subset of input cell IDs.")
  }
  invisible(TRUE)
}
