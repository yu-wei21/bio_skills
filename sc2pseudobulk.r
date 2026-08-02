## 输入：Seurat RDS（RNA counts/counts.*为非负整数，基因名为唯一标准SYMBOL；orig.ident和指定细胞大类列无NA/空值）。
## 输出：output/sc2pseudobulk/下以输入RDS文件名开头的all-cell及逐细胞大类gene×sample CSV。
## 说明：样本×细胞大类组合中无细胞的样本，会从对应细胞大类的CSV列中剔除（不保留0填充列），
##       以避免下游DESeq2因样本总counts为0而失败；剔除情况会在运行信息中说明。
##       DESeq2.r 会自动取"表达矩阵列 ∩ group.csv 样本"的交集进行分析，
##       因此 group.csv 无需按细胞大类裁剪，保留全部样本即可。
## 示例：Rscript sc2pseudobulk.r data/pbmc.rds celltype

## CSV列名保留原始orig.ident。

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
})

## ----------------------------- Parameters ---------------------------------
args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 2) {
  stop(
    "Usage: Rscript sc2pseudobulk.r <input_rds> <celltype_col>"
  )
}

input_rds <- args[1]
celltype_col <- args[2]

if (!file.exists(input_rds)) {
  stop("Input RDS file does not exist: ", input_rds)
}
if (!nzchar(celltype_col)) {
  stop("celltype_col must not be empty.")
}

input_rds <- normalizePath(input_rds, mustWork = TRUE)
input_basename <- basename(tools::file_path_sans_ext(input_rds))
output_dir <- file.path("output", "sc2pseudobulk")

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
if (!dir.exists(output_dir)) {
  stop("Failed to create output directory: ", output_dir)
}
output_dir <- normalizePath(output_dir, mustWork = TRUE)

## ----------------------------- Helpers ------------------------------------
sanitize_filename <- function(label) {
  filename <- gsub("[/\\\\]", "_", label)
  filename <- gsub("[<>:\"|?*]", "_", filename)
  filename <- gsub("[[:cntrl:]]", "_", filename)
  filename <- gsub("[. ]+$", "_", filename)

  if (!nzchar(filename) || filename %in% c(".", "..")) {
    stop("Cell-class label cannot be converted to a valid filename: ", label)
  }
  filename
}

aggregate_by_sample <- function(counts, cell_indices, sample_ids, sample_levels) {
  group_matrix <- Matrix::sparseMatrix(
    i = seq_along(cell_indices),
    j = match(sample_ids[cell_indices], sample_levels),
    x = 1,
    dims = c(length(cell_indices), length(sample_levels)),
    dimnames = list(colnames(counts)[cell_indices], sample_levels)
  )

  pseudobulk <- counts[, cell_indices, drop = FALSE] %*% group_matrix
  rownames(pseudobulk) <- rownames(counts)
  colnames(pseudobulk) <- sample_levels
  pseudobulk
}

write_pseudobulk_csv <- function(pseudobulk, output_csv, expected_samples) {
  if (!identical(colnames(pseudobulk), expected_samples)) {
    stop("Internal error: pseudobulk sample names or order changed.")
  }

  output_df <- as.data.frame(as.matrix(pseudobulk), check.names = FALSE)
  output_df <- cbind(
    data.frame(gene = rownames(pseudobulk), check.names = FALSE),
    output_df
  )
  colnames(output_df) <- c("gene", expected_samples)
  utils::write.csv(output_df, output_csv, row.names = FALSE, quote = TRUE)
}

## ----------------------------- Read and Validate ---------------------------
message("Reading Seurat object: ", input_rds)
sc_all <- readRDS(input_rds)

if (!inherits(sc_all, "Seurat")) {
  stop("Input RDS does not contain a Seurat object.")
}
if (!"RNA" %in% Seurat::Assays(sc_all)) {
  stop("Assay 'RNA' not found in the Seurat object.")
}
if (!"orig.ident" %in% colnames(sc_all@meta.data)) {
  stop("Column 'orig.ident' not found in Seurat metadata.")
}
if (!celltype_col %in% colnames(sc_all@meta.data)) {
  stop("Cell-class column not found in Seurat metadata: ", celltype_col)
}

## 兼容 Assay(V4) 与 Assay5(V5)：Layers(search=) 在两种 assay 类上语义不一致
## （V4 字面量 intersect，V5 精确匹配+正则兜底），因此自行 grep 层名。
count_layers <- grep("^counts", SeuratObject::Layers(sc_all[["RNA"]]), value = TRUE)
if (length(count_layers) == 0) {
  stop("No counts layer found in assay 'RNA'.")
}

if (length(count_layers) > 1) {
  message(
    "Joining ", length(count_layers),
    " RNA count layers in memory: ", paste(count_layers, collapse = ", ")
  )
  sc_all <- SeuratObject::JoinLayers(
    sc_all,
    assay = "RNA",
    layers = "counts",
    new = "counts"
  )
  count_layer <- "counts"
} else {
  count_layer <- count_layers[[1]]
}

message("Extracting raw counts from assay 'RNA', layer '", count_layer, "'.")
counts <- tryCatch(
  SeuratObject::LayerData(sc_all, assay = "RNA", layer = count_layer),
  error = function(e) {
    stop("Failed to read RNA counts layer: ", conditionMessage(e))
  }
)

if (nrow(counts) == 0 || ncol(counts) == 0) {
  stop("RNA counts matrix is empty.")
}
if (is.null(rownames(counts)) || anyNA(rownames(counts)) || any(!nzchar(rownames(counts)))) {
  stop("RNA counts matrix contains missing or empty gene names.")
}
if (anyDuplicated(rownames(counts))) {
  stop("RNA counts matrix contains duplicated gene names.")
}
gene_ids <- rownames(counts)
ensembl_ids <- gene_ids[grepl("^ENSG[0-9]+(?:\\.[0-9]+)?$", gene_ids, perl = TRUE)]
entrez_ids <- gene_ids[grepl("^[0-9]+$", gene_ids)]
if (length(ensembl_ids) > 0 || length(entrez_ids) > 0) {
  stop(
    "RNA feature names must be standard gene SYMBOL for the DESeq2/GSEA/ORA pipeline; ",
    "Ensembl or pure numeric ENTREZID values were detected."
  )
}
if (!setequal(colnames(counts), rownames(sc_all@meta.data))) {
  stop("RNA counts cells do not exactly match the Seurat metadata cells.")
}

meta <- sc_all@meta.data[colnames(counts), , drop = FALSE]
sample_ids <- as.character(meta[["orig.ident"]])
celltype_ids <- as.character(meta[[celltype_col]])

if (anyNA(sample_ids) || any(!nzchar(sample_ids)) || any(!nzchar(trimws(sample_ids)))) {
  stop("orig.ident contains NA, empty, or whitespace-only sample identifiers.")
}
if (anyNA(celltype_ids) || any(!nzchar(celltype_ids)) || any(!nzchar(trimws(celltype_ids)))) {
  stop(celltype_col, " contains NA, empty, or whitespace-only cell-class labels.")
}

stored_counts <- if (inherits(counts, "sparseMatrix")) counts@x else as.vector(counts)
if (any(!is.finite(stored_counts))) {
  stop("RNA counts contain non-finite values.")
}
if (any(stored_counts < 0)) {
  stop("RNA counts contain negative values.")
}
if (any(abs(stored_counts - round(stored_counts)) > sqrt(.Machine$double.eps))) {
  stop("RNA counts contain non-integer values; raw counts are required.")
}

sample_levels <- unique(sample_ids)
celltype_levels <- unique(celltype_ids)

message(
  "Input dimensions: ", nrow(counts), " genes x ", ncol(counts), " cells; ",
  length(sample_levels), " samples; ", length(celltype_levels), " cell classes."
)

sample_celltype_counts <- table(
  factor(sample_ids, levels = sample_levels),
  factor(celltype_ids, levels = celltype_levels)
)
empty_combinations <- which(sample_celltype_counts == 0, arr.ind = TRUE)

if (nrow(empty_combinations) > 0) {
  empty_labels <- apply(
    empty_combinations,
    1,
    function(index) {
      paste0(
        "sample='", rownames(sample_celltype_counts)[index[1]],
        "', cell_class='", colnames(sample_celltype_counts)[index[2]], "'"
      )
    }
  )
  message(
    "No cells were found for the following sample x cell-class combinations; ",
    "these samples will be EXCLUDED from the corresponding cell-class pseudobulk tables ",
    "(0-filled columns would make DESeq2 fail on zero-total-count samples):\n  ",
    paste(empty_labels, collapse = "\n  ")
  )
}

safe_celltype_names <- vapply(celltype_levels, sanitize_filename, character(1))
celltype_filenames <- paste0(input_basename, "_pseudobulk_", safe_celltype_names, ".csv")
all_cells_filename <- paste0(input_basename, "_pseudobulk_all_cells.csv")

all_output_filenames <- c(all_cells_filename, celltype_filenames)
if (anyDuplicated(tolower(all_output_filenames))) {
  duplicated_names <- unique(
    all_output_filenames[
      duplicated(tolower(all_output_filenames)) |
        duplicated(tolower(all_output_filenames), fromLast = TRUE)
    ]
  )
  stop(
    "Output filenames collide after replacing filename-unfriendly characters: ",
    paste(duplicated_names, collapse = ", ")
  )
}

renamed_files <- safe_celltype_names != celltype_levels
if (any(renamed_files)) {
  rename_messages <- paste0(
    "'", celltype_levels[renamed_files], "' -> '",
    celltype_filenames[renamed_files], "'"
  )
  message(
    "Replaced filename-unfriendly characters in cell-class filenames:\n  ",
    paste(rename_messages, collapse = "\n  ")
  )
}

## ----------------------------- Pseudobulk ---------------------------------
message("Creating all-cell pseudobulk counts.")
all_pseudobulk <- aggregate_by_sample(
  counts = counts,
  cell_indices = seq_len(ncol(counts)),
  sample_ids = sample_ids,
  sample_levels = sample_levels
)

all_cells_csv <- file.path(output_dir, all_cells_filename)
write_pseudobulk_csv(all_pseudobulk, all_cells_csv, sample_levels)
message("Wrote: ", all_cells_csv)

sum_by_celltype <- Matrix::Matrix(
  0,
  nrow = nrow(all_pseudobulk),
  ncol = ncol(all_pseudobulk),
  sparse = TRUE,
  dimnames = dimnames(all_pseudobulk)
)

for (celltype_index in seq_along(celltype_levels)) {
  celltype_label <- celltype_levels[celltype_index]
  cell_indices <- which(celltype_ids == celltype_label)

  ## 剔除该细胞大类下无细胞的样本（原因见上方空组合说明），
  ## 避免保留全0列导致下游 DESeq2 对零总counts样本报错。
  sample_counts <- table(factor(sample_ids[cell_indices], levels = sample_levels))
  keep_samples <- sample_levels[sample_counts > 0]

  if (length(keep_samples) < length(sample_levels)) {
    message(
      "Cell class '", celltype_label, "': excluding sample(s) ",
      paste(setdiff(sample_levels, keep_samples), collapse = ", "),
      " from its pseudobulk table (no cells for this cell class). ",
      "These samples may remain in group.csv; DESeq2 takes the intersection ",
      "of expression-matrix columns and group.csv samples."
    )
  }
  if (length(keep_samples) == 0) {
    message(
      "Cell class '", celltype_label, "' has no cells at all; skipping its table."
    )
    next
  }

  celltype_pseudobulk <- aggregate_by_sample(
    counts = counts,
    cell_indices = cell_indices,
    sample_ids = sample_ids,
    sample_levels = keep_samples
  )

  celltype_csv <- file.path(output_dir, celltype_filenames[celltype_index])
  write_pseudobulk_csv(celltype_pseudobulk, celltype_csv, keep_samples)
  message(
    "Wrote cell class '", celltype_label, "' (", length(cell_indices),
    " cells, ", length(keep_samples), " samples): ", celltype_csv
  )

  ## 对齐回全部样本水平（被剔除列视为0），用于"各类之和==全细胞表"验证
  celltype_full <- Matrix::Matrix(
    0,
    nrow = nrow(all_pseudobulk),
    ncol = ncol(all_pseudobulk),
    sparse = TRUE,
    dimnames = dimnames(all_pseudobulk)
  )
  celltype_full[, keep_samples] <- celltype_pseudobulk
  sum_by_celltype <- sum_by_celltype + celltype_full
}

if (Matrix::nnzero(all_pseudobulk - sum_by_celltype) != 0) {
  stop(
    "Validation failed: the sum of cell-class pseudobulk matrices does not ",
    "equal the all-cell pseudobulk matrix."
  )
}

message(
  "Validation passed: cell-class matrices sum exactly to the all-cell matrix."
)
message("Done. Output directory: ", output_dir)
message("")
message("Next step (group.csv is an additional required input):")
message(
  "  Rscript DESeq2.r ", shQuote(all_cells_csv),
  " ", shQuote("<group.csv>"),
  " ", shQuote("<numerator_group>"),
  " ", shQuote("<denominator_group>"),
  " ", shQuote("<gene_annotation.csv>"),
  " ", shQuote("<gmt_dir>")
)
message("  The same DESeq2 command applies to each cell-class CSV.")
