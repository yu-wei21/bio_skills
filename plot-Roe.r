#!/usr/bin/env Rscript

# Pooled-cell Ro/e analysis
#
# Usage:
#   Rscript Roe.r <counts.csv> [output_dir]
#
# Input CSV:
#   - rows: cell types
#   - first column: unique cell type symbols/labels
#   - remaining columns: study groups
#   - values: non-negative integer cell counts
#
# Outputs:
#   <input_stem>.roe.xlsx
#   <input_stem>.roe.heatmap.pdf
#   <input_stem>.roe.heatmap.png

required_packages <- c("ComplexHeatmap", "openxlsx", "ragg")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_packages) > 0) {
  stop(
    "Missing required R package(s): ",
    paste(missing_packages, collapse = ", "),
    ". Please install them before running this script.",
    call. = FALSE
  )
}

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1 || length(args) > 2) {
  stop(
    "Usage: Rscript Roe.r <counts.csv> [output_dir]",
    call. = FALSE
  )
}

input_file <- normalizePath(args[1], mustWork = TRUE)
output_dir <- if (length(args) == 2) args[2] else dirname(input_file)

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
}
output_dir <- normalizePath(output_dir, mustWork = TRUE)

input_data <- read.csv(
  input_file,
  header = TRUE,
  check.names = FALSE,
  stringsAsFactors = FALSE,
  na.strings = c("NA", "NaN", "")
)

if (nrow(input_data) < 2) {
  stop("The input must contain at least two cell types (rows).", call. = FALSE)
}
if (ncol(input_data) < 3) {
  stop(
    "The input must contain one cell-type column and at least two group columns.",
    call. = FALSE
  )
}

cell_types <- as.character(input_data[[1]])
group_names <- colnames(input_data)[-1]

if (anyNA(cell_types) || any(trimws(cell_types) == "")) {
  stop("Cell type labels in the first column cannot be missing or empty.", call. = FALSE)
}
if (any(cell_types != trimws(cell_types))) {
  stop("Cell type labels cannot contain leading or trailing whitespace.", call. = FALSE)
}
if (anyDuplicated(cell_types)) {
  duplicated_names <- unique(cell_types[duplicated(cell_types)])
  stop(
    "Cell type labels must be unique. Duplicated label(s): ",
    paste(duplicated_names, collapse = ", "),
    call. = FALSE
  )
}

if (anyNA(group_names) || any(trimws(group_names) == "")) {
  stop("Group column names cannot be missing or empty.", call. = FALSE)
}
if (any(group_names != trimws(group_names))) {
  stop("Group column names cannot contain leading or trailing whitespace.", call. = FALSE)
}
if (anyDuplicated(group_names)) {
  duplicated_names <- unique(group_names[duplicated(group_names)])
  stop(
    "Group column names must be unique. Duplicated name(s): ",
    paste(duplicated_names, collapse = ", "),
    call. = FALSE
  )
}

count_columns <- input_data[-1]
counts_list <- lapply(seq_along(count_columns), function(index) {
  original_values <- count_columns[[index]]
  numeric_values <- suppressWarnings(as.numeric(as.character(original_values)))

  invalid_numeric <- is.na(numeric_values) & !is.na(original_values)
  if (any(invalid_numeric)) {
    bad_rows <- cell_types[invalid_numeric]
    stop(
      "Non-numeric count(s) found in group '", group_names[index],
      "' for cell type(s): ", paste(bad_rows, collapse = ", "),
      call. = FALSE
    )
  }

  numeric_values
})

counts <- do.call(cbind, counts_list)
dimnames(counts) <- list(cell_types, group_names)

if (anyNA(counts) || any(!is.finite(counts))) {
  stop("Cell counts cannot contain NA, NaN, or infinite values.", call. = FALSE)
}
if (any(counts < 0)) {
  stop("Cell counts must be non-negative.", call. = FALSE)
}
if (any(abs(counts - round(counts)) > sqrt(.Machine$double.eps))) {
  stop("Cell counts must be integers.", call. = FALSE)
}
counts <- round(counts)
storage.mode(counts) <- "numeric"

row_totals <- rowSums(counts)
column_totals <- colSums(counts)

if (any(row_totals == 0)) {
  stop(
    "All-zero cell type row(s) are not allowed: ",
    paste(names(row_totals)[row_totals == 0], collapse = ", "),
    call. = FALSE
  )
}
if (any(column_totals == 0)) {
  stop(
    "All-zero group column(s) are not allowed: ",
    paste(names(column_totals)[column_totals == 0], collapse = ", "),
    call. = FALSE
  )
}

grand_total <- sum(counts)
expected <- outer(row_totals, column_totals) / grand_total
dimnames(expected) <- dimnames(counts)

roe <- counts / expected
dimnames(roe) <- dimnames(counts)

if (any(!is.finite(expected)) || any(!is.finite(roe))) {
  stop("Non-finite expected or Ro/e values were generated.", call. = FALSE)
}

input_stem <- tools::file_path_sans_ext(basename(input_file))
output_prefix <- file.path(output_dir, paste0(input_stem, ".roe"))
excel_file <- paste0(output_prefix, ".xlsx")
pdf_file <- paste0(output_prefix, ".heatmap.pdf")
png_file <- paste0(output_prefix, ".heatmap.png")

matrix_to_table <- function(matrix_data) {
  data.frame(
    celltype = rownames(matrix_data),
    matrix_data,
    check.names = FALSE,
    row.names = NULL
  )
}

workbook <- openxlsx::createWorkbook()
sheet_data <- list(
  observed = matrix_to_table(counts),
  expected = matrix_to_table(expected),
  roe = matrix_to_table(roe)
)

for (sheet_name in names(sheet_data)) {
  openxlsx::addWorksheet(workbook, sheet_name)
  openxlsx::writeData(
    workbook,
    sheet = sheet_name,
    x = sheet_data[[sheet_name]],
    keepNA = TRUE
  )
  openxlsx::freezePane(
    workbook,
    sheet = sheet_name,
    firstRow = TRUE,
    firstCol = TRUE
  )
  openxlsx::setColWidths(
    workbook,
    sheet = sheet_name,
    cols = seq_len(ncol(sheet_data[[sheet_name]])),
    widths = "auto"
  )
}

openxlsx::saveWorkbook(workbook, excel_file, overwrite = TRUE)

# Match the reference pooled-cell workflow: assign each cell type to the group
# with its largest Ro/e value, then place later input groups first. Ties are
# resolved in favor of the first group, matching which.max().
preferred_group <- max.col(roe, ties.method = "first")
row_order <- order(preferred_group, decreasing = TRUE, method = "radix")
plot_roe_raw <- roe[row_order, , drop = FALSE]
plot_roe <- pmin(plot_roe_raw, 3)

preference_colors <- c("#FEE6CE", "#FDC08C", "#F5904B", "#E6550D")

preference_class <- function(value) {
  if (value < 1) {
    1L
  } else if (value < 1.5) {
    2L
  } else if (value < 3) {
    3L
  } else {
    4L
  }
}

preference_symbol <- function(value) {
  if (value < 1) {
    "±"
  } else if (value < 1.5) {
    "+"
  } else if (value < 3) {
    "++"
  } else {
    "+++"
  }
}

cell_fun <- function(j, i, x, y, width, height, fill) {
  value <- plot_roe_raw[i, j]
  category <- preference_class(value)
  grid::grid.rect(
    x = x,
    y = y,
    width = width,
    height = height,
    gp = grid::gpar(
      fill = preference_colors[category],
      col = "white",
      lwd = 1
    )
  )
  grid::grid.text(
    preference_symbol(value),
    x = x,
    y = y,
    gp = grid::gpar(
      col = if (category == 4L) "white" else "black",
      fontsize = 10,
      fontface = "bold"
    )
  )
}

heatmap <- ComplexHeatmap::Heatmap(
  plot_roe,
  name = "Ro/e",
  col = grDevices::colorRampPalette(preference_colors)(100),
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  show_heatmap_legend = FALSE,
  rect_gp = grid::gpar(type = "none"),
  cell_fun = cell_fun,
  row_names_side = "left",
  column_names_rot = 45,
  row_names_gp = grid::gpar(fontsize = 10),
  column_names_gp = grid::gpar(fontsize = 10),
  column_title = "Cell Type Preference (pooled-cell Ro/e)",
  column_title_gp = grid::gpar(fontsize = 13, fontface = "bold")
)

preference_legend <- ComplexHeatmap::Legend(
  title = "Ro/e",
  labels = c("<1", "1–<1.5", "1.5–<3", "≥3"),
  legend_gp = grid::gpar(fill = preference_colors, col = "white"),
  grid_height = grid::unit(5, "mm"),
  grid_width = grid::unit(5, "mm")
)

draw_heatmap <- function() {
  ComplexHeatmap::draw(
    heatmap,
    heatmap_legend_list = list(preference_legend),
    heatmap_legend_side = "right"
  )
}

plot_width <- max(6, 2.8 + 0.65 * ncol(plot_roe))
plot_height <- max(5, 2.5 + 0.38 * nrow(plot_roe))

if (capabilities("cairo")) {
  grDevices::cairo_pdf(pdf_file, width = plot_width, height = plot_height)
} else {
  grDevices::pdf(pdf_file, width = plot_width, height = plot_height)
}
draw_heatmap()
grDevices::dev.off()

ragg::agg_png(
  png_file,
  width = plot_width,
  height = plot_height,
  units = "in",
  res = 150,
  background = "white"
)
draw_heatmap()
grDevices::dev.off()

message("Ro/e analysis completed.")
message("Input orientation: rows = cell types; columns = groups.")
message("Cell types: ", nrow(counts), "; groups: ", ncol(counts), "; cells: ", grand_total)
message("Excel: ", excel_file)
message("PDF: ", pdf_file)
message("PNG: ", png_file)
