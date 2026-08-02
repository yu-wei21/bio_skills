#!/usr/bin/env Rscript

## 输入：Seurat RDS 文件。对象的 metadata 中必须包含由 --sample-col 和
## --celltype-col 指定的两列；两列均不可含 NA 或空字符串。样本和细胞类型
## 的因子水平（若已设置）决定图中与结果表中的排列顺序，否则按首次出现顺序排列。
## 输出：默认写入 output/plot-celltype-barplot/。以输入 RDS 文件名（去除 .rds）
## 为前缀，生成样本内细胞类型比例 CSV、PNG 和 PDF 堆叠柱状图。
## 示例：Rscript plot-celltype-barplot.r --input data/example.rds \\
##   --sample-col sample --celltype-col major_celltype
## 可选：--output-dir output/custom-directory

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(scales)
})

usage <- paste(
  "Usage:",
  "  Rscript plot-celltype-barplot.r --input <input.rds> --sample-col <column> \\",
  "    --celltype-col <column> [--output-dir <directory>]",
  sep = "\n"
)

parse_args <- function(args) {
  allowed <- c("input", "sample-col", "celltype-col", "output-dir")
  values <- setNames(rep(NA_character_, length(allowed)), allowed)

  i <- 1L
  while (i <= length(args)) {
    arg <- args[[i]]
    if (arg %in% c("-h", "--help")) {
      cat(usage, "\n")
      quit(save = "no", status = 0)
    }
    if (!startsWith(arg, "--")) {
      stop("Unexpected argument: ", arg, "\n", usage, call. = FALSE)
    }

    key <- substring(arg, 3L)
    if (!key %in% allowed) {
      stop("Unknown argument: ", arg, "\n", usage, call. = FALSE)
    }
    if (!is.na(values[[key]])) {
      stop("Argument supplied more than once: ", arg, call. = FALSE)
    }
    if (i == length(args) || startsWith(args[[i + 1L]], "--")) {
      stop("Missing value for argument: ", arg, call. = FALSE)
    }

    values[[key]] <- args[[i + 1L]]
    i <- i + 2L
  }

  required <- c("input", "sample-col", "celltype-col")
  missing <- required[is.na(values[required])]
  if (length(missing) > 0L) {
    stop("Missing required argument(s): --", paste(missing, collapse = ", --"),
      "\n", usage,
      call. = FALSE
    )
  }
  values
}

observed_levels <- function(x) {
  observed <- unique(as.character(x))
  if (is.factor(x)) {
    levels(x)[levels(x) %in% observed]
  } else {
    observed
  }
}

args <- parse_args(commandArgs(trailingOnly = TRUE))
input_file <- args[["input"]]
sample_column <- args[["sample-col"]]
celltype_column <- args[["celltype-col"]]
output_dir <- args[["output-dir"]]
if (is.na(output_dir)) {
  output_dir <- file.path("output", "plot-celltype-barplot")
}

if (!file.exists(input_file)) {
  stop("Input file does not exist: ", input_file, call. = FALSE)
}

obj <- readRDS(input_file)
if (!inherits(obj, "Seurat")) {
  stop("Input RDS is not a Seurat object: ", input_file, call. = FALSE)
}

metadata <- obj@meta.data
required_columns <- c(sample_column, celltype_column)
missing_columns <- setdiff(required_columns, colnames(metadata))
if (length(missing_columns) > 0L) {
  stop("Metadata column(s) not found: ", paste(missing_columns, collapse = ", "),
    call. = FALSE
  )
}
if (!identical(rownames(metadata), colnames(obj))) {
  stop("Metadata row names do not exactly match Seurat cell barcodes.", call. = FALSE)
}

sample_values <- metadata[[sample_column]]
celltype_values <- metadata[[celltype_column]]
sample_text <- as.character(sample_values)
celltype_text <- as.character(celltype_values)

missing_sample <- is.na(sample_text) | !nzchar(trimws(sample_text))
missing_celltype <- is.na(celltype_text) | !nzchar(trimws(celltype_text))
if (any(missing_sample) || any(missing_celltype)) {
  stop(
    "Missing or empty metadata values: ",
    sample_column, " = ", sum(missing_sample), "; ",
    celltype_column, " = ", sum(missing_celltype),
    call. = FALSE
  )
}

sample_levels <- observed_levels(sample_values)
celltype_levels <- observed_levels(celltype_values)
if (length(sample_levels) == 0L || length(celltype_levels) == 0L) {
  stop("No observed samples or cell types remain after validation.", call. = FALSE)
}

count_matrix <- table(
  factor(sample_text, levels = sample_levels),
  factor(celltype_text, levels = celltype_levels)
)
if (sum(count_matrix) != nrow(metadata)) {
  stop("Cell counts do not equal the number of metadata rows.", call. = FALSE)
}

proportion_df <- expand.grid(
  sample = sample_levels,
  celltype = celltype_levels,
  KEEP.OUT.ATTRS = FALSE,
  stringsAsFactors = FALSE
)
proportion_df$n_cells <- as.integer(count_matrix)
sample_totals <- rowSums(count_matrix)
proportion_df$n_cells_in_sample <- unname(sample_totals[proportion_df$sample])
proportion_df$proportion <- proportion_df$n_cells / proportion_df$n_cells_in_sample

if (any(!is.finite(proportion_df$proportion))) {
  stop("Could not calculate finite within-sample proportions.", call. = FALSE)
}
if (any(abs(tapply(proportion_df$proportion, proportion_df$sample, sum) - 1) > 1e-12)) {
  stop("Within-sample proportions do not sum to 1.", call. = FALSE)
}

proportion_df$sample <- factor(proportion_df$sample, levels = sample_levels)
proportion_df$celltype <- factor(proportion_df$celltype, levels = celltype_levels)

celltype_colors <- setNames(
  hue_pal(l = 65, c = 100)(length(celltype_levels)),
  celltype_levels
)
plot_width <- max(8, 0.38 * length(sample_levels) + 3)

plot <- ggplot(proportion_df, aes(x = sample, y = proportion, fill = celltype)) +
  geom_col(width = 0.85) +
  scale_fill_manual(values = celltype_colors, name = "Cell type") +
  scale_y_continuous(labels = label_percent(accuracy = 1), expand = c(0, 0)) +
  labs(x = "Sample", y = "Cell proportion") +
  theme_classic(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    legend.position = "right"
  )

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
input_prefix <- tools::file_path_sans_ext(basename(input_file))
output_csv <- file.path(output_dir, paste0(input_prefix, "_celltype_proportions.csv"))
output_png <- file.path(output_dir, paste0(input_prefix, "_celltype_proportions.png"))
output_pdf <- file.path(output_dir, paste0(input_prefix, "_celltype_proportions.pdf"))

write.csv(proportion_df, output_csv, row.names = FALSE)
ggsave(output_png, plot, width = plot_width, height = 7, dpi = 300, bg = "white")
ggsave(output_pdf, plot, width = plot_width, height = 7, bg = "white")

message("Wrote: ", output_csv)
message("Wrote: ", output_png)
message("Wrote: ", output_pdf)
