#!/usr/bin/env Rscript

## 输入为 CSV 或制表符分隔的文本文件：每行一个样本；样本列和分组列由参数指定；
## 其余列均视为细胞类型丰度。细胞类型在图中的顺序与输入列顺序一致，分组顺序与
## 输入文件中首次出现的顺序一致。
## 输出：默认写入 output/plot-celltype-boxplot/；生成以输入文件名为前缀的长表 CSV、
##       分组细胞类型丰度箱线图 PNG 与 PDF。
##
## 示例：
## Rscript plot-celltype-boxplot.r \
##   --input celltype_abundance.csv \
##   --sample-col sample_id \
##   --group-col group \
##   --output-dir output/celltype-boxplot

suppressPackageStartupMessages({
  library(ggplot2)
  library(scales)
})

usage <- paste(
  "Usage:",
  "  Rscript plot-celltype-boxplot.r --input <table.csv|table.tsv> \\",
  "    --sample-col <column> --group-col <column> [--output-dir <directory>] \\",
  "    [--title <plot title>] [--y-label <axis label>]",
  "",
  "The input must contain one row per sample. All columns other than --sample-col",
  "and --group-col must be non-negative numeric cell-type abundance columns.",
  sep = "\n"
)

parse_args <- function(args) {
  allowed <- c("input", "sample-col", "group-col", "output-dir", "title", "y-label")
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

  required <- c("input", "sample-col", "group-col")
  missing <- required[is.na(values[required])]
  if (length(missing) > 0L) {
    stop(
      "Missing required argument(s): --", paste(missing, collapse = ", --"),
      "\n", usage,
      call. = FALSE
    )
  }
  values
}

read_input_table <- function(input_file) {
  extension <- tolower(tools::file_ext(input_file))
  separator <- if (extension == "csv") "," else if (extension %in% c("tsv", "txt")) "\t" else NA_character_
  if (is.na(separator)) {
    stop("Input must have a .csv, .tsv, or .txt extension: ", input_file, call. = FALSE)
  }

  read.table(
    input_file,
    header = TRUE,
    sep = separator,
    quote = "\"",
    comment.char = "",
    check.names = FALSE,
    stringsAsFactors = FALSE,
    na.strings = c("", "NA", "NaN")
  )
}

args <- parse_args(commandArgs(trailingOnly = TRUE))
input_file <- args[["input"]]
sample_column <- args[["sample-col"]]
group_column <- args[["group-col"]]
output_dir <- args[["output-dir"]]
if (is.na(output_dir)) {
  output_dir <- file.path("output", "plot-celltype-boxplot")
}
plot_title <- args[["title"]]
if (is.na(plot_title)) {
  plot_title <- "Cell-type abundance by group"
}
y_label <- args[["y-label"]]
if (is.na(y_label)) {
  y_label <- "Cell abundance"
}

if (!file.exists(input_file)) {
  stop("Input file does not exist: ", input_file, call. = FALSE)
}

abundance_wide <- read_input_table(input_file)
if (nrow(abundance_wide) == 0L) {
  stop("Input contains no sample rows.", call. = FALSE)
}
if (anyDuplicated(colnames(abundance_wide))) {
  stop("Input contains duplicated column names.", call. = FALSE)
}

required_columns <- c(sample_column, group_column)
missing_columns <- setdiff(required_columns, colnames(abundance_wide))
if (length(missing_columns) > 0L) {
  stop(
    "Column(s) not found: ", paste(missing_columns, collapse = ", "),
    call. = FALSE
  )
}
if (identical(sample_column, group_column)) {
  stop("--sample-col and --group-col must name different columns.", call. = FALSE)
}

sample_id <- trimws(as.character(abundance_wide[[sample_column]]))
group <- trimws(as.character(abundance_wide[[group_column]]))
if (anyNA(sample_id) || any(!nzchar(sample_id))) {
  stop("Sample column contains missing or empty values: ", sample_column, call. = FALSE)
}
if (anyDuplicated(sample_id)) {
  stop("Sample IDs must be unique; duplicated IDs: ",
       paste(unique(sample_id[duplicated(sample_id)]), collapse = ", "), call. = FALSE)
}
if (anyNA(group) || any(!nzchar(group))) {
  stop("Group column contains missing or empty values: ", group_column, call. = FALSE)
}

celltype_columns <- setdiff(colnames(abundance_wide), required_columns)
if (length(celltype_columns) == 0L) {
  stop("No cell-type abundance columns remain after excluding sample and group columns.", call. = FALSE)
}

non_numeric <- celltype_columns[!vapply(
  abundance_wide[celltype_columns], is.numeric, logical(1)
)]
if (length(non_numeric) > 0L) {
  stop(
    "Cell-type abundance column(s) must be numeric: ",
    paste(non_numeric, collapse = ", "),
    call. = FALSE
  )
}

abundance_matrix <- as.matrix(abundance_wide[celltype_columns])
if (any(!is.finite(abundance_matrix))) {
  stop("Cell-type abundance columns contain missing or non-finite values.", call. = FALSE)
}
if (any(abundance_matrix < 0)) {
  stop("Cell-type abundance values must be non-negative.", call. = FALSE)
}

group_levels <- unique(group)
long_data <- data.frame(
  sample_id = rep(sample_id, times = length(celltype_columns)),
  group = rep(group, times = length(celltype_columns)),
  celltype = rep(celltype_columns, each = length(sample_id)),
  abundance = unlist(abundance_wide[celltype_columns], use.names = FALSE),
  stringsAsFactors = FALSE
)
long_data$group <- factor(long_data$group, levels = group_levels)
long_data$celltype <- factor(long_data$celltype, levels = celltype_columns)

group_colors <- setNames(hue_pal(l = 65, c = 100)(length(group_levels)), group_levels)
plot_width <- max(8, 0.7 * length(celltype_columns) + 3)

plot <- ggplot(long_data, aes(x = celltype, y = abundance, fill = group, color = group)) +
  geom_boxplot(
    position = position_dodge(width = 0.78),
    width = 0.64,
    outlier.shape = NA,
    alpha = 0.75,
    linewidth = 0.45
  ) +
  geom_point(
    position = position_jitterdodge(
      jitter.width = 0.12,
      dodge.width = 0.78,
      seed = 1
    ),
    size = 1.5,
    alpha = 0.70
  ) +
  scale_fill_manual(values = group_colors, name = group_column) +
  scale_color_manual(values = group_colors, guide = "none") +
  labs(title = plot_title, x = "Cell type", y = y_label) +
  theme_classic(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.position = "top",
    legend.title = element_text(face = "bold"),
    plot.margin = margin(8, 12, 8, 8)
  )

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
input_prefix <- tools::file_path_sans_ext(basename(input_file))
output_csv <- file.path(output_dir, paste0(input_prefix, "_celltype_abundance_long.csv"))
output_png <- file.path(output_dir, paste0(input_prefix, "_celltype_abundance_boxplot.png"))
output_pdf <- file.path(output_dir, paste0(input_prefix, "_celltype_abundance_boxplot.pdf"))

write.csv(long_data, output_csv, row.names = FALSE)
ggsave(output_png, plot, width = plot_width, height = 7, units = "in", dpi = 300, bg = "white")
ggsave(output_pdf, plot, width = plot_width, height = 7, units = "in", bg = "white")

message("Samples: ", length(sample_id))
message("Groups: ", paste(group_levels, collapse = ", "))
message("Cell types: ", length(celltype_columns))
message("Wrote: ", output_csv)
message("Wrote: ", output_png)
message("Wrote: ", output_pdf)
