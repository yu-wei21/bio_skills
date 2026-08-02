#!/usr/bin/env Rscript

# 输入：Seurat RDS；marker CSV/TSV 至少含 panel、marker 两列，所有 marker 必须在指定
#       assay 中存在，分组列由 --group-by 指定。
# 输出：以 --output-prefix 为前缀生成 marker DotPlot PDF/PNG 与作图数据 CSV；默认前缀为
#       output/marker_dotplot。
# 示例：Rscript plot_marker_dotplot.r --seurat-rds input/seurat.rds --markers marker_panels.csv
#
# Seurat marker DotPlot：命令行版本
#
# marker 文件为 CSV/TSV，至少含两列：panel,marker。
# 文件行顺序决定分面顺序及各分面内 marker 顺序；可选 panel_order、marker_order
# 两列进行显式排序。每个 marker 只能出现一次。
#
# 最小示例：同上。
# 可选参数 --group-by、--assay、--output-prefix 未提供时，分别默认为
# celltype、RNA、output/marker_dotplot。
#
# 查看全部参数：Rscript plot_marker_dotplot.r --help

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(scales)
})

usage <- function() {
  paste(
    "Usage:",
    "  Rscript plot_marker_dotplot.r --seurat-rds FILE --markers FILE [options]",
    "",
    "Required:",
    "  --seurat-rds FILE       Seurat .rds file.",
    "  --markers FILE          CSV/TSV with columns: panel,marker.",
    "",
    "Optional:",
    "  --assay NAME            Assay to use (default: RNA).",
    "  --group-by COLUMN       Metadata column used for rows (default: celltype).",
    "  --output-prefix PREFIX  Output prefix (default: output/marker_dotplot).",
    "  --group-order A,B,C     Comma-separated top-to-bottom row order.",
    "  --expression-min NUM    Lower color limit (default: -2).",
    "  --expression-max NUM    Upper color limit (default: 2).",
    "  --point-max-size NUM    Maximum dot diameter (default: 4.8).",
    "  --width NUM             Figure width in inches; inferred if omitted.",
    "  --height NUM            Figure height in inches; inferred if omitted.",
    "  --help                  Show this help message.",
    sep = "\n"
  )
}

assert_that <- function(condition, message) {
  if (!isTRUE(condition)) stop(message, call. = FALSE)
}

parse_cli_args <- function(args) {
  defaults <- list(
    seurat_rds = NULL,
    markers = NULL,
    output_prefix = "output/marker_dotplot",
    assay = "RNA",
    group_by = "celltype",
    group_order = NULL,
    expression_min = -2,
    expression_max = 2,
    point_max_size = 4.8,
    width = NULL,
    height = NULL
  )
  allowed <- names(defaults)
  i <- 1L

  while (i <= length(args)) {
    arg <- args[[i]]
    if (arg %in% c("--help", "-h")) {
      cat(usage(), "\n")
      quit(status = 0)
    }
    assert_that(startsWith(arg, "--"), paste0("Invalid argument: ", arg, "\n\n", usage()))

    if (grepl("=", arg, fixed = TRUE)) {
      split_at <- regexpr("=", arg, fixed = TRUE)[1]
      key <- substring(arg, 3L, split_at - 1L)
      value <- substring(arg, split_at + 1L)
    } else {
      key <- substring(arg, 3L)
      assert_that(i < length(args), paste0("Missing value for --", key))
      i <- i + 1L
      value <- args[[i]]
    }

    key <- gsub("-", "_", key, fixed = TRUE)
    assert_that(key %in% allowed, paste0("Unknown argument: --", key, "\n\n", usage()))
    assert_that(nzchar(value), paste0("Empty value for --", key))
    defaults[[key]] <- value
    i <- i + 1L
  }

  required <- c("seurat_rds", "markers")
  missing <- required[vapply(defaults[required], is.null, logical(1))]
  assert_that(length(missing) == 0,
              paste0("Missing required argument(s): --", gsub("_", "-", missing), "\n\n", usage()))

  numeric_args <- c("expression_min", "expression_max", "point_max_size", "width", "height")
  for (name in numeric_args) {
    if (!is.null(defaults[[name]])) {
      defaults[[name]] <- suppressWarnings(as.numeric(defaults[[name]]))
      assert_that(length(defaults[[name]]) == 1 && is.finite(defaults[[name]]),
                  paste0("--", gsub("_", "-", name), " must be a finite number."))
    }
  }
  assert_that(defaults$expression_min < defaults$expression_max,
              "--expression-min must be smaller than --expression-max.")
  assert_that(defaults$point_max_size > 0, "--point-max-size must be greater than zero.")
  for (name in c("width", "height")) {
    if (!is.null(defaults[[name]])) {
      assert_that(defaults[[name]] > 0, paste0("--", name, " must be greater than zero."))
    }
  }

  if (!is.null(defaults$group_order)) {
    defaults$group_order <- trimws(unlist(strsplit(defaults$group_order, ",", fixed = TRUE)))
    defaults$group_order <- defaults$group_order[nzchar(defaults$group_order)]
    assert_that(length(defaults$group_order) > 0,
                "--group-order must contain at least one non-empty group name.")
  }
  defaults
}

read_marker_panels <- function(path) {
  assert_that(file.exists(path), paste0("Marker file not found: ", path))
  extension <- tolower(tools::file_ext(path))
  assert_that(extension %in% c("csv", "tsv", "txt"),
              "Marker file must be CSV, TSV, or TXT.")
  separator <- if (extension == "csv") "," else "\t"
  marker_df <- read.table(
    path, header = TRUE, sep = separator, quote = "\"", comment.char = "",
    stringsAsFactors = FALSE, check.names = FALSE, na.strings = c("", "NA")
  )
  required_columns <- c("panel", "marker")
  missing_columns <- setdiff(required_columns, colnames(marker_df))
  assert_that(length(missing_columns) == 0,
              paste0("Marker file is missing column(s): ", paste(missing_columns, collapse = ", ")))

  marker_df$panel <- trimws(as.character(marker_df$panel))
  marker_df$marker <- trimws(as.character(marker_df$marker))
  assert_that(nrow(marker_df) > 0 && all(nzchar(marker_df$panel)) && all(nzchar(marker_df$marker)),
              "Columns panel and marker must not contain empty values.")
  duplicated_markers <- unique(marker_df$marker[duplicated(marker_df$marker)])
  assert_that(length(duplicated_markers) == 0,
              paste0("A marker can belong to only one panel. Duplicated marker(s): ",
                     paste(duplicated_markers, collapse = ", ")))

  marker_df$.input_order <- seq_len(nrow(marker_df))
  if ("panel_order" %in% colnames(marker_df)) {
    marker_df$.panel_order <- suppressWarnings(as.numeric(marker_df$panel_order))
    assert_that(all(is.finite(marker_df$.panel_order)), "panel_order must contain finite numbers.")
  } else {
    panel_names <- unique(marker_df$panel)
    marker_df$.panel_order <- match(marker_df$panel, panel_names)
  }
  if ("marker_order" %in% colnames(marker_df)) {
    marker_df$.marker_order <- suppressWarnings(as.numeric(marker_df$marker_order))
    assert_that(all(is.finite(marker_df$.marker_order)), "marker_order must contain finite numbers.")
  } else {
    marker_df$.marker_order <- marker_df$.input_order
  }

  marker_df <- marker_df[order(
    marker_df$.panel_order, marker_df$.marker_order, marker_df$.input_order
  ), , drop = FALSE]
  panel_names <- unique(marker_df$panel)
  panels <- lapply(panel_names, function(panel) marker_df$marker[marker_df$panel == panel])
  names(panels) <- panel_names
  panels
}

validate_configuration <- function(obj, assay, group_column, groups, panels) {
  assert_that(assay %in% Assays(obj), paste0("Assay not found: ", assay))
  assert_that(group_column %in% colnames(obj@meta.data),
              paste0("Metadata column not found: ", group_column))

  markers <- unlist(panels, use.names = FALSE)
  missing_markers <- setdiff(markers, rownames(obj[[assay]]))
  assert_that(length(missing_markers) == 0,
              paste0("Marker(s) not found in assay '", assay, "': ",
                     paste(missing_markers, collapse = ", ")))

  absent_groups <- setdiff(groups, unique(as.character(obj@meta.data[[group_column]])))
  assert_that(length(absent_groups) == 0,
              paste0("--group-order contains group(s) absent from '", group_column, "': ",
                     paste(absent_groups, collapse = ", ")))
}

get_group_order <- function(obj, group_column, configured_order) {
  observed <- as.character(obj@meta.data[[group_column]])
  observed <- observed[!is.na(observed) & nzchar(observed)]
  assert_that(length(observed) > 0,
              paste0("No non-missing values found in metadata column: ", group_column))
  if (!is.null(configured_order)) return(configured_order)

  original <- obj@meta.data[[group_column]]
  if (is.factor(original)) return(levels(droplevels(original)))
  sort(unique(observed))
}

make_dotplot <- function(obj, assay, group_column, groups, panels, expression_range, max_size) {
  DefaultAssay(obj) <- assay
  obj$dotplot_group <- factor(as.character(obj@meta.data[[group_column]]), levels = groups)
  Idents(obj) <- "dotplot_group"

  marker_order <- unlist(panels, use.names = FALSE)
  marker_annotation <- data.frame(
    marker = marker_order,
    marker_panel = rep(names(panels), lengths(panels)),
    stringsAsFactors = FALSE
  )
  dot_data <- DotPlot(
    object = obj, features = marker_order, group.by = "dotplot_group",
    col.min = expression_range[1], col.max = expression_range[2]
  )$data

  dot_data$marker <- as.character(dot_data$features.plot)
  dot_data$group <- as.character(dot_data$id)
  dot_data$marker_panel <- marker_annotation$marker_panel[
    match(dot_data$marker, marker_annotation$marker)
  ]
  dot_data$avg.exp.scaled <- pmax(
    pmin(dot_data$avg.exp.scaled, expression_range[2]), expression_range[1]
  )
  dot_data$marker <- factor(dot_data$marker, levels = marker_order)
  dot_data$group <- factor(dot_data$group, levels = groups)
  dot_data$marker_panel <- factor(dot_data$marker_panel, levels = names(panels))

  plot <- ggplot(dot_data, aes(x = marker, y = group)) +
    geom_point(aes(size = pct.exp, color = avg.exp.scaled), alpha = 0.95) +
    facet_grid(. ~ marker_panel, scales = "free_x", space = "free_x") +
    scale_y_discrete(limits = rev(groups), drop = FALSE) +
    scale_size_area(
      name = "Percent expressed", limits = c(0, 100),
      breaks = c(0, 25, 50, 75, 100), max_size = max_size
    ) +
    scale_color_gradient2(
      name = "Average expression", low = "#3FA7D6", mid = "white", high = "#C05A4C",
      midpoint = 0, limits = expression_range, oob = squish,
      breaks = seq(expression_range[1], expression_range[2], by = 1)
    ) +
    labs(x = NULL, y = NULL) +
    theme_bw(base_size = 10) +
    theme(
      panel.grid.major = element_line(color = "grey86", linewidth = 0.25),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(color = "grey55", fill = NA, linewidth = 0.35),
      strip.background = element_rect(fill = "grey82", color = "grey55", linewidth = 0.35),
      strip.text.x = element_text(size = 9, face = "bold"),
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 10),
      axis.text.y = element_text(size = 10, color = "black"),
      axis.ticks = element_line(color = "grey55", linewidth = 0.25),
      legend.position = "right", legend.title = element_text(size = 9),
      legend.text = element_text(size = 9), plot.margin = margin(4, 8, 4, 4)
    )
  list(plot = plot, data = dot_data)
}

suggest_figure_size <- function(n_markers, n_groups) {
  c(
    width = min(max(5.5 + 0.34 * n_markers, 8), 22),
    height = min(max(1.5 + 0.30 * n_groups, 3.5), 10)
  )
}

options <- parse_cli_args(commandArgs(trailingOnly = TRUE))
assert_that(file.exists(options$seurat_rds), paste0("Seurat RDS file not found: ", options$seurat_rds))

marker_panels <- read_marker_panels(options$markers)
obj <- readRDS(options$seurat_rds)
group_levels <- get_group_order(obj, options$group_by, options$group_order)
validate_configuration(obj, options$assay, options$group_by, group_levels, marker_panels)

result <- make_dotplot(
  obj = obj, assay = options$assay, group_column = options$group_by,
  groups = group_levels, panels = marker_panels,
  expression_range = c(options$expression_min, options$expression_max),
  max_size = options$point_max_size
)

output_dir <- dirname(options$output_prefix)
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
figure_size <- suggest_figure_size(
  n_markers = length(unlist(marker_panels, use.names = FALSE)), n_groups = length(group_levels)
)
if (!is.null(options$width)) figure_size[["width"]] <- options$width
if (!is.null(options$height)) figure_size[["height"]] <- options$height

ggsave(
  filename = paste0(options$output_prefix, ".pdf"), plot = result$plot,
  width = figure_size[["width"]], height = figure_size[["height"]],
  units = "in", limitsize = FALSE
)
ggsave(
  filename = paste0(options$output_prefix, ".png"), plot = result$plot,
  width = figure_size[["width"]], height = figure_size[["height"]],
  units = "in", dpi = 300, limitsize = FALSE
)
write.csv(result$data, paste0(options$output_prefix, "_data.csv"), row.names = FALSE)

message("Saved: ", paste0(options$output_prefix, ".pdf"))
message("Saved: ", paste0(options$output_prefix, ".png"))
message("Saved: ", paste0(options$output_prefix, "_data.csv"))
