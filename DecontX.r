## 输入：含单一gene×cell原始UMI counts层的Seurat RDS；可选标准10X raw目录或含Seurat/SCE/count矩阵的背景RDS，基因ID须一致且包含人类T/B/Mono/NK marker symbol用于QC。
## 输出：output/DecontX/下以输入RDS文件名开头的校正RDS、metadata CSV、QC PDF和sessionInfo TXT。
## 示例：Rscript DecontX.r sample.rds --background=sample/outs/raw_feature_bc_matrix

## 推荐按单个样本/单个 10X channel 分别运行，不建议直接对合并对象运行。

## DecontX 运行阶段:
## 1. Cell Ranger/STARsolo 等完成 cell calling 之后
## 2. Seurat 正式归一化、整合、聚类、注释和差异分析之前

suppressPackageStartupMessages({
    library(Seurat)
    library(SingleCellExperiment)
    library(decontX)
    library(patchwork)
    library(Matrix)
})


usage <- paste(
    "Usage:",
    "  Rscript DecontX.r input_seurat.rds",
    "    [--background=raw_10x_dir_or_rds]",
    "    [--assay=RNA]",
    "    [--seed=12345]",
    sep = "\n"
)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
    stop(usage, call. = FALSE)
}

opts <- list(
    input_rds = args[1],
    background = NULL,
    assay = "RNA",
    seed = 12345
)

extra_args <- if (length(args) > 1) args[-1] else character(0)
for (arg in extra_args) {
    if (startsWith(arg, "--background=")) {
        opts$background <- sub("^--background=", "", arg)
    } else if (startsWith(arg, "--assay=")) {
        opts$assay <- sub("^--assay=", "", arg)
    } else if (startsWith(arg, "--seed=")) {
        opts$seed <- as.integer(sub("^--seed=", "", arg))
    } else if (is.null(opts$background)) {
        opts$background <- arg
    } else {
        stop("无法识别参数: ", arg, "\n\n", usage, call. = FALSE)
    }
}

null_like <- function(x) {
    is.null(x) || x %in% c("", "NA", "NULL", "null", "None", "none")
}

if (null_like(opts$background)) {
    opts$background <- NULL
}
if (is.na(opts$seed)) {
    stop("--seed 必须是整数。", call. = FALSE)
}

output_dir <- file.path("output", "DecontX")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
input_tag <- tools::file_path_sans_ext(basename(opts$input_rds))
output_prefix <- file.path(output_dir, input_tag)
output_rds <- paste0(output_prefix, "_decontX.rds")
output_cell_metadata <- paste0(output_prefix, "_decontX_cell_metadata.csv")
output_qc_pdf <- paste0(output_prefix, "_decontX_qc.pdf")
output_session <- paste0(output_prefix, "_decontX_sessionInfo.txt")

get_counts_from_seurat <- function(seu, assay = "RNA") {
    ## 显式 SeuratObject::Assays，避免被 decontX 间接加载的 SummarizedExperiment::Assays 遮蔽
    if (!assay %in% SeuratObject::Assays(seu)) {
        stop("Seurat 对象中未找到 assay: ", assay, call. = FALSE)
    }

    assay_obj <- seu[[assay]]
    if (inherits(assay_obj, "Assay5")) {
        count_layers <- Layers(assay_obj, search = "^counts$")
        if (length(count_layers) != 1) {
            stop(
                assay, " assay 中应有且仅有一个 counts layer；当前找到 ",
                length(count_layers), " 个。若是多样本合并对象，建议按样本拆开后分别运行 DecontX。",
                call. = FALSE
            )
        }
        counts <- LayerData(seu, assay = assay, layer = count_layers)
    } else {
        counts <- GetAssayData(seu, assay = assay, slot = "counts")
    }

    if (nrow(counts) == 0 || ncol(counts) == 0) {
        stop("counts 矩阵为空。", call. = FALSE)
    }
    if (is.null(rownames(counts)) || is.null(colnames(counts))) {
        stop("counts 矩阵必须包含 gene names 和 cell barcodes。", call. = FALSE)
    }

    as(counts, "dgCMatrix")
}

select_gene_expression <- function(x) {
    if (is.list(x)) {
        if ("Gene Expression" %in% names(x)) {
            return(x[["Gene Expression"]])
        }
        return(x[[1]])
    }
    x
}

read_background_counts <- function(path, assay = "RNA") {
    if (is.null(path)) {
        return(NULL)
    }
    if (!file.exists(path)) {
        stop("background 路径不存在: ", path, call. = FALSE)
    }

    message("读取 background/raw droplet counts: ", path)
    if (dir.exists(path)) {
        counts <- select_gene_expression(Read10X(path))
    } else {
        obj <- readRDS(path)
        if (inherits(obj, "Seurat")) {
            counts <- get_counts_from_seurat(obj, assay = assay)
        } else if (inherits(obj, "SingleCellExperiment")) {
            counts <- SummarizedExperiment::assay(obj, "counts")
        } else {
            counts <- select_gene_expression(obj)
        }
    }

    if (nrow(counts) == 0 || ncol(counts) == 0) {
        stop("background counts 矩阵为空。", call. = FALSE)
    }

    as(counts, "dgCMatrix")
}

create_decontx_assay <- function(counts) {
    if ("CreateAssay5Object" %in% getNamespaceExports("SeuratObject")) {
        SeuratObject::CreateAssay5Object(counts = counts)
    } else {
        SeuratObject::CreateAssayObject(counts = counts)
    }
}

filter_markers <- function(markers, features) {
    markers <- lapply(markers, function(x) intersect(x, features))
    markers <- markers[lengths(markers) > 0]
    if (length(markers) == 0) {
        stop("表达矩阵中未找到 T/B/Monocyte/NK marker genes，无法生成 marker QC 图。", call. = FALSE)
    }
    markers
}

save_qc_pdf <- function(sce, path) {
    if (!"decontX_UMAP" %in% reducedDimNames(sce)) {
        stop("SCE 中未找到 decontX_UMAP，无法生成 DecontX QC PDF。", call. = FALSE)
    }

    umap <- reducedDim(sce, "decontX_UMAP")
    markers <- list(
        T_cells = c("CD3D", "CD3E", "TRAC"),
        B_cells = c("MS4A1", "CD79A", "CD79B"),
        Monocytes = c("LYZ", "S100A8", "S100A9"),
        NK_cells = c("GNLY", "NKG7")
    )
    markers <- filter_markers(markers, rownames(sce))

    p_cluster <- celda::plotDimReduceCluster(
        x = sce$decontX_clusters,
        dim1 = umap[, 1],
        dim2 = umap[, 2]
    ) +
        ggplot2::ggtitle("DecontX broad clusters")

    p_contamination <- decontX::plotDecontXContamination(sce) +
        ggplot2::ggtitle("DecontX contamination")

    p_marker <- decontX::plotDecontXMarkerPercentage(
        sce,
        markers = markers,
        assayName = c("counts", "decontXcounts")
    ) +
        ggplot2::ggtitle("Marker detection before and after DecontX")

    qc_plot <- (p_cluster | p_contamination) / p_marker +
        plot_layout(heights = c(1, 1.8))

    ggplot2::ggsave(
        filename = path,
        plot = qc_plot,
        width = 12,
        height = 11,
        units = "in"
    )
}

set.seed(opts$seed)

message("读取 Seurat 对象: ", opts$input_rds)
seu <- readRDS(opts$input_rds)
if (!inherits(seu, "Seurat")) {
    stop("输入文件必须是 Seurat RDS 对象。", call. = FALSE)
}
DefaultAssay(seu) <- opts$assay

message("提取 filtered cells 的原始 counts 矩阵 ...")
counts <- get_counts_from_seurat(seu, assay = opts$assay)
sce <- SingleCellExperiment(list(counts = counts))

background_counts <- read_background_counts(opts$background, assay = opts$assay)
background_sce <- NULL
if (!is.null(background_counts)) {
    background_sce <- SingleCellExperiment(list(counts = background_counts))
}

message("运行 DecontX 内置 broad clustering 和污染估计 ...")
if (is.null(background_sce)) {
    sce <- decontX::decontX(sce, seed = opts$seed)
} else {
    sce <- decontX::decontX(sce, background = background_sce, seed = opts$seed)
}

message("写回 Seurat metadata 和 decontXcounts assay ...")
decont_counts <- round(decontX::decontXcounts(sce))
seu[["decontXcounts"]] <- create_decontx_assay(decont_counts)
seu$decontX_contamination <- sce$decontX_contamination
seu$decontX_clusters <- sce$decontX_clusters

cell_metadata <- data.frame(
    cell = colnames(seu),
    decontX_contamination = seu$decontX_contamination,
    decontX_clusters = seu$decontX_clusters,
    stringsAsFactors = FALSE
)

message(sprintf(
    "污染比例统计: min=%.3f  median=%.3f  max=%.3f",
    min(seu$decontX_contamination),
    median(seu$decontX_contamination),
    max(seu$decontX_contamination)
))

message("生成 DecontX QC PDF ...")
save_qc_pdf(sce, output_qc_pdf)

message("保存输出 ...")
saveRDS(seu, output_rds)
write.csv(cell_metadata, output_cell_metadata, row.names = FALSE)
writeLines(capture.output(sessionInfo()), output_session)

message("完成。输出文件:")
message("  Seurat: ", output_rds)
message("  cell metadata: ", output_cell_metadata)
message("  QC PDF: ", output_qc_pdf)
message("  sessionInfo: ", output_session)
