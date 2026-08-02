## 输入：基因集XLSX（每个sheet仅读第一列；第一列列名是基因集名称，下面各行为与gene_format一致的SYMBOL/ENTREZID）、SYMBOL格式GMT目录（文件名须符合gmt_registry），以及背景CSV（第一列SYMBOL）或含gene_name/gene_biotype的注释CSV。
## 输出：output/ORA/下以基因集XLSX文件名开头的ORA结果XLSX及图形PDF/PNG。
## 示例：Rscript ORA.r data/gene_sets.xlsx SYMBOL data/universe.csv --gmt-dir=data/gmt

## 脚本目的：对多个基因集合进行ORA过表达富集分析（compareCluster）
## 数据库：GO (BP/CC/MF 合并), KEGG, Reactome, Hallmark, MSigDB c1-c9
## 全部使用本地 GMT + clusterProfiler::enricher，零联网依赖
## 未提供背景CSV时，使用GTF中的protein-coding genes。

suppressPackageStartupMessages({
  library(clusterProfiler)
  library(openxlsx)
  library(enrichplot)
  library(tidyverse)
  library(ggplot2)
  library(cowplot)
  library(data.table)
})

## ---- 参数解析 ----
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop(
    "Usage: Rscript ORA.r <input.xlsx> <gene_format> [universe.csv] ",
    "--gmt-dir=DIR [--gtf=FILE]\n",
    "  gene_format: SYMBOL 或 ENTREZID\n",
    "  universe.csv: 可选，第一列为本次实际检验的 SYMBOL"
  )
}

input_file <- args[1]
gene_format <- toupper(trimws(args[2]))
universe_file <- NULL
gmt_dir <- NULL
gtf_file <- NULL

extra_args <- if (length(args) > 2) args[-c(1, 2)] else character(0)
for (arg in extra_args) {
  if (startsWith(arg, "--gmt-dir=")) {
    gmt_dir <- sub("^--gmt-dir=", "", arg)
  } else if (startsWith(arg, "--gtf=")) {
    gtf_file <- sub("^--gtf=", "", arg)
  } else if (startsWith(arg, "--")) {
    stop("无法识别参数: ", arg)
  } else if (is.null(universe_file)) {
    universe_file <- arg
  } else {
    stop("无法识别参数: ", arg)
  }
}

if (is.null(gmt_dir) || !nzchar(gmt_dir)) {
  stop("必须通过 --gmt-dir=DIR 指定GMT目录。")
}
if (!file.exists(input_file)) {
  stop("基因集XLSX不存在: ", input_file)
}
if (!dir.exists(gmt_dir)) {
  stop("GMT目录不存在: ", gmt_dir)
}
if (is.null(universe_file) && (is.null(gtf_file) || !nzchar(gtf_file))) {
  stop("未提供universe.csv时，必须通过 --gtf=FILE 指定回退GTF注释CSV。")
}
if (!is.null(universe_file) && !file.exists(universe_file)) {
  stop("背景基因CSV不存在: ", universe_file)
}
if (is.null(universe_file) && !file.exists(gtf_file)) {
  stop("回退GTF注释CSV不存在: ", gtf_file)
}

if (!gene_format %in% c("SYMBOL", "ENTREZID")) {
  stop("gene_format 必须为 SYMBOL 或 ENTREZID，当前值: ", gene_format)
}

output_dir <- file.path("output", "ORA")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
input_tag <- tools::file_path_sans_ext(basename(input_file))
output_prefix <- file.path(output_dir, input_tag)
message("Input  : ", input_file)
message("Format : ", gene_format)
message("Universe: ", if (is.null(universe_file)) "GTF protein-coding genes" else universe_file)
message("Output : ", output_prefix, "_ora_results.xlsx / _ora_plots.pdf / _ora_plots_*.png")

## ---- 读取背景基因（universe）----
if (!is.null(universe_file)) {
  message("Loading tested-gene universe ...")
  universe_df <- fread(universe_file, data.table = FALSE)
  universe_genes <- unique(trimws(as.character(universe_df[[1]])))
  universe_genes <- universe_genes[!is.na(universe_genes) & universe_genes != ""]
} else {
  message("Loading fallback universe (all GTF protein-coding genes) ...")
  gtf <- fread(gtf_file, data.table = FALSE)
  required_gtf_cols <- c("gene_name", "gene_biotype")
  missing_gtf_cols <- setdiff(required_gtf_cols, colnames(gtf))
  if (length(missing_gtf_cols) > 0) {
    stop("回退GTF注释CSV缺少必需列: ", paste(missing_gtf_cols, collapse = ", "))
  }
  universe_genes <- unique(trimws(as.character(
    gtf$gene_name[gtf$gene_biotype == "protein_coding"]
  )))
  universe_genes <- universe_genes[!is.na(universe_genes) & nzchar(universe_genes)]
}
if (length(universe_genes) == 0) {
  stop("背景基因集合为空。")
}
message("  Universe size: ", length(universe_genes))

## ---- 读取 Excel ----
## 每个sheet：第一行（列标题）= 基因集名称，余下行 = 基因
sheet_names <- getSheetNames(input_file)
message("Sheets : ", paste(sheet_names, collapse = ", "))

gene_cluster <- lapply(sheet_names, function(sheet) {
  df <- read.xlsx(
    input_file,
    sheet = sheet,
    colNames = TRUE,
    check.names = FALSE,
    sep.names = " "
  )
  group_name <- colnames(df)[1] # 第一列列标题 = 集合名称
  genes <- unique(trimws(as.character(df[[1]]))) # 第一列数据 = 基因
  genes <- genes[!is.na(genes) & genes != ""]
  list(name = group_name, genes = genes)
})

gene_cluster <- setNames(
  lapply(gene_cluster, `[[`, "genes"),
  sapply(gene_cluster, `[[`, "name")
)
if (anyNA(names(gene_cluster)) || any(!nzchar(names(gene_cluster)))) {
  stop("每个sheet第一列都必须有非空列名，作为基因集名称。")
}
if (anyDuplicated(names(gene_cluster))) {
  stop(
    "不同sheet的第一列列名（基因集名称）不得重复: ",
    paste(unique(names(gene_cluster)[duplicated(names(gene_cluster))]), collapse = ", ")
  )
}

## ---- 基因ID转换 ----
## 本地 GMT 全部使用 SYMBOL，若输入为 ENTREZID 则先转换
if (gene_format == "ENTREZID") {
  message("Converting ENTREZID to SYMBOL ...")
  library(org.Hs.eg.db)
  gene_cluster <- lapply(gene_cluster, function(genes) {
    mapped <- bitr(genes,
      fromType = "ENTREZID",
      toType   = "SYMBOL",
      OrgDb    = org.Hs.eg.db,
      drop     = TRUE
    )
    mapped$SYMBOL
  })
}

genes_before_universe <- sum(lengths(gene_cluster))
gene_cluster <- lapply(gene_cluster, function(genes) {
  base::intersect(genes, universe_genes)
})
genes_after_universe <- sum(lengths(gene_cluster))
message(
  "Candidate genes overlapping universe: ",
  genes_after_universe, " / ", genes_before_universe
)
if (genes_after_universe == 0) {
  stop("候选基因与背景基因没有交集，请核查基因ID体系。")
}

## 过滤空基因集
n_before <- length(gene_cluster)
gene_cluster <- Filter(function(x) length(x) >= 3, gene_cluster)
if (length(gene_cluster) < n_before) {
  message("警告：部分基因集基因数 < 3，已跳过")
}
if (length(gene_cluster) == 0) {
  stop("所有基因集均少于3个有效基因，无法进行ORA。")
}
message(
  "有效基因集 (", length(gene_cluster), "): ",
  paste(names(gene_cluster), collapse = ", ")
)

## ---- 加载本地 GMT 文件 ----
gmt_registry <- list(
  GO = c(
    "GO_Biological_Process_2026.gmt",
    "GO_Cellular_Component_2026.gmt",
    "GO_Molecular_Function_2026.gmt"
  ),
  KEGG = "KEGG_2026.gmt",
  Reactome = "Reactome_Pathways_2024.gmt",
  Hallmark = "hallmark.gmt",
  c1 = "c1.gmt",
  c2 = "c2.gmt",
  c3 = "c3.gmt",
  c4 = "c4.gmt",
  c5 = "c5.gmt",
  c6 = "c6.gmt",
  c7 = "c7.gmt",
  c8 = "c8.gmt",
  c9 = "c9.gmt"
)

load_local_gmt <- function(filenames, label) {
  t2g_list <- lapply(filenames, function(f) {
    fpath <- file.path(gmt_dir, f)
    if (!file.exists(fpath)) {
      message("  ", label, ": GMT 文件不存在: ", fpath)
      return(NULL)
    }
    read.gmt(fpath)
  })
  t2g_list <- Filter(Negate(is.null), t2g_list)
  if (length(t2g_list) == 0) {
    return(NULL)
  }
  do.call(rbind, t2g_list)
}

gmt_loaded <- lapply(names(gmt_registry), function(db) {
  filenames <- unlist(gmt_registry[[db]])
  t2g <- load_local_gmt(filenames, db)
  if (is.null(t2g)) {
    message("  跳过: ", db)
    return(NULL)
  }
  message(
    "  ", db, ": ", length(unique(t2g$term)), " terms, ",
    length(unique(t2g$gene)), " genes"
  )
  t2g
})
names(gmt_loaded) <- names(gmt_registry)
if (all(vapply(gmt_loaded, is.null, logical(1)))) {
  stop("未加载到任何gmt_registry中定义的GMT文件。")
}

## ---- compareCluster 分析（全部 enricher + 本地 GMT + universe）----
run_compare <- function(label, t2g) {
  if (is.null(t2g)) {
    return(NULL)
  }
  message("Running ", label, " ORA ...")
  tryCatch(
    compareCluster(
      geneCluster   = gene_cluster,
      fun           = "enricher",
      TERM2GENE     = t2g,
      universe      = universe_genes,
      pAdjustMethod = "BH",
      pvalueCutoff  = 0.05,
      qvalueCutoff  = 0.2
    ),
    error = function(e) {
      message("  ", label, " 失败: ", e$message)
      NULL
    }
  )
}

results <- lapply(names(gmt_loaded), function(db) {
  run_compare(db, gmt_loaded[[db]])
})
names(results) <- names(gmt_loaded)

## ---- 后处理 ----
## Hallmark: 去掉 HALLMARK_ 前缀，下划线换空格，转 Title Case
fmt_hallmark <- function(x) {
  sub("^HALLMARK_", "", x) %>%
    gsub("_", " ", .) %>%
    str_to_title()
}
if (!is.null(results[["Hallmark"]]) && nrow(as.data.frame(results[["Hallmark"]])) > 0) {
  results[["Hallmark"]]@compareClusterResult$Description <-
    fmt_hallmark(results[["Hallmark"]]@compareClusterResult$Description)
}

## Reactome: 去掉 REACTOME_ 前缀（如有），下划线换空格，转 Title Case
if (!is.null(results[["Reactome"]]) && nrow(as.data.frame(results[["Reactome"]])) > 0) {
  results[["Reactome"]]@compareClusterResult$Description <-
    results[["Reactome"]]@compareClusterResult$Description %>%
    sub("^REACTOME_", "", .) %>%
    gsub("_", " ", .) %>%
    str_to_title()
}

## ---- 保存 Excel ----
wb <- createWorkbook()
n_result_sheets <- 0L
for (db in names(results)) {
  res <- results[[db]]
  if (!is.null(res)) {
    df <- as.data.frame(res)
    if (nrow(df) > 0) {
      addWorksheet(wb, db)
      writeData(wb, db, df)
      n_result_sheets <- n_result_sheets + 1L
    }
  }
}
if (n_result_sheets == 0L) {
  addWorksheet(wb, "no_significant_results")
  writeData(
    wb,
    "no_significant_results",
    data.frame(message = "No significant ORA results", stringsAsFactors = FALSE)
  )
}

excel_out <- paste0(output_prefix, "_ora_results.xlsx")
saveWorkbook(wb, excel_out, overwrite = TRUE)
message("Excel saved: ", excel_out)

## ---- 可视化 ----
n_cluster <- length(gene_cluster)
base_fig_w <- max(10, 7 + n_cluster * 1.2)
fig_w <- base_fig_w * 0.6
fig_h <- 10 * 2 / 3
n_show <- 8

make_dotplot <- function(result, title, n_show = 8) {
  if (is.null(result)) {
    return(NULL)
  }
  df <- as.data.frame(result)
  if (nrow(df) == 0) {
    return(NULL)
  }

  tryCatch(
    {
      p <- dotplot(result, showCategory = n_show) +
        ggtitle(title) +
        scale_y_discrete(labels = function(x) str_wrap(x, width = 55)) +
        theme_cowplot(font_size = 11) +
        theme(
          axis.text.x     = element_text(angle = 45, hjust = 1, size = 9),
          axis.text.y     = element_text(size = 9),
          legend.position = "right",
          plot.title      = element_text(size = 14, face = "bold"),
          strip.text      = element_text(size = 10, face = "bold")
        )
      p
    },
    error = function(e) {
      message("  dotplot 失败 (", title, "): ", e$message)
      NULL
    }
  )
}

plots <- lapply(names(results), function(db) {
  make_dotplot(results[[db]], paste(db, "ORA"), n_show = n_show)
})
names(plots) <- names(results)
plots <- Filter(Negate(is.null), plots)
if (length(plots) == 0) {
  plots <- list(
    no_significant_results = ggplot() +
      annotate("text", x = 0, y = 0, label = "No significant ORA results") +
      xlim(-1, 1) +
      ylim(-1, 1) +
      theme_void()
  )
}

pdf_out <- paste0(output_prefix, "_ora_plots.pdf")
png_out <- paste0(output_prefix, "_ora_plots_%03d.png")
pdf(pdf_out, width = fig_w, height = fig_h)
for (db in names(plots)) {
  print(plots[[db]])
}
dev.off()
png(png_out, width = fig_w, height = fig_h, units = "in", res = 300)
for (db in names(plots)) {
  print(plots[[db]])
}
dev.off()
message("PDF saved : ", pdf_out)
message("PNG saved : ", sub("%03d", "*", png_out, fixed = TRUE))
message("Done.")
