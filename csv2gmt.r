## 输入：CSV文件，每列一个基因集，唯一列名为基因集名称，列内为同一种基因ID；缺失成员留空。
## 输出：output/csv2gmt/下以输入CSV文件名开头的GMT文件。
## 示例：Rscript csv2gmt.r data/gene_sets.csv


##  read file
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1) {
    stop("Usage: Rscript csv2gmt.r <gene_sets.csv>")
}
file <- args[1]
cat(file)
geneset <- read.csv(file, na.string = "")
output_dir <- file.path("output", "csv2gmt")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
input_tag <- tools::file_path_sans_ext(basename(file))
output_file <- file.path(output_dir, paste0(input_tag, "_gene_sets.gmt"))

## 转换成list
name <- colnames(geneset)
geneset <- lapply(colnames(geneset), function(x) {
    data <- geneset[, x]
    data <- data[!is.na(data)]
})
names(geneset) <- name

## 转换成gmt
## list to gmt
write.list2gmt <- function(gs, file) {
    sink(file)
    lapply(names(gs), function(i) {
        cat(paste(c(i, "tmp", gs[[i]]), collapse = "\t"))
        cat("\n")
    })
    sink()
}
write.list2gmt(geneset, file = output_file)

## end
message("successful !")
