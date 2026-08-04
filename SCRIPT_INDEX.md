# 生物信息分析脚本索引

本项目是个人可复用的生物信息分析脚本集合。开始新任务时，先读 [AGENTS.md](AGENTS.md)，再用本页按场景选择脚本；每个脚本文件头是输入格式与参数的最终准则。

## 使用约定

- 命令通常为 `Rscript <脚本名> ...`；先使用脚本头部的示例命令。
- 常规结果位于 `output/<脚本名>/`；不要把原始数据写入 `output/`。
- 基因 ID、物种和表达层级必须与下游 GMT/注释文件一致。差异表达中的 count 方法使用原始整数 counts。
- 附带资源：[GRCh38.v98.gtf.csv](GRCh38.v98.gtf.csv)（人类注释）和 [gmt/](gmt/)（本地基因集）。

## 单细胞预处理与互操作

| 脚本 | 功能与最小输入 | 主输出 | 图形产物 |
|---|---|---|---|
| [DecontX.r](scrna-major-annotation/scripts/DecontX.r) | `seurat-qc-metrics.R` 实际生成的项目级 RDS（当前为 `<project>.qc-metrics.rds`）；按 `orig.ident` 单样本建模并行，可选逐样本 raw 背景表 | original/corrected 两个独立 Seurat RDS（各自 RNA 为实际分析 counts，均含污染比例）、双分支 manifest/review、逐细胞/样本摘要、sessionInfo | 每样本 marker 校正前/后 PDF/PNG（含上皮与成纤维 markers）；污染比例可进入 Round 1 人工阈值审批但无自动 cutoff；两个分支仅并行完成 Round 1，人工选择后单分支继续 |
| [scDblFinder.r](scDblFinder.r) | 含原始 counts 的 Seurat RDS、捕获批次列；可选初步聚类列 | 含 score/class 的 Seurat RDS、逐细胞/汇总 CSV、sessionInfo | score 分布 PDF/PNG |
| [seurat-to-h5ad.r](seurat-to-h5ad.r) | Seurat v4/v5 RDS；可选 assay | AnnData `.h5ad` | 无 |

## 注释、标记与谱系可视化

| 脚本 | 功能与最小输入 | 主输出 | 图形产物 |
|---|---|---|---|
| [seurat.4.annotation.r](seurat.4.annotation.r) | 已归一化 Seurat RDS，含 `RNA_snn_res.0.5` | marker 结果 | marker panel 与 top-marker heatmap PDF/PNG |
| [plot_marker_dotplot.r](plot_marker_dotplot.r) | Seurat RDS、含 `panel,marker` 的 CSV/TSV | 作图数据 CSV | marker DotPlot PDF/PNG |
| [plot-lineage-umap.r](plot-lineage-umap.r) | 已注释 Seurat RDS、celltype→display_name 映射 CSV（行序=集群编号） | 编号图例 | lineage UMAP 两页 PDF + PNG |

## 细胞组成与比例展示

| 脚本 | 功能与最小输入 | 主输出 | 图形产物 |
|---|---|---|---|
| [plot-celltype-barplot.r](plot-celltype-barplot.r) | Seurat RDS、样本列、细胞类型列 | 每样本比例 CSV | 堆叠柱图 PNG/PDF |
| [plot-celltype-boxplot.r](plot-celltype-boxplot.r) | 样本×细胞类型丰度表、样本列、组别列 | 长表 CSV | 分组箱线图 PNG/PDF |
| [propeller.r](propeller.r) | 已注释 Seurat RDS、样本/分组/细胞类型列（每组至少 2 个生物学样本） | 细胞类型比例差异检验、样本比例、组样本数 CSV、sessionInfo | 无 |
| [plot-Roe.r](plot-Roe.r) | 细胞类型×组别的非负整数细胞计数 CSV | 默认写入输入文件目录的 Ro/e XLSX | Ro/e 热图 PNG/PDF |

## 样本级表达与差异分析

| 脚本 | 功能与最小输入 | 主输出 | 图形产物 |
|---|---|---|---|
| [sc2pseudobulk.r](sc2pseudobulk.r) | Seurat RNA counts、`orig.ident` 与细胞类型列 | all-cell/逐细胞类 gene×sample counts CSV | 无 |
| [DESeq2.r](DESeq2.r) | 整数 count 矩阵、两组分组表、GTF、GMT 目录 | DEG、ORA 背景、GSEA 排序表 | PCA/VST PDF/PNG |
| [presto.r](presto.r) | Seurat RDS、细胞类型列、GTF | 单细胞 DEG/ORA XLSX、背景和 GSEA 排序表 | 无 |

推荐衔接：`sc2pseudobulk.r` → `DESeq2.r` → `ORA.r` 或 `GSEA.r`。

## 通路富集与活性评分

| 脚本 | 功能与最小输入 | 主输出 | 图形产物 |
|---|---|---|---|
| [ORA.r](ORA.r) | 多基因集 XLSX、背景/GTF、本地 GMT | 多基因集 ORA XLSX | 富集 PDF/PNG |
| [GSEA.r](GSEA.r) | gene symbol 排序表、本地 GMT 目录 | GSEA XLSX | 富集 PDF/PNG |
| [plot-GSEA-dotplot.r](plot-GSEA-dotplot.r) | GSEA XLSX；内置 FDG 主题通路 | 内置 FDG 通路的 NES 配对点图 PDF/PNG | NES 配对点图 PDF/PNG |
| [GSVA.r](GSVA.r) | gene×sample 表达矩阵、表达类型、GMT | GSVA 评分 CSV | 无 |
| [GSVA_2_limma.r](GSVA_2_limma.r) | GSVA 评分与含 `IDbn`/`FDG_groupe` 的临床 CSV；项目特异分组映射 | limma 差异通路 CSV | 无 |

`plot-GSEA-dotplot.r` 是项目主题图；使用前确认其内置通路列表符合当前问题。

## 格式转换

| 脚本 | 功能与最小输入 | 主输出 | 图形产物 |
|---|---|---|---|
| [csv2gmt.r](csv2gmt.r) | 每列一个基因集的 CSV | GMT 文件 | 无 |

## 验证状态与维护

- 2026-08-04：`scrna-major-annotation/scripts/DecontX.r` 已按 Bioconductor decontX 官方教程/API 复核，并接入该 skill 的 QC 状态链。现输出 original/corrected 双对象，矫正分支重算 count 相关 QC，两个分支独立完成 Round 1 后由人工选择；Seurat v5 聚类显式检查 split layers，并在 marker/保存前检查 JoinLayers。已完成 R 语法与 workflow 静态验证；当前环境未安装 `decontX`，尚未完成真实数据端到端运行。
- 2026-08-03 曾对当时除 soupX.r 外的全部脚本基于子集数据实测（子集对象 5188 细胞，FDG 分组）；测试报告见 `output/test_run/TEST_REPORT.md`（测试产物由用户手动清理）。
- 实测中发现并修复的问题：GSVA.r 旧 API（gsva() 旧签名在 GSVA≥2.0 defunct）与 check.names 样本 ID 改写、GSVA_2 越界调试打印与 check.names、DecontX.r 的 Assays 命名空间遮蔽。详情见测试报告。
- 2026-08-04 每日扫描：新增 `scDblFinder.r`（双细胞判定）与 `propeller.r`（样本级细胞组成差异检验）；两者依据官方 vignette 的当前 API 编写，并完成静态语法/帮助文本检查，尚未以真实数据实测。随后按用户要求删除 `soupX.r`；当前共 19 个可执行 R 脚本。
- 已完成 `scDblFinder` 与 `speckle` 候选流程的官方教程 API 复核；其他既有脚本尚未逐包复核。
- 每日 05:00 自动扫描会更新本页，最多提议 2 个新增通用流程；新建候选脚本由用户决定是否保留。
- 运行环境请参阅 [分析环境配置.md](分析环境配置.md)；注释参考见 [细胞分群注释.md](细胞分群注释.md)；跨脚本背景见 [CONTEXT.md](CONTEXT.md)。
