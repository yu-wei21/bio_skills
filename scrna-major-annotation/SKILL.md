---
name: scrna-major-annotation
description: Run a human/GRCh38 Seurat v5 workflow that creates a counts-bearing object, computes QC with per-sample scDblFinder and DecontX ambient-RNA estimates, performs repeatable LogNormalize/Harmony clustering, compares resolutions, filters cells through Agent recommendation plus human approval, and writes final major-cell annotations. Use for single-cell RNA-seq major-lineage QC/annotation from Cell Ranger filtered matrices or a counts-bearing Seurat RDS, including round 1–3 reclustering decisions.
---

# Run major-cell QC and annotation

Use only the nine user-facing scripts in this skill. Treat `orig.ident` as immutable. Accept human/GRCh38 gene symbols only. Never modify an input RDS, call an AI service from R, set human approval fields, or apply an unapproved decision.

Load [review-protocol.md](references/review-protocol.md) when preparing or checking reviews. Use [marker-panels.tsv](references/marker-panels.tsv) as the fixed annotation panel.

## Set paths

```bash
SKILL_DIR="/absolute/path/to/scrna-major-annotation"
OUT="/absolute/path/to/output"
PROJECT="prostate"
cd "$SKILL_DIR"
```

Derived outputs and review CSVs overwrite by default. The root and each count-version branch maintain separate `workflow_state.json` files; do not edit them manually.

## Prepare the QC object

For Cell Ranger output, create a sample sheet whose first two columns are ordered `sample_id,matrix_path`, then run:

```bash
Rscript scripts/create-seurat-object.R --sample-sheet samples.csv --project "$PROJECT" --output-dir "$OUT"
Rscript scripts/seurat-qc-metrics.R --input "$OUT/$PROJECT.raw.rds" --project "$PROJECT" --output-dir "$OUT"
Rscript scripts/DecontX.r --input "$OUT/$PROJECT.qc-metrics.rds" --project "$PROJECT" --output-dir "$OUT"
```

For an existing counts-bearing Seurat v5 RDS with non-empty `orig.ident`, skip creation:

```bash
Rscript scripts/seurat-qc-metrics.R --input /absolute/path/input.rds --project "$PROJECT" --output-dir "$OUT"
Rscript scripts/DecontX.r --input "$OUT/$PROJECT.qc-metrics.rds" --project "$PROJECT" --output-dir "$OUT"
```

The QC preparation stage joins RNA layers, recalculates `nCount_RNA` and `nFeature_RNA` from counts, fails on zero-count cells, uses LogNormalize, computes cell-cycle scores, and runs random-mode `scDblFinder` per `orig.ident`. It keeps only `scDblFinder.score` and `scDblFinder.class` in metadata. `DecontX.r` then takes the actual project-level RDS produced by `seurat-qc-metrics.R` (currently named `<project>.qc-metrics.rds`; `qc-metrics.rds` is not a literal fixed filename) and runs one independent DecontX model per `orig.ident`, parallelizing across samples.

DecontX writes two first-round seeds with identical cells and genes:

```bash
ORIGINAL_OUT="$OUT/count-versions/original"
CORRECTED_OUT="$OUT/count-versions/corrected"
ORIGINAL_INPUT="$ORIGINAL_OUT/$PROJECT.decontx-original.rds"
CORRECTED_INPUT="$CORRECTED_OUT/$PROJECT.decontx-corrected.rds"
COUNT_VERSION_REVIEW="$OUT/$PROJECT.count-version-review.csv"
```

- The original object keeps observed UMI counts in `RNA` and adds `decontX_contamination` and `decontX_clusters`.
- The corrected object stores rounded, non-negative DecontX-corrected counts directly in `RNA`, adds the same contamination metadata, and recalculates count-dependent QC, LogNormalize data, and cell-cycle scores.
- `scDblFinder` is not rerun on corrected counts. Its score/class and provenance are copied from the original-count QC object so doublet evidence remains based on observed counts.

If matched raw droplets are available, pass a two-column `sample_id,background_path` CSV/TSV with `--background-sheet`. Every `orig.ident` must have exactly one same-channel background; by default, every filtered-cell barcode must be found and removed from that raw matrix. Use `--background-is-empty-only` only for a background already restricted to empty droplets. DecontX uses its internal broad clustering for `z`; review each sample's before/after marker plots because the official method notes that cluster choice affects contamination estimates. The plot marker groups include immune lineages plus epithelial (`EPCAM`, `KRT8`, `KRT18`) and fibroblast (`DCN`, `LUM`, `COL1A1`) markers.

## Compare count versions in Round 1

Run the complete first-round clustering and evidence bundle once in each isolated branch. `orig.ident` must identify the independent capture/batch unit, not a biological grouping invented after merging.

```bash
Rscript scripts/seurat-cluster.R --input "$ORIGINAL_INPUT" --project "$PROJECT" --round 1 --output-dir "$ORIGINAL_OUT"
Rscript scripts/seurat-find-resolution.R --input "$ORIGINAL_OUT/round1/$PROJECT.round1.clustered.rds" --project "$PROJECT" --round 1 --output-dir "$ORIGINAL_OUT"

Rscript scripts/seurat-cluster.R --input "$CORRECTED_INPUT" --project "$PROJECT" --round 1 --output-dir "$CORRECTED_OUT"
Rscript scripts/seurat-find-resolution.R --input "$CORRECTED_OUT/round1/$PROJECT.round1.clustered.rds" --project "$PROJECT" --round 1 --output-dir "$CORRECTED_OUT"
```

Inspect and approve one resolution separately in each branch, then generate both QC and provisional annotation reports in both branches:

```bash
Rscript scripts/seurat-qc-report.R --input "$ORIGINAL_OUT/round1/$PROJECT.round1.clustered.rds" --resolution-review "$ORIGINAL_OUT/round1/resolution-review.csv" --project "$PROJECT" --round 1 --output-dir "$ORIGINAL_OUT"
Rscript scripts/seurat-annotation-report.R --mode provisional --input "$ORIGINAL_OUT/round1/$PROJECT.round1.clustered.rds" --resolution-review "$ORIGINAL_OUT/round1/resolution-review.csv" --project "$PROJECT" --round 1 --output-dir "$ORIGINAL_OUT"

Rscript scripts/seurat-qc-report.R --input "$CORRECTED_OUT/round1/$PROJECT.round1.clustered.rds" --resolution-review "$CORRECTED_OUT/round1/resolution-review.csv" --project "$PROJECT" --round 1 --output-dir "$CORRECTED_OUT"
Rscript scripts/seurat-annotation-report.R --mode provisional --input "$CORRECTED_OUT/round1/$PROJECT.round1.clustered.rds" --resolution-review "$CORRECTED_OUT/round1/resolution-review.csv" --project "$PROJECT" --round 1 --output-dir "$CORRECTED_OUT"
```

Compare the two Round 1 UMAPs, cluster stability, QC distributions, contamination patterns, marker specificity, and provisional annotation interpretability. Fill only `agent_*` in the root `count-version-review.csv`; pause until a human fills `human_version` (`original` or `corrected`), `human_reason`, reviewer/time, and `approved=1`. Set `BRANCH_OUT` to exactly one selected directory (`$ORIGINAL_OUT` or `$CORRECTED_OUT`). From this point onward, do not continue the unselected branch.

Seurat v5 layer handling is enforced in every selected branch: clustering first joins any existing RNA layers, then splits counts/data by `orig.ident` before per-layer normalization and Harmony integration, and finally calls `JoinLayers` before marker analysis or saving downstream objects.

## Run each round

For the selected branch, start from the already completed Round 1 clustered object; do not rerun its resolution/QC/annotation reports because doing so overwrites their review templates.

```bash
ROUND=1
CLUSTERED="$BRANCH_OUT/round1/$PROJECT.round1.clustered.rds"
```

Both reports are mandatory in every round and are bound to the same clustered RDS and approved resolution. The QC UMAP and cluster/sample violin bundles include `decontX_contamination`. In Round 1 it is available as a human-approved cell-level threshold metric, but no default or automatic contamination cutoff is applied. Use the selected branch's Round 1 bundles to fill only the `agent_*` fields in `qc-review.csv`:

1. Round 1 may recommend cell-level thresholds, per-sample `scDblFinder.class` removal, and cluster actions.
2. Round 2/3 may recommend cluster actions only. Do not add threshold or doublet-removal rows; QC values and doublet plots remain evidence for cluster-level decisions.

Run the selected branch's Round 1 impact preview:

```bash
Rscript scripts/seurat-filter-cells.R --mode preview --input "$CLUSTERED" --resolution-review "$BRANCH_OUT/round$ROUND/resolution-review.csv" --qc-review "$BRANCH_OUT/round$ROUND/qc-review.csv" --count-version-review "$COUNT_VERSION_REVIEW" --project "$PROJECT" --round 1 --output-dir "$BRANCH_OUT"
```

Pause for human review. In Round 1, require `human_action` on every cluster/doublet row and `human_min`/`human_max` on every threshold row. In Round 2/3, require only cluster `human_action`. Every dropped cluster requires `human_reason`; every row requires `approved=1`. Rerun the same preview command, inspect the exact union impact, then apply:

```bash
Rscript scripts/seurat-filter-cells.R --mode apply --input "$CLUSTERED" --resolution-review "$BRANCH_OUT/round$ROUND/resolution-review.csv" --qc-review "$BRANCH_OUT/round$ROUND/qc-review.csv" --count-version-review "$COUNT_VERSION_REVIEW" --project "$PROJECT" --round 1 --output-dir "$BRANCH_OUT"
FILTERED="$BRANCH_OUT/round$ROUND/$PROJECT.round$ROUND.filtered.rds"
```

For each required Round 2/3, use the preceding filtered output as `ROUND_INPUT`, then run clustering, resolution comparison, and both reports only in the selected branch:

```bash
Rscript scripts/seurat-cluster.R --input "$ROUND_INPUT" --project "$PROJECT" --round "$ROUND" --output-dir "$BRANCH_OUT"
CLUSTERED="$BRANCH_OUT/round$ROUND/$PROJECT.round$ROUND.clustered.rds"
Rscript scripts/seurat-find-resolution.R --input "$CLUSTERED" --project "$PROJECT" --round "$ROUND" --output-dir "$BRANCH_OUT"
```

Inspect the resolution bundle. Fill only `agent_resolution`, `agent_reason`, and `agent_evidence`; pause until a human fills `human_resolution` and `approved=1`. Then generate both reports:

```bash
Rscript scripts/seurat-qc-report.R --input "$CLUSTERED" --resolution-review "$BRANCH_OUT/round$ROUND/resolution-review.csv" --project "$PROJECT" --round "$ROUND" --output-dir "$BRANCH_OUT"
Rscript scripts/seurat-annotation-report.R --mode provisional --input "$CLUSTERED" --resolution-review "$BRANCH_OUT/round$ROUND/resolution-review.csv" --project "$PROJECT" --round "$ROUND" --output-dir "$BRANCH_OUT"
```

For Round 2/3 filtering, use the same preview/apply pattern with `--round "$ROUND"` and branch-local paths, but omit `--count-version-review`.

Read `filter-summary.csv` and `workflow_state.json`:

1. Round 1 always sets `ROUND_INPUT="$FILTERED"`, increments to round 2, and reruns clustering.
2. Round 2 cluster-level removal with global removal ≥5% or any sample removal ≥10% sets `ROUND_INPUT="$FILTERED"`, increments to round 3, and reruns clustering.
3. Round 2 or 3 below the trigger proceeds to final annotation.
4. Round 3 at or above the trigger stops at `needs_manual_review`; never create round 4 or claim completion.

## Apply final annotation

Regenerate marker evidence on the final filtered cell set:

```bash
Rscript scripts/seurat-annotation-report.R --mode final --input "$FILTERED" --resolution-review "$BRANCH_OUT/round$ROUND/resolution-review.csv" --project "$PROJECT" --round "$ROUND" --output-dir "$BRANCH_OUT"
```

Fill only `agent_label`, `agent_reason`, and `agent_evidence`. Pause until a human explicitly labels every cluster and sets every row to `approved=1`. Then run:

```bash
Rscript scripts/seurat-annotate-cells.R --input "$FILTERED" --resolution-review "$BRANCH_OUT/round$ROUND/resolution-review.csv" --annotation-review "$BRANCH_OUT/round$ROUND/annotation-review.csv" --project "$PROJECT" --round "$ROUND" --output-dir "$BRANCH_OUT"
```

Report the final RDS, final celltype UMAP, annotation summary, final round, removal rates, any custom labels, and the final `workflow_state.json` status. If any command fails, stop and report the exact error; do not bypass validation or silently drop cells.
