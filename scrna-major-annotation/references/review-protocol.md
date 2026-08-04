# Review protocol

Use this reference only when filling or validating review CSVs. R scripts generate the files and evidence; the Agent fills only `agent_*`; a human fills `human_*` and `approved`.

## Count-version review

`<project>.count-version-review.csv` is generated once by `DecontX.r` and is bound by path and MD5 to the parent QC RDS plus the original/corrected branch seed RDS files.

1. Both branches must independently finish Round 1 clustering, resolution approval, QC report, and provisional annotation report before a selection is valid.
2. Compare UMAP/cluster structure, QC distributions, `decontX_contamination`, marker specificity, and biological interpretability. Do not select corrected counts solely because they reduce marker detection.
3. The Agent may fill `agent_version` with `original` or `corrected`, plus `agent_reason` and auditable `agent_evidence`; these values are advisory only.
4. A human explicitly fills `human_version`, `human_reason`, `reviewer`, `reviewed_at`, and `approved=1`.
5. Round 1 filtering validates both branch histories and permits only the selected branch. The unselected branch stops after Round 1 comparison; Round 2/3 use only the selected branch.

## Shared binding

Every row carries `project`, `round`, `input_rds`, `input_md5`, and `cell_set_md5`. Never copy a review between rounds or objects. Do not alter binding columns. Regenerating a report overwrites its existing review CSV, including an approved one; workflow state retains the audit trail.

## Resolution review

`resolution-review.csv` contains exactly one row and is also bound to `resolution-candidates.csv` by `candidate_md5`.

1. The Agent inspects the 15-panel UMAP, clustree, marker heatmaps, candidate table, and marker table.
2. The Agent fills one `agent_resolution` in 0.1–1.5, plus `agent_reason` and auditable `agent_evidence`.
3. Prefer stable broad-lineage structure with distinct markers and limited fragmentation. A small pDC, Mast, Neutrophil, or other plausible rare cluster is not automatically invalid.
4. A human explicitly fills `human_resolution`, `reviewer`, `reviewed_at`, and `approved=1`. The executable value never falls back to `agent_resolution`.

## QC review

`qc-review.csv` is round-specific.

- Round 1 combines `threshold`, `cluster`, and `doublet` decisions.
- Round 2/3 contain `cluster` decisions only. Threshold rows and per-cell doublet-removal rows are invalid; QC distributions, doublet UMAP, and singlet/doublet proportions are evidence for deciding whether a whole cluster should be retained or removed.

1. `threshold`: `scope` is `global` or one exact `orig.ident`; `profile` is `general` or `neutrophil`; `metric` is only `nFeature_RNA`, `nCount_RNA`, `percent_mito`, or `decontX_contamination`. Agent and human bounds use a finite number or literal `none`. There is no automatic contamination cutoff. A value is removed only when `value < min` or `value > max`; equality is retained. A sample-specific row overrides the matching global row.
   The Agent may append a sample-specific threshold row only when the evidence requires it; copy all binding columns exactly and keep the decision key unique.
2. `cluster`: the Agent may recommend `keep`, `drop`, `neutrophil`, or `review`. The human must choose `keep`, `drop`, or `neutrophil`. A human `drop` requires `human_reason` equal to `low_quality`, `doublet`, or `contaminant`.
3. `doublet`: one row per `orig.ident`. The Agent may recommend `keep`, `remove`, or `review`; the human must choose `keep` or `remove`. `remove` deletes cells whose per-sample `scDblFinder.class` is `doublet`; never impose a manual score cutoff or compare `scDblFinder.score` between samples.
4. Doublet proportions alone cannot establish correctness. Combine them with QC distributions, mixed-lineage marker evidence, cluster structure, and provisional annotation.
5. First run `preview` with Agent fields. After human edits every executable field and sets every row to `approved=1`, rerun `preview`; only then run `apply`. Apply checks the exact review hash.
6. The removal trigger uses the deduplicated union across all allowed approved sources. Round 1 always continues to Round 2. In Round 2/3 the union comes only from deleted clusters. Round 2 global removal at least 5% or any sample removal at least 10% triggers Round 3. Round 3 never triggers Round 4; the same threshold yields `needs_manual_review`.

## Annotation review

1. The Agent label is one of `Epithelial`, `T_NK_cell`, `B_cell`, `Plasma`, `Myeloid`, `pDC`, `Mast`, `Neutrophil`, `Fibroblast`, `Endothelial`, `Mural`, `Proliferation`, or `Ambiguous`.
2. Use these marker-set mappings: `T_cell`/`NK_cell` → `T_NK_cell`; `Pan_myeloid`/`Monocyte`/`macro` → `Myeloid`; `Pericyte`/`Smooth_muscle` → `Mural`.
3. `Proliferation` remains valid only when a whole cluster is coherently proliferative after S/G2M regression and removal of cell-cycle genes from HVGs. Use `MKI67`, `TOP2A`, `UBE2C`, `CENPF`, `TYMS`, and `STMN1` together with lineage evidence.
4. Platelet, erythrocyte, and immunoglobulin panels are contamination evidence, not final labels. Immunoglobulin alone is insufficient outside B/Plasma context.
5. A human explicitly fills one non-empty `human_label` for every retained cluster and sets every row to `approved=1`. `Ambiguous` is forbidden as a final label. A custom human label is accepted but reported as lacking a bundled downstream panel.
