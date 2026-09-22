# Copy manifest — files still to be added from the author's local environment

Confirmed against the author's local directory listings for both
`F:\MASLD\Rewrite 20260415 v1 轉投\` (main project folder) and
`F:\MASLD\Rewrite 20260415 v1 轉投\Independent Validation\`
(2026-07-01). Copy each file below into the indicated repository path, then
delete this manifest.

## → scripts/05_cross_cohort_validation/  (top level)
Source: `Independent Validation\` (top level of that folder)
- [x] fig_ValidationCrossCohort_6datasets.py            (CORRECTED 2026-07-01: was
      previously marked "already copied" in error — it had NOT actually been
      added. Now genuinely included, verified against uploaded file content.)
- [x] quantify_separation_vs_early_intervention.py      (same correction as above;
      genuinely included — this is the script that produces Fig. 9,
      `Fig_SeparationIndex_vs_AUCratio_scatter.png`, the manuscript's main-text
      scatter plot. If you have that PNG, add it too — see note below.)
- [x] quantify_all_outputs_separation_vs_intervention.py (same correction as above;
      genuinely included — produces Supplementary Fig. S5.)
- [x] check_paths_6datasets.py                          (utility script; genuinely
      included.)

**Note on Fig. 9 (main text) — RESOLVED (2026-07-02):** `quantify_separation_vs_early_intervention.py`
saves its figure as `Fig_SeparationIndex_vs_AUCratio_scatter.png`. This PNG
and its underlying `Table_SeparationIndex_vs_AUCratio.xlsx` have been run
and provided; both are already included in `figures/` and `data/` in this
repository, and the manuscript's Fig. 9 has been updated from a placeholder
to the real figure (verified: N = 18, Mann-Whitney P = 0.0059, median ratio
1.80 vs 0.96, Spearman r = 0.47, P = 0.048 — all matching the manuscript text
exactly).

**Note on the superseded 3-dataset prototype:** an earlier script,
`fig_ValidationCrossCohort_DiseaseProgression_and_Silymarin_3datasets.py`
(covering only GSE126848/GSE48452/GSE89632), was reviewed and intentionally
**excluded** from this repository — it is fully superseded by
`fig_ValidationCrossCohort_6datasets.py` above, which covers all six cohorts.
Keeping only the current version avoids confusion about which script produced
the manuscript's actual figures.

## → scripts/05_cross_cohort_validation/shared/
Source: `Independent Validation\GSE130970_GSE162694_GSE213621_validation_pipeline\`
- [ ] node_ratio_core.py

## → scripts/05_cross_cohort_validation/GSE130970_GSE162694_GSE213621/
Source: same folder as above
- [ ] compute_node_ratios_GSE130970.py
- [ ] compute_node_ratios_GSE162694.py
- [ ] compute_node_ratios_GSE213621.py

## → scripts/05_cross_cohort_validation/GSE48452/
Source: `Independent Validation\GSE48452_validation_pipeline\`
- [ ] step0_setup_environment.py
- [ ] step1_parse_series_matrix.py
- [ ] step2_annotate_probes.py
- [ ] step3_DEG.py
- [ ] step4_concordance.py
- [ ] step5_model_output_validation.py
- [ ] step6_silymarin_validation.py
- [ ] compute_node_ratios_from_expression.py
- [ ] compare_node_ratios_vs_Ahrens2013_TableS2_v2.py   (use v2, the later revision)

## → scripts/05_cross_cohort_validation/GSE89632/
Source: `Independent Validation\GSE89632_validation_pipeline\`
- [ ] step0_setup_environment.py
- [ ] step1_parse_series_matrix.py
- [ ] step2_annotate_probes.py
- [ ] step3_DEG.py
- [ ] step4_concordance.py
- [ ] step5_model_output_validation.py
- [ ] step6_silymarin_validation.py
- [ ] compare_GSE89632_vs_Arendt2015.py

## → scripts/05_cross_cohort_validation/diagnostics/  (optional; kept for transparency, not required to reproduce the manuscript's results)
- [ ] diagnose_sample_id_matching.py           (from GSE130970_GSE162694_GSE213621_validation_pipeline\)
- [ ] diagnose_gse162694_gene_id.py            (from GSE130970_GSE162694_GSE213621_validation_pipeline\)
- [ ] diagnose_gse48452_89632_sample_n.py      (identical copy exists in all three pipeline folders; take any one)

## → scripts/01_differential_expression/  (journal-format companions)
Source: `Rewrite 20260415 v1 轉投\` (main project folder, top level)
- [ ] GSE126848_Fig2_journal.py           (journal-spec companion to the volcano-plot script; produces Fig2_volcano_journal.png actually used in the manuscript)
- [ ] GSE126848_Fig3_journal.py           (journal-spec companion to the enrichment dot-plot script; produces Fig3_dotplot_journal.png actually used in the manuscript)

## → scripts/02_ode_model_construction/  (journal-format companions)
Source: same folder as above
- [ ] render_Fig4_journal.bat
- [ ] Fig4_ODE_flat_v2.dot

## → scripts/03_simulation/  (journal-format companion)
Source: same folder as above
- [ ] masld_all_12P_overview_journal.py   (journal-spec companion; produces All_12_P_outputs_journal.png actually used in the manuscript as Fig. 5)

## Topological-distance script — RESOLVED (reconstructed & independently verified)
The original `quantify_topological_distance_effect.py` could not be located
(searched both the `Independent Validation\` tree and the main project
folder — most likely run interactively and not saved). A reconstruction,
`quantify_topological_distance_effect_RECONSTRUCTED.py`, was written from
the documented methodology (BFS shortest path from the 8 silymarin target
nodes to each P_* output, using the KEGG hsa04932-derived interaction edge
list) and is included directly in this folder (already copied, no action
needed).

**Independently verified against the original results (2026-07-01):** the
author ran the reconstructed script locally and it reproduced the exact
values on record: N = 60, Pearson r = -0.149 (p = 0.2563), Spearman
r = -0.034 (p = 0.7992) — matching Table_TopologicalDistance_vs_AUCratio.xlsx
to the reported decimal places. The per-output BFS distances also matched
exactly for all 12 pathway outputs, including the two unreachable nodes
(P_Adipogenesis, P_Improvement_of_NAFLD). This script is therefore a
validated, functionally equivalent replacement for the original.
- [x] `quantify_topological_distance_effect_RECONSTRUCTED.py` — included, ready to run
- [x] `Fig_TopologicalDistance_vs_AUCratio.png` — already in `figures/`
- [x] `Table_TopologicalDistance_vs_AUCratio.xlsx` — already in `data/`

## Self-containment update (2026-07-02): ODE model files + ratio tables now REQUIRED
**Superseded decision:** earlier revisions of this manifest said the
`ode_model_pydeseq2_*_vs_Normal_v10_mean.py` files were redundant and did
not need to be added, since a similar-purpose script already existed under
`scripts/02_ode_model_construction/`. That guidance is now **reversed**: the
four cross-cohort scripts (`fig_ValidationCrossCohort_6datasets.py`,
`quantify_separation_vs_early_intervention.py`,
`quantify_all_outputs_separation_vs_intervention.py`,
`check_paths_6datasets.py`) have been rewritten to import these exact model
files at runtime — they are load-bearing, not optional reference material.
To make this repository fully self-contained (runnable without any external
folder outside the repo), the following are now required:

- [x] `data/models/ode_model_pydeseq2_NASH_vs_Normal_v10_mean.py`
- [x] `data/models/ode_model_pydeseq2_NAFL_vs_Normal_v10_mean.py`
- [ ] `data/models/ode_model_pydeseq2_Obese_vs_Normal_v10_mean.py` (optional —
      only the Obese curve is skipped if this is missing)
- [x] `data/ratios/GSE48452_node_initial_ratios.xlsx`
- [x] `data/ratios/GSE89632_node_initial_ratios.xlsx`
- [x] `data/ratios/GSE130970_node_initial_ratios.xlsx`
- [x] `data/ratios/GSE162694_node_initial_ratios.xlsx`
- [x] `data/ratios/GSE213621_node_initial_ratios.xlsx`

All eight are copied automatically by the updated `organize_repo_files.py`
(source: `DROPBOX_ROOT/models_pydeseq2_mean/` for the three model files, and
the respective `*_validation_pipeline\` folders for the five ratio files).

## Explicitly NOT recommended to add (redundant / not needed)
- Top-level `step1_parse_series_matrix.py` … `step4_concordance.py` (the
  4-file version directly under `Independent Validation\`) — superseded by
  the more complete 7-file (step0–step6) versions inside each cohort's own
  `*_validation_pipeline\` subfolder, which include the model-output and
  silymarin-validation steps.
- Raw GEO data files (`*_series_matrix.txt[.gz]`, `*_expression_probe.csv`,
  `*_raw_counts.csv[.gz]`, `*_FPKMs_allsamples.txt[.gz]`,
  `*_salmon_tximport_*.csv[.gz]`) — publicly available from GEO; do not
  redistribute in this repository (see main README, Data availability).
- `__pycache__/` folders, regenerable heatmap/dynamics PNGs — build artifacts.
- `node_name_table_hsa04932.xlsx` copies inside each pipeline folder — a
  single copy is already in `data/`.
