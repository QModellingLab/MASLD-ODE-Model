# MASLD ODE Model: A Transcriptomics-Constrained Systems Model of Non-Alcoholic Fatty Liver Disease

Code and analysis pipeline accompanying:

> Tseng, Y.-Y. A transcriptomics-constrained ODE model of MASLD signaling
> identifies injury drivers and a reproducible stage-dependent predicted
> response to intervention.
> *npj Systems Biology and Applications* (under revision).

## Overview

This repository contains the code used to build a 70-node, 140-state
ordinary differential equation (ODE) model of the KEGG NAFLD signaling
pathway (hsa04932), constrained by liver transcriptomes from six independent
patient cohorts, and used to simulate disease-stage progression and *in
silico* silymarin intervention.

## Revision analyses (`analyses/`)

Eleven control analyses added during peer review. Each folder is self-contained
(`README.txt`, script, and the outputs it produces) and is runnable on its own.

| Folder | Addresses | Supplementary table |
|---|---|---|
| `patient_bootstrap_resampling/` | Inter-patient variability (200 stratified bootstrap replicates) | S11 |
| `alt_node_aggregation/` | Four multi-gene node aggregation rules | S7 |
| `alt_initial_conditions/` | 3 x 3 grid of active-fraction and output-baseline assumptions | S8 |
| `auc_window_sensitivity/` | AUC windows of 100, 300, 600 and 1000 h | S9 |
| `joint_n_ksp_grid/` | Joint 3 x 3 Hill-coefficient / half-saturation grid | S4 |
| `global_sensitivity_analysis/` | Morris screening over all 58 nodes; topology comparison | — |
| `random_target_control/` | Random-target and matched-control null distributions (N = 2,000 each) | S10 |
| `leave_one_out_targets/` | Each silymarin target omitted in turn | S16 |
| `edge_scramble_control/` | Degree-preserving and degree-unconstrained topology nulls | S5 |
| `static_hill_curve_control/` | Algebraic steady state and single-reaction Hill reserve | S6 |
| `equilibrate_then_inhibit/` | Intervention applied from the pre-equilibrated state | S15 |
| `cohort_permutation_test/` | Cohort-level sign-flip and cohort-block permutation tests | S12 |

All analyses use fixed random seeds and are exactly reproducible. Raw per-replicate
records are included as `.jsonl` where the analysis is resampling-based, so the
reported statistics can be recomputed independently of the summary tables.

## Node naming

The fibrosis pathway-output node was spelled `P_Figrosis` in earlier versions of the
model code. It is spelled `P_Fibrosis` throughout this repository. The rename is
confined to the identifier; re-running the analyses before and after the change
reproduces every reported value exactly.

## Repository structure

```
MASLD-ODE-Model/
├── data/                              Derived reference tables (node↔gene mapping,
│                                       initial-condition ratios) and, under
│                                       data/GSE126848/, the primary-cohort count
│                                       matrices used by the analyses that refit
│                                       PyDESeq2. Other cohorts are downloaded from
│                                       GEO; see "Data availability" below.
├── scripts/
│   ├── 01_differential_expression/    PyDESeq2 analysis, volcano plot, GO/KEGG
│   │                                  enrichment dot plot, + journal-spec
│   │                                  companions (Figs. 2, 3)
│   ├── 02_ode_model_construction/     KEGG hsa04932 → ODE network translation,
│   │                                  network diagram + journal-spec Graphviz
│   │                                  rendering (Fig. 4)
│   ├── 03_simulation/                 Disease-stage dynamics, 12-output overview
│   │                                  (+ journal-spec companion), silymarin
│   │                                  intervention (Figs. 5, 6)
│   ├── 04_sensitivity_analysis/       Local response-coefficient sensitivity
│   │                                  analysis (Fig. 7) and an optional Morris
│   │                                  global sensitivity analysis script
│   └── 05_cross_cohort_validation/    Independent validation across 5 additional
│       ├── shared/                    cohorts (Table 3, Table 4, Fig. 9,
│       │                              Supplementary Figs. S2, S3, S5).
│       │                              node_ratio_core.py (gene-symbol/Entrez/
│       │                              Ensembl ID support, shared by all 3 new
│       │                              RNA-seq cohorts)
│       ├── GSE130970_GSE162694_GSE213621/  Per-cohort node-ratio computation
│       ├── GSE48452/                  Microarray pipeline (step0–step6)
│       ├── GSE89632/                  Microarray pipeline (step0–step6)
│       └── diagnostics/               Optional sample/gene-ID diagnostic
│                                      scripts (not required to reproduce
│                                      manuscript results)
└── figures/                           (empty; populated by the scripts above)
```

## Data availability

All transcriptomic datasets analyzed in this study are publicly available
from the NCBI Gene Expression Omnibus (GEO):

| Accession   | Reference             | Platform                                   |
|-------------|------------------------|---------------------------------------------|
| GSE126848   | Suppli et al., 2019    | RNA-seq (Illumina NextSeq 500)               |
| GSE48452    | Ahrens et al., 2013    | Microarray (Affymetrix Human Gene 1.1 ST)    |
| GSE89632    | Arendt et al., 2015    | Microarray (Illumina HumanHT-12 WG-DASL V4)  |
| GSE130970   | Hoang et al., 2019     | RNA-seq (Illumina HiSeq 2500)                |
| GSE162694   | Pantano et al., 2021   | RNA-seq (Illumina HiSeq 3000)                |
| GSE213621   | Chen et al., 2023      | RNA-seq (Illumina HiSeq 2500)                |

Raw count/expression matrices should be downloaded directly from GEO; they
are not redistributed in this repository.

## Environment

- Python 3.11
- Key packages: `pydeseq2`, `numpy`, `scipy`, `pandas`, `matplotlib`,
  `seaborn`, `openpyxl`, `networkx`, `graphviz`
- See `requirements.txt` for a full pinned list (fill in exact versions from
  your `deseq2_env` conda environment before publishing, e.g. via
  `pip freeze > requirements.txt`).

```bash
conda create -n deseq2_env python=3.11
conda activate deseq2_env
pip install -r requirements.txt
```

## Model summary

- **Pathway backbone:** KEGG hsa04932 (Non-alcoholic fatty liver disease)
- **Scale:** 70 nodes (58 molecular + 12 pathway-output), 140 state variables,
  140 reactions, 106 regulatory interactions
- **Kinetics:** Hill-function activation/inactivation, all rate and affinity
  constants fixed at 2.0
- **Initialization:** node-level initial conditions from PyDESeq2
  condition-vs-normal fold-change ratios (arithmetic mean across genes for
  multi-gene nodes)
- **Intervention:** silymarin modeled as a uniform inhibitory multiplier on
  eight target nodes (CASP3, CASP7, CASP8, CYP2E1, CXCL8/IL-8, TNF-α, NF-κB,
  TGF-β1)

## Citation

If you use this code, please cite the manuscript above. A Zenodo DOI will be
minted for the version-of-record upon publication.

## Contact

Yu-Yao Tseng (York) — york@g2.usc.edu.tw
Department of Food Science, Nutrition, and Nutraceutical Biotechnology,
Shih Chien University, Taipei, Taiwan
ORCID: 0000-0002-9944-9953

## License

Code: MIT License (suggested — replace with your preferred license before
publishing; see `LICENSE`).
