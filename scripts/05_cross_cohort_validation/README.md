# MASLD ODE Model: A Transcriptomics-Constrained Systems Model of Non-Alcoholic Fatty Liver Disease

Code and analysis pipeline accompanying:

> Tseng, Y.-Y. Integrating transcriptomic big data with ODE modeling to
> identify regulatory mechanisms and an early-intervention window in MASLD.
> *npj Systems Biology and Applications* (submitted).

## Overview

This repository contains the code used to build a 70-node, 140-state
ordinary differential equation (ODE) model of the KEGG NAFLD signaling
pathway (hsa04932), constrained by liver transcriptomes from six independent
patient cohorts, and used to simulate disease-stage progression and *in
silico* silymarin intervention.

**Status note:** This repository skeleton was assembled programmatically from
the author's working files. See `COPY_MANIFEST.md` for the exact list of
files still to be copied in from the author's local environment (with
confirmed source paths). The topological-distance sensitivity analysis
(Supplementary Fig. S4) was reconstructed after the original script could
not be located; the reconstruction has been independently verified to
reproduce the original results exactly — see `COPY_MANIFEST.md` for details.

**Self-containment (2026-07-02):** the five scripts under
`scripts/05_cross_cohort_validation/` were rewritten so that all data
dependencies (the 3 ODE model files and 5 cohort ratio tables) are read from
`data/models/` and `data/ratios/` inside this repository, resolved via
`os.path.dirname(__file__)`. No external drive paths or personal folder
names need to be edited — once the files listed in `COPY_MANIFEST.md` are
copied in, this repository runs standalone.

## Repository structure

```
MASLD-ODE-Model/
├── COPY_MANIFEST.md                   Checklist of files still to be copied in
│                                       from the author's local environment,
│                                       with confirmed source paths. Delete
│                                       once all items are checked off.
├── data/                               Fully self-contained data:
│   ├── models/                         3 ODE model definitions (used by
│   │                                   the cross-cohort scripts)
│   ├── ratios/                         5 cohort node-initial-ratio tables
│   └── (root)                          small derived reference tables
│                                       (node↔gene mapping, initial ratios).
│                                       Raw GEO count matrices are NOT
│                                       included; see "Data availability"
│                                       below. No external folders outside
│                                       this repo are required to run any
│                                       of the analysis scripts.
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
