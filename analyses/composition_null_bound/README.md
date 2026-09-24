# composition_null_bound

Composition-null bound for the stage-dependent predicted response
(Supplementary Table S17; response to Reviewer 3, Major Comment 4 and Reviewer 1, Major Comment 2).

## Question

Bulk liver transcriptomes mix hepatocyte-intrinsic expression with shifts in cell-type composition.
Seven model nodes are transcribed in substantial part by non-parenchymal cells:
TNF-α, FAS ligand, TGF-β1, IL-1, IL-6, CXCL8/IL-8 and CASP7.
Could composition alone produce the NAFL-over-NASH advantage?

## Design

Each of the seven nodes is re-initialized as `observed ratio ** f`, where `f` is the share of the
observed log₂ fold change assumed to be hepatocyte-intrinsic (`f = 0`: no intrinsic change).

| Scenario | Nodes | f |
|---|---|---|
| Joint sweep | all seven | 1, 0.75, 0.5, 0.25, 0.1, 0.05, 0.02, 0 |
| Single-node null | one at a time | 0 |

Everything else follows the standard protocol: kᵢ = 0.3 on the eight silymarin targets,
t = 0–300 h (3,001 points), `scipy.integrate.odeint` (LSODA), trapezoidal AUC of the active form.
The analysis is deterministic.

## Files

| File | Content |
|---|---|
| `composition_null_bound.py` | Main script |
| `masld_model.py` | Self-contained hsa04932 v10 Mean ODE model (same equations as `interaction_to_ode_keggid_v10_Mean.py`) |
| `inputs/result_v4_symbol_interactions.xlsx` | Network topology (106 interactions) |
| `inputs/node_name_table_hsa04932.xlsx` | Node → gene mapping (58 molecular nodes) |
| `inputs/DEG_pydeseq2_3comparisons.xlsx` | pyDESeq2 ratios, GSE126848 |
| `outputs/SuppTable_S17_composition_null_bound.xlsx` | Created on run |

## Run

```
conda activate deseq2_env
python composition_null_bound.py
```

Requires Python ≥ 3.9, numpy, scipy, pandas, openpyxl. Runtime is about 1–2 minutes.

The script first reproduces the primary analysis and prints `Baseline check: OK` when the
advantage ratios match the manuscript (P_Cell_death 2.443, P_Hepatocyte_injury 2.105,
P_Inflammation 3.034).

## Expected result

| Scenario | P_Cell_death | P_Hepatocyte_injury | P_Inflammation |
|---|---|---|---|
| Primary analysis (f = 1) | 2.443 | 2.105 | 3.034 |
| Joint, f = 0.25 | 1.684 | 1.094 | 3.088 |
| Joint composition-null (f = 0) | 1.086 | 1.017 | 1.000 |
| Single-node null: FAS ligand | 0.984 | 2.101 | 3.034 |
| Single-node null: CXCL8/IL-8 | 2.443 | 2.105 | 0.646 |
