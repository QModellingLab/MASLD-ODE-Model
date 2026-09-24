#!/usr/bin/env python3
"""
composition_null_bound.py
=========================
Composition-null bound for nodes whose bulk-liver transcripts may derive in
substantial part from non-parenchymal cells (primary cohort GSE126848).

Question
--------
Can shifts in cell-type composition alone produce the stage-dependent
predicted response (NAFL/NASH advantage ratio > 1)?

Design
------
Seven nodes (TNF-a, FAS ligand, TGF-b1, IL-1, IL-6, CXCL8/IL-8, CASP7) are
re-initialized assuming that only a fraction f of each node's observed log2
fold change is hepatocyte-intrinsic; the remainder is attributed to
composition:

    initial ratio = observed ratio ** f        (f = 0 -> no intrinsic change)

  1. Joint sweep:   all seven nodes, f = 1, 0.75, 0.5, 0.25, 0.1, 0.05, 0.02, 0
  2. Single-node:   f = 0 for one node at a time, others at observed ratios

All other settings follow the standard protocol: k_i = 0.3 on the eight
silymarin targets, t = 0-300 h (3,001 points), LSODA (odeint), trapezoidal AUC
of the active-form trajectory. The model is deterministic (no random seed).

Output
------
outputs/SuppTable_S17_composition_null_bound.xlsx
    Sheet 'Advantage_ratio' : one row per scenario (Supplementary Table S17)
    Sheet 'Long_format'     : per-stage reductions for every scenario/output
    Sheet 'Initial_ratios'  : observed NAFL/NASH ratios of the seven nodes
    Sheet 'Metadata'        : settings

Usage
-----
    python composition_null_bound.py
(or open in Spyder and press Run; paths are resolved relative to this file)
"""

import os
import sys
import datetime
import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from masld_model import MasldModel, stage_advantage, CORE_OUTPUTS  # noqa: E402

# =============================================================================
# Configuration
# =============================================================================
INPUT_DIR = os.path.join(HERE, "inputs")
OUTPUT_DIR = os.path.join(HERE, "outputs")
INTERACTIONS_FILE = os.path.join(INPUT_DIR, "result_v4_symbol_interactions.xlsx")
NODE_TABLE_FILE = os.path.join(INPUT_DIR, "node_name_table_hsa04932.xlsx")
DEG_FILE = os.path.join(INPUT_DIR, "DEG_pydeseq2_3comparisons.xlsx")
OUTPUT_FILE = os.path.join(OUTPUT_DIR, "SuppTable_S17_composition_null_bound.xlsx")

KI = 0.3
NPC_NODES = ["TNFa", "FasL", "TGF_b1", "IL_1", "IL_6", "IL_8", "CASP7"]
NODE_LABEL = {"TNFa": "TNF-α", "FasL": "FAS ligand", "TGF_b1": "TGF-β1",
              "IL_1": "IL-1", "IL_6": "IL-6", "IL_8": "CXCL8/IL-8", "CASP7": "CASP7"}
JOINT_FRACTIONS = [0.75, 0.50, 0.25, 0.10, 0.05, 0.02, 0.0]

# Reference values of the primary analysis (manuscript Table 4 / Table S16)
EXPECTED_BASELINE = {"P_Cell_death": 2.443, "P_Hepatocyte_injury": 2.105,
                     "P_Inflammation": 3.034}


# =============================================================================
# Helpers
# =============================================================================
def scaled_initial_states(model, nodes, fraction):
    """Return (Y0_NAFL, Y0_NASH) with the chosen nodes' ratios raised to `fraction`."""
    y_nafl = model.initial_state("NAFL")
    y_nash = model.initial_state("NASH")
    for node in nodes:
        i = model.inactive_index(node)
        y_nafl[i] = y_nafl[i] ** fraction
        y_nash[i] = y_nash[i] ** fraction
    return y_nafl, y_nash


def main():
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    for f in (INTERACTIONS_FILE, NODE_TABLE_FILE, DEG_FILE):
        if not os.path.exists(f):
            sys.exit(f"[ERROR] Input not found: {f}")

    model = MasldModel(INTERACTIONS_FILE, NODE_TABLE_FILE, DEG_FILE)
    print(f"Model: {len(model.genes)} molecular nodes, {len(model.pathways)} pathway outputs, "
          f"{2 * len(model.nodes)} state variables")
    missing = [n for n in NPC_NODES if n not in model.idx]
    if missing:
        sys.exit(f"[ERROR] Nodes not in model: {missing}")

    scenarios = [("Primary analysis (f = 1)", NPC_NODES, 1.0)]
    scenarios += [(f"Joint, f = {f:.2f}", NPC_NODES, f) for f in JOINT_FRACTIONS[:-1]]
    scenarios += [("Joint composition-null (f = 0)", NPC_NODES, 0.0)]
    scenarios += [(f"Single-node null: {NODE_LABEL[n]}", [n], 0.0) for n in NPC_NODES]

    wide, long_rows = [], []
    for label, nodes, frac in scenarios:
        y_nafl, y_nash = scaled_initial_states(model, nodes, frac)
        res = stage_advantage(model, y_nafl, y_nash, ki=KI)
        row = {"Scenario": label}
        for o in CORE_OUTPUTS:
            row[o] = round(res[o]["ratio"], 3)
            long_rows.append({"Scenario": label, "Nodes": "+".join(nodes), "f": frac,
                              "Output": o,
                              "NAFL_reduction_pct": res[o]["NAFL"],
                              "NASH_reduction_pct": res[o]["NASH"],
                              "Advantage_ratio": res[o]["ratio"]})
        wide.append(row)
        print(f"  {label:<38s} " + "  ".join(f"{row[o]:7.3f}" for o in CORE_OUTPUTS))

    # --- reproducibility check against the primary analysis ---
    base = wide[0]
    for o, ref in EXPECTED_BASELINE.items():
        if abs(base[o] - ref) > 0.002:
            print(f"[WARN] Baseline {o} = {base[o]:.3f} differs from manuscript value {ref:.3f}")
    print("Baseline check:", "OK" if all(abs(base[o] - r) <= 0.002
                                         for o, r in EXPECTED_BASELINE.items()) else "MISMATCH")

    init = pd.DataFrame([{"Node": NODE_LABEL[n],
                          "NAFL_ratio": model.initial_state("NAFL")[model.inactive_index(n)],
                          "NASH_ratio": model.initial_state("NASH")[model.inactive_index(n)]}
                         for n in NPC_NODES])
    meta = pd.DataFrame([
        ["Script", os.path.basename(__file__)],
        ["Run date", datetime.datetime.now().strftime("%Y-%m-%d %H:%M")],
        ["Cohort", "GSE126848 (primary), pyDESeq2 ratios, v10 Mean aggregation"],
        ["Nodes", ", ".join(NODE_LABEL[n] for n in NPC_NODES)],
        ["Attribution", "initial ratio = observed ratio ** f (f = hepatocyte-intrinsic share of log2FC)"],
        ["Intervention", f"k_i = {KI} on CASP3, CASP7, CASP8, CYP2E1, IL-8, TNF-a, NF-kB, TGF-b1"],
        ["Simulation", "t = 0-300 h, 3,001 points, scipy.integrate.odeint (LSODA), mxstep = 10,000"],
        ["Metric", "trapezoidal AUC of active form; ratio = NAFL % reduction / NASH % reduction"],
        ["Kinetics", "Vmax = kcat = ksp = 2.0, n = 2.0"],
    ], columns=["Parameter", "Value"])

    with pd.ExcelWriter(OUTPUT_FILE, engine="openpyxl") as w:
        pd.DataFrame(wide).to_excel(w, sheet_name="Advantage_ratio", index=False)
        pd.DataFrame(long_rows).to_excel(w, sheet_name="Long_format", index=False)
        init.to_excel(w, sheet_name="Initial_ratios", index=False)
        meta.to_excel(w, sheet_name="Metadata", index=False)
    print(f"Saved: {OUTPUT_FILE}")


if __name__ == "__main__":
    main()
