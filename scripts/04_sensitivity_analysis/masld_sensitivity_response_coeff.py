#!/usr/bin/env python3
"""
MASLD ODE — Initial-Concentration Response Coefficient Sensitivity Analysis
===========================================================================
Method (consistent with Tseng 2024, Nutrition & Metabolism, Eq. 3):

    R_i = (x0_i / O) * (dO / dx0_i)

i.e. the normalised local sensitivity of a pathway output O to a perturbation
of the initial concentration x0_i of each molecular node. After normalisation
x0_i cancels, so this equals  d(ln O) / d(ln x0_i)  and is computed by a
two-sided (central) finite difference:

    R_i = ( O(x0_i*(1+delta)) - O(x0_i*(1-delta)) ) / ( 2*delta*O_base )

Difference from N&M 2024:
  The three MASLD core outputs (P_Cell_death, P_Hepatocyte_injury,
  P_Inflammation) are SATURATING. Their steady-state levels converge to a
  common ceiling, so a steady-state response coefficient is near-zero and
  uninformative. We therefore use AUC (0-300 h, Simpson integration) as the
  output functional O — the same primary metric used throughout the
  manuscript. Steady-state-based coefficients are also reported (secondary
  columns) for completeness.

Perturbation target : initial INACTIVE-form concentration of each molecular
                      (gene) node — i.e. the value that carries the
                      transcriptomic expression ratio (active forms start at 0).
                      Pathway-output (P_*) nodes are NOT perturbed (they are the
                      readouts, matching N&M 2024 which perturbs enzymes and
                      measures the output).
Disease condition   : NASH by default (most severe state; analogous to the HFD
                      model used for the response coefficient in N&M 2024).
                      Add 'NAFL' to CONDITIONS to also run the steatosis state.

Outputs:
  - Sensitivity_ResponseCoefficients.xlsx  (one sheet per condition + metadata)
  - Fig_Sensitivity_tornado.{png,pdf}      (1x3 tornado, top-N nodes per output)

Run from the same directory as masld_fig5_simulation_v10_mean.py
(i.e. the parent of models_pydeseq2_mean/). Environment: deseq2_env.

Author: Yu-Yao Tseng
"""

import os
import sys
import importlib.util
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.integrate import odeint, simpson

# ============================================================
# CONFIG
# ============================================================
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))

MODEL_FILES = {
    "NAFL": "models_pydeseq2_mean/ode_model_pydeseq2_NAFL_vs_Normal_v10_mean.py",
    "NASH": "models_pydeseq2_mean/ode_model_pydeseq2_NASH_vs_Normal_v10_mean.py",
}

# Disease condition(s) on which to compute the response coefficient.
# Default: NASH only (matches the single-disease-condition choice of N&M 2024).
# To also analyse the steatosis stage, set: CONDITIONS = ["NASH", "NAFL"]
CONDITIONS = ["NASH"]

TARGET_OUTPUTS = ["P_Cell_death", "P_Hepatocyte_injury", "P_Inflammation"]

T_END     = 300        # h  (matches manuscript AUC window)
T_POINTS  = 7501
DELTA     = 0.01       # +/-1% relative perturbation for the central difference
TOP_N     = 15         # nodes shown per tornado panel

# Colour scheme: positive RC = output-increasing (red); negative = output-
# decreasing (blue). Red/blue is colour-blind-safe (unlike the red/green of
# N&M 2024 Fig. 3); this also satisfies the npj colour-accessibility guidance.
COL_POS = "#C0392B"
COL_NEG = "#2C6FB3"

plt.rcParams.update({
    "font.family":      "sans-serif",
    "font.sans-serif":  ["Calibri", "Arial", "DejaVu Sans"],
    "font.size":        7,
    "axes.labelsize":   7.5,
    "axes.titlesize":   8,
    "axes.linewidth":   0.5,
    "xtick.labelsize":  6.5,
    "ytick.labelsize":  6.5,
    "savefig.dpi":      300,
})


# ============================================================
# HELPERS
# ============================================================
def load_model(filepath):
    name = os.path.splitext(os.path.basename(filepath))[0]
    spec = importlib.util.spec_from_file_location(name, filepath)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def simulate(mod, y0, t):
    """Baseline kinetics (params=None -> model PARAMS), no drug."""
    return odeint(mod.ode_system, y0, t, args=(None,), mxstep=10000)


def output_metrics(y, out_idx, t, last_n):
    """Return dict of AUC and steady-state for each target output index."""
    auc = {name: float(simpson(y[:, i], x=t)) for name, i in out_idx.items()}
    ss  = {name: float(np.mean(y[-last_n:, i])) for name, i in out_idx.items()}
    return auc, ss


# ============================================================
# CORE: response coefficients for one condition
# ============================================================
def response_coefficients(mod, t):
    SV = mod.STATE_VARS
    idx = {s: k for k, s in enumerate(SV)}
    last_n = max(1, int(len(t) * 0.1))

    # Output active-form indices
    out_idx = {}
    for name in TARGET_OUTPUTS:
        key = name + "_active"
        if key not in idx:
            print(f"  WARNING: output {key} not found in STATE_VARS")
            continue
        out_idx[name] = idx[key]

    # Molecular (gene) nodes to perturb = inactive form of every GENE_METADATA node
    gene_nodes = list(mod.GENE_METADATA.keys())
    pert = []  # (node, symbol, inactive_index, x0)
    for node in gene_nodes:
        ina = node + "_inactive"
        if ina not in idx:
            continue
        j = idx[ina]
        x0 = float(mod.Y0[j])
        if x0 == 0.0:
            continue  # cannot take a relative perturbation of zero
        sym = mod.GENE_METADATA[node].get("symbol", node)
        pert.append((node, sym, j, x0))

    # Baseline
    y_base = simulate(mod, mod.Y0.copy(), t)
    auc_base, ss_base = output_metrics(y_base, out_idx, t, last_n)

    rows = []
    n = len(pert)
    print(f"  Perturbing {n} molecular nodes (+/-{DELTA*100:.0f}%) "
          f"against {len(out_idx)} outputs ...")
    for k, (node, sym, j, x0) in enumerate(pert, 1):
        # +delta
        yp = mod.Y0.copy(); yp[j] = x0 * (1.0 + DELTA)
        auc_p, ss_p = output_metrics(simulate(mod, yp, t), out_idx, t, last_n)
        # -delta
        ym = mod.Y0.copy(); ym[j] = x0 * (1.0 - DELTA)
        auc_m, ss_m = output_metrics(simulate(mod, ym, t), out_idx, t, last_n)

        row = {"Node": node, "Symbol": sym, "x0_ratio": round(x0, 4)}
        max_abs = 0.0
        for name in out_idx:
            rc_auc = ((auc_p[name] - auc_m[name]) /
                      (2.0 * DELTA * auc_base[name])) if auc_base[name] != 0 else 0.0
            rc_ss = ((ss_p[name] - ss_m[name]) /
                     (2.0 * DELTA * ss_base[name])) if ss_base[name] != 0 else 0.0
            row[f"RC_AUC_{name}"] = round(rc_auc, 5)
            row[f"RC_SS_{name}"]  = round(rc_ss, 5)
            max_abs = max(max_abs, abs(rc_auc))
        row["max_abs_RC_AUC"] = round(max_abs, 5)
        rows.append(row)
        if k % 10 == 0 or k == n:
            print(f"    {k}/{n}")

    rows.sort(key=lambda r: r["max_abs_RC_AUC"], reverse=True)
    return rows, out_idx, auc_base, ss_base


# ============================================================
# FIGURE
# ============================================================
def tornado_figure(results_by_cond, primary_cond):
    rows = results_by_cond[primary_cond]
    fig, axes = plt.subplots(1, 3, figsize=(7.205, 4.7))
    fig.subplots_adjust(top=0.90, bottom=0.16, left=0.115, right=0.985, wspace=0.55)

    for ci, name in enumerate(TARGET_OUTPUTS):
        ax = axes[ci]
        col = f"RC_AUC_{name}"
        ranked = sorted(rows, key=lambda r: abs(r.get(col, 0.0)), reverse=True)[:TOP_N]
        ranked = sorted(ranked, key=lambda r: r.get(col, 0.0))  # ascending for barh
        labels = [r["Symbol"] for r in ranked]
        vals   = [r.get(col, 0.0) for r in ranked]
        colors = [COL_POS if v >= 0 else COL_NEG for v in vals]
        ypos = np.arange(len(ranked))
        ax.barh(ypos, vals, color=colors, edgecolor="black", linewidth=0.3, height=0.72)
        ax.set_yticks(ypos)
        ax.set_yticklabels(labels, fontsize=6)
        ax.axvline(0, color="0.3", linewidth=0.5)
        ax.set_title(name, fontsize=8, fontweight="bold")
        ax.set_xlabel("Response coefficient (AUC)")
        ax.grid(True, axis="x", linestyle=":", color="0.85", linewidth=0.3)
        for s in ax.spines.values():
            s.set_linewidth(0.4)

    # legend
    from matplotlib.patches import Patch
    handles = [Patch(facecolor=COL_POS, edgecolor="black", linewidth=0.3,
                     label="Increases output (positive)"),
               Patch(facecolor=COL_NEG, edgecolor="black", linewidth=0.3,
                     label="Decreases output (negative)")]
    fig.legend(handles=handles, loc="lower center", ncol=2, frameon=False,
               fontsize=7.5, bbox_to_anchor=(0.5, 0.015))
    fig.suptitle(f"Initial-concentration response coefficients ({primary_cond} model)",
                 fontsize=9, fontweight="bold", y=0.97)

    for ext in ["png", "pdf"]:
        out = os.path.join(SCRIPT_DIR, f"Fig_Sensitivity_tornado.{ext}")
        fig.savefig(out, dpi=300 if ext == "png" else None,
                    bbox_inches="tight", pad_inches=0.08)
        print(f"  Saved: {out}")
    plt.close(fig)


# ============================================================
# MAIN
# ============================================================
def main():
    t = np.linspace(0, T_END, T_POINTS)
    results_by_cond = {}
    base_by_cond = {}

    for cond in CONDITIONS:
        fpath = os.path.join(SCRIPT_DIR, MODEL_FILES[cond])
        if not os.path.exists(fpath):
            print(f"ERROR: model not found: {fpath}")
            sys.exit(1)
        print(f"\n=== Condition: {cond} ===")
        mod = load_model(fpath)
        rows, out_idx, auc_base, ss_base = response_coefficients(mod, t)
        results_by_cond[cond] = rows
        base_by_cond[cond] = (auc_base, ss_base)

        # Console: top drivers per output
        print(f"\n  Top {min(8, len(rows))} drivers per output ({cond}, RC on AUC):")
        for name in TARGET_OUTPUTS:
            col = f"RC_AUC_{name}"
            top = sorted(rows, key=lambda r: abs(r.get(col, 0.0)), reverse=True)[:8]
            disp = ", ".join(f"{r['Symbol']}({r.get(col,0.0):+.3f})" for r in top)
            print(f"    {name}: {disp}")

    # Excel
    try:
        import pandas as pd
        xlsx = os.path.join(SCRIPT_DIR, "Sensitivity_ResponseCoefficients.xlsx")
        with pd.ExcelWriter(xlsx, engine="openpyxl") as writer:
            for cond in CONDITIONS:
                pd.DataFrame(results_by_cond[cond]).to_excel(
                    writer, sheet_name=f"RC_{cond}", index=False)
            meta_rows = [
                ["Script", os.path.basename(__file__)],
                ["Method", "Normalised initial-concentration response coefficient (Tseng 2024, Eq. 3)"],
                ["Definition", "R_i = d(ln O)/d(ln x0_i); central finite difference"],
                ["Output functional O", "AUC over 0-300 h (Simpson); steady-state (last 10%) reported as secondary"],
                ["Rationale (AUC vs SS)", "Core outputs are saturating; steady-state RC is near-zero/uninformative, so AUC (manuscript primary metric) is used"],
                ["Perturbation", f"+/-{DELTA*100:.0f}% of each molecular node inactive-form initial value"],
                ["Perturbed set", "GENE_METADATA nodes with nonzero initial ratio; P_* outputs excluded"],
                ["Conditions", ", ".join(CONDITIONS)],
                ["Outputs", ", ".join(TARGET_OUTPUTS)],
                ["Simulation", f"t=0..{T_END} h, {T_POINTS} points, odeint/LSODA, params=PARAMS"],
            ]
            for cond in CONDITIONS:
                ab, sb = base_by_cond[cond]
                for name in TARGET_OUTPUTS:
                    meta_rows.append([f"Baseline AUC [{cond}] {name}", round(ab.get(name, 0.0), 4)])
            pd.DataFrame(meta_rows, columns=["Parameter", "Value"]).to_excel(
                writer, sheet_name="Metadata", index=False)
        print(f"\n  Excel saved: {xlsx}")
    except ImportError:
        print("  (pandas/openpyxl not available; skipped Excel export)")

    # Figure (primary condition = first in CONDITIONS)
    print("\nGenerating tornado figure ...")
    tornado_figure(results_by_cond, CONDITIONS[0])
    print("\nDone.")


if __name__ == "__main__":
    main()
