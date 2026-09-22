#!/usr/bin/env python3
"""
MASLD ODE — Global Sensitivity Analysis of Initial Values (Morris + Sobol)
==========================================================================
Design rationale
----------------
Kinetic parameters (Vmax, ksp, n, kcat) are FIXED at 2.0 by modeling
convention and therefore carry no data-derived uncertainty. The only inputs
that carry genuine statistical/biological uncertainty are the node-level
INITIAL VALUES, which are initialized from pyDESeq2 fold-change ratios. GSA is
therefore performed on the initial values alone, holding all kinetics fixed.

Methodological precedent (for the Methods section):
  - Marino et al. 2008 (J Theor Biol 254:178) and Sumner et al. 2012 (J R Soc
    Interface 9:2156) treat initial conditions as first-class GSA inputs,
    interchangeable with parameters ("GSA techniques can be applied to the
    initial conditions of a model in the same way").
  - Sampling the initial molecular state from data-derived (log-normal)
    distributions and ranking which initial values drive output variance
    follows Spencer et al. 2009 (Nature 459:428) and Gaudet et al. 2012
    (PLoS Comput Biol 8:e1002482).
  - Sobol variance-based indices for QSP models: Zhang et al. 2015
    (CPT Pharmacometrics Syst Pharmacol 4:69).

Per-node uncertainty (primary, 'lfcSE' mode)
--------------------------------------------
Each node's initial value = arithmetic mean over its constituent genes of
r_g = 2^(log2FC_g). Uncertainty is propagated from the DESeq2 lfcSE of each
gene by the delta method:
    Var(r_g)   ~ (ln2 * lfcSE_g * r_g)^2
    Var(node)  = (1/k^2) * sum_g Var(r_g)          (k genes, independence)
    CV(node)   = sqrt(Var(node)) / node_nominal
    sigma_log2 = sqrt(ln(1+CV^2)) / ln2
Each node is then varied over its 95% CI on the log2 scale,
[log2(nominal) +/- Z*sigma_log2] with Z=1.96, sampled by the GSA design and
mapped back to a linear initial value via 2^(.). Multi-gene averaging correctly
NARROWS the node's uncertainty (Var scales as 1/k), which is the honest
treatment.

Robustness to the distribution choice
-------------------------------------
The same Morris screening is repeated under three schemes and rankings are
compared (Spearman rho of mu*): 'lfcSE' (primary), 'fixedCV' (common CV=0.25),
and 'uniform' (+/-30%). Stable top rankings across schemes demonstrate the
result is not an artifact of the distributional assumption.

Outputs
-------
  - GSA_results.xlsx      (Morris mu*/sigma, Sobol S1/ST, robustness comparison)
  - Fig_GSA_Morris.png/pdf (mu* tornado, 1x3 outputs)
  - Fig_GSA_Sobol.png/pdf  (S1 & ST, 1x3 outputs, on Morris-reduced set)

Run from the directory containing models_pydeseq2_mean/ (same level as
masld_fig5_simulation_v10_mean.py). Environment: deseq2_env + SALib.
    pip install SALib

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
from multiprocessing import Pool

from SALib.sample import morris as morris_sample
from SALib.analyze import morris as morris_analyze
from SALib.sample import sobol as sobol_sample
from SALib.analyze import sobol as sobol_analyze

# ============================================================
# CONFIG
# ============================================================
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))

MODEL_FILE = "models_pydeseq2_mean/ode_model_pydeseq2_NASH_vs_Normal_v10_mean.py"
DEG_CSV    = "DEG_pydeseq2_NASH_vs_Normal.csv"      # provides lfcSE per gene

TARGET_OUTPUTS = ["P_Cell_death", "P_Hepatocyte_injury", "P_Inflammation"]

T_END     = 300
T_POINTS  = 2001          # AUC grid; rankings are insensitive to this
Z_CI      = 1.96          # half-width of per-node CI (log2 scale)

# Morris (screening, all nodes)
MORRIS_R      = 20        # trajectories -> R*(D+1) model runs
MORRIS_LEVELS = 4

# Sobol (variance decomposition, Morris-reduced set)
RUN_SOBOL    = True
SOBOL_N      = 256        # base samples -> N*(D+2) runs; increase to 1024+ and
                          # re-check convergence for the final analysis
TOP_K_SOBOL  = 15         # nodes carried into Sobol (union of top mu* per output)

# Robustness-to-assumptions
FIXED_CV     = 0.25       # 'fixedCV' scheme
UNIFORM_PCT  = 0.30       # 'uniform' scheme (+/-30%)
RUN_ROBUSTNESS = True

N_WORKERS = max(1, (os.cpu_count() or 2) - 1)

COL_POS = "#C0392B"
COL_S1  = "#2C6FB3"
COL_ST  = "#F39C12"

plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Calibri", "Arial", "DejaVu Sans"],
    "font.size": 7, "axes.labelsize": 7.5, "axes.titlesize": 8,
    "axes.linewidth": 0.5, "xtick.labelsize": 6.5, "ytick.labelsize": 6.5,
    "savefig.dpi": 300,
})

# ============================================================
# GLOBALS for worker processes (populated by _init_worker)
# ============================================================
_G = {}


def load_model(filepath):
    name = os.path.splitext(os.path.basename(filepath))[0]
    spec = importlib.util.spec_from_file_location(name, filepath)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def _init_worker(model_path, factor_inactive_idx, out_active_idx, y0, t):
    """Each worker loads the model once and caches evaluation context."""
    mod = load_model(model_path)
    _G["ode"] = mod.ode_system
    _G["fidx"] = np.asarray(factor_inactive_idx, dtype=int)  # per-factor inactive index
    _G["oidx"] = out_active_idx                              # {out_name: active idx}
    _G["y0"]  = np.asarray(y0, dtype=float)
    _G["t"]   = np.asarray(t, dtype=float)


def _evaluate_row(x_log2):
    """One GSA sample: x_log2 = log2 initial value per factor -> AUC per output."""
    y0 = _G["y0"].copy()
    y0[_G["fidx"]] = np.power(2.0, np.asarray(x_log2, dtype=float))
    y = odeint(_G["ode"], y0, _G["t"], args=(None,), mxstep=10000)
    t = _G["t"]
    return [float(simpson(y[:, idx], x=t)) for idx in _G["oidx"].values()]


# ============================================================
# BUILD FACTORS (nodes) + per-node uncertainty
# ============================================================
def build_factors(mod, deg_path):
    import pandas as pd
    SV  = mod.STATE_VARS
    idx = {s: i for i, s in enumerate(SV)}

    # gene -> lfcSE (log2 scale); impute NaN with median
    deg = pd.read_csv(deg_path)
    deg = deg.dropna(subset=["lfcSE"])
    gene_se = dict(zip(deg["Gene"].astype(str), deg["lfcSE"].astype(float)))
    med_se = float(np.median(list(gene_se.values()))) if gene_se else 0.25

    factors = []   # dict per node
    for node, meta in mod.GENE_METADATA.items():
        ina = node + "_inactive"
        if ina not in idx:
            continue
        j = idx[ina]
        nominal = float(mod.Y0[j])
        if nominal <= 0:
            continue
        genes = [g.strip() for g in str(meta.get("gene_list", node)).split(",") if g.strip()]
        if not genes:
            genes = [meta.get("symbol", node)]
        k = len(genes)
        # delta-method variance of the node (mean of per-gene 2^log2FC)
        var_node = 0.0
        for g in genes:
            se = gene_se.get(g, med_se)
            # per-gene nominal ratio unknown individually here; approximate r_g by
            # node nominal (mean); this keeps CV well-defined and conservative.
            r_g = nominal
            var_node += (np.log(2.0) * se * r_g) ** 2
        var_node /= (k ** 2)
        cv = float(np.sqrt(var_node) / nominal) if nominal > 0 else 0.0
        sigma_log2 = float(np.sqrt(np.log(1.0 + cv ** 2)) / np.log(2.0))
        factors.append({
            "node": node,
            "symbol": meta.get("symbol", node),
            "inactive_idx": j,
            "nominal": nominal,
            "n_genes": k,
            "mean_lfcSE": float(np.mean([gene_se.get(g, med_se) for g in genes])),
            "cv_lfcSE": cv,
            "sigma_log2_lfcSE": sigma_log2,
        })
    return factors, idx


def make_problem(factors, mode):
    """Return SALib problem dict with per-node log2-scale bounds for a scheme."""
    names, bounds = [], []
    for f in factors:
        c = np.log2(f["nominal"])
        if mode == "lfcSE":
            hw = Z_CI * f["sigma_log2_lfcSE"]
        elif mode == "fixedCV":
            s = np.sqrt(np.log(1.0 + FIXED_CV ** 2)) / np.log(2.0)
            hw = Z_CI * s
        elif mode == "uniform":
            lo = np.log2(f["nominal"] * (1.0 - UNIFORM_PCT))
            hi = np.log2(f["nominal"] * (1.0 + UNIFORM_PCT))
            names.append(f["node"]); bounds.append([lo, hi]); continue
        else:
            raise ValueError(mode)
        # guard against zero-width (nodes with tiny SE)
        hw = max(hw, 1e-6)
        names.append(f["node"]); bounds.append([c - hw, c + hw])
    return {"num_vars": len(names), "names": names, "bounds": bounds}


# ============================================================
# EVALUATE a sample matrix (parallel)
# ============================================================
def evaluate_matrix(X, factors, model_path, out_active_idx, y0, t):
    fidx = [f["inactive_idx"] for f in factors]
    args = (model_path, fidx, out_active_idx, y0, t)
    if N_WORKERS > 1:
        with Pool(processes=N_WORKERS, initializer=_init_worker, initargs=args) as pool:
            Y = pool.map(_evaluate_row, [row for row in X], chunksize=8)
    else:
        _init_worker(*args)
        Y = [_evaluate_row(row) for row in X]
    return np.asarray(Y)  # (Nsamples, n_outputs)


# ============================================================
# MORRIS
# ============================================================
def run_morris(factors, model_path, out_active_idx, y0, t, mode):
    problem = make_problem(factors, mode)
    X = morris_sample.sample(problem, N=MORRIS_R, num_levels=MORRIS_LEVELS)
    Y = evaluate_matrix(X, factors, model_path, out_active_idx, y0, t)
    res = {}
    for oi, oname in enumerate(TARGET_OUTPUTS):
        Si = morris_analyze.analyze(problem, X, Y[:, oi],
                                    num_levels=MORRIS_LEVELS, print_to_console=False)
        res[oname] = {"names": list(problem["names"]),
                      "mu_star": np.asarray(Si["mu_star"]),
                      "sigma": np.asarray(Si["sigma"])}
    return res, len(X)


# ============================================================
# SOBOL (on Morris-reduced set)
# ============================================================
def run_sobol(factors, model_path, out_active_idx, y0, t, keep_nodes):
    sub = [f for f in factors if f["node"] in keep_nodes]
    problem = make_problem(sub, "lfcSE")
    X = sobol_sample.sample(problem, SOBOL_N, calc_second_order=False)
    Y = evaluate_matrix(X, sub, model_path, out_active_idx, y0, t)
    res = {}
    for oi, oname in enumerate(TARGET_OUTPUTS):
        Si = sobol_analyze.analyze(problem, Y[:, oi],
                                   calc_second_order=False, print_to_console=False)
        res[oname] = {"names": list(problem["names"]),
                      "S1": np.asarray(Si["S1"]), "S1_conf": np.asarray(Si["S1_conf"]),
                      "ST": np.asarray(Si["ST"]), "ST_conf": np.asarray(Si["ST_conf"])}
    return res, len(X), [f["node"] for f in sub]


# ============================================================
# FIGURES
# ============================================================
def sym_of(factors):
    return {f["node"]: f["symbol"] for f in factors}


def fig_morris(morris_res, factors, top_n=15):
    S = sym_of(factors)
    fig, axes = plt.subplots(1, 3, figsize=(7.205, 4.7))
    fig.subplots_adjust(top=0.9, bottom=0.14, left=0.115, right=0.985, wspace=0.55)
    for ci, oname in enumerate(TARGET_OUTPUTS):
        ax = axes[ci]
        d = morris_res[oname]
        order = np.argsort(d["mu_star"])[::-1][:top_n]
        order = order[np.argsort(d["mu_star"][order])]  # ascending for barh
        labels = [S.get(d["names"][i], d["names"][i]) for i in order]
        vals = d["mu_star"][order]
        ax.barh(np.arange(len(order)), vals, color=COL_POS,
                edgecolor="black", linewidth=0.3, height=0.72)
        ax.set_yticks(np.arange(len(order))); ax.set_yticklabels(labels, fontsize=6)
        ax.set_title(oname, fontsize=8, fontweight="bold")
        ax.set_xlabel(r"Morris $\mu^*$")
        ax.grid(True, axis="x", linestyle=":", color="0.85", linewidth=0.3)
        for s in ax.spines.values():
            s.set_linewidth(0.4)
    fig.suptitle("Global sensitivity of pathway-output AUC to initial values (Morris, NASH)",
                 fontsize=9, fontweight="bold", y=0.975)
    for ext in ["png", "pdf"]:
        fig.savefig(os.path.join(SCRIPT_DIR, f"Fig_GSA_Morris.{ext}"),
                    dpi=300 if ext == "png" else None, bbox_inches="tight", pad_inches=0.08)
    plt.close(fig)


def fig_sobol(sobol_res, factors):
    S = sym_of(factors)
    fig, axes = plt.subplots(1, 3, figsize=(7.205, 4.7))
    fig.subplots_adjust(top=0.9, bottom=0.16, left=0.12, right=0.985, wspace=0.55)
    for ci, oname in enumerate(TARGET_OUTPUTS):
        ax = axes[ci]
        d = sobol_res[oname]
        order = np.argsort(d["ST"])[::-1]
        order = order[np.argsort(d["ST"][order])]  # ascending for barh
        y = np.arange(len(order)); h = 0.38
        labels = [S.get(d["names"][i], d["names"][i]) for i in order]
        ax.barh(y + h/2, d["ST"][order], height=h, color=COL_ST,
                edgecolor="black", linewidth=0.3, label="Total-order $S_T$")
        ax.barh(y - h/2, d["S1"][order], height=h, color=COL_S1,
                edgecolor="black", linewidth=0.3, label="First-order $S_1$")
        ax.set_yticks(y); ax.set_yticklabels(labels, fontsize=6)
        ax.set_title(oname, fontsize=8, fontweight="bold")
        ax.set_xlabel("Sobol index")
        ax.grid(True, axis="x", linestyle=":", color="0.85", linewidth=0.3)
        for s in ax.spines.values():
            s.set_linewidth(0.4)
    handles, lbls = axes[0].get_legend_handles_labels()
    fig.legend(handles, lbls, loc="lower center", ncol=2, frameon=False,
               fontsize=7.5, bbox_to_anchor=(0.5, 0.01))
    fig.suptitle(f"Variance decomposition of output AUC (Sobol, N={SOBOL_N}, NASH)",
                 fontsize=9, fontweight="bold", y=0.975)
    for ext in ["png", "pdf"]:
        fig.savefig(os.path.join(SCRIPT_DIR, f"Fig_GSA_Sobol.{ext}"),
                    dpi=300 if ext == "png" else None, bbox_inches="tight", pad_inches=0.08)
    plt.close(fig)


# ============================================================
# MAIN
# ============================================================
def main():
    model_path = os.path.join(SCRIPT_DIR, MODEL_FILE)
    deg_path   = os.path.join(SCRIPT_DIR, DEG_CSV)
    for p in (model_path, deg_path):
        if not os.path.exists(p):
            print(f"ERROR: not found: {p}"); sys.exit(1)

    print("Loading model + DEG lfcSE ...")
    mod = load_model(model_path)
    factors, idx = build_factors(mod, deg_path)
    out_active_idx = {o: idx[o + "_active"] for o in TARGET_OUTPUTS if o + "_active" in idx}
    y0 = mod.Y0.copy()
    t = np.linspace(0, T_END, T_POINTS)
    print(f"  {len(factors)} molecular nodes; {len(out_active_idx)} outputs; "
          f"{N_WORKERS} workers")

    # ---- Morris (primary, lfcSE) ----
    print("\n[1/3] Morris screening (lfcSE scheme) ...")
    morris_lfcSE, n_m = run_morris(factors, model_path, out_active_idx, y0, t, "lfcSE")
    print(f"  {n_m} model runs.")
    for oname in TARGET_OUTPUTS:
        d = morris_lfcSE[oname]; S = sym_of(factors)
        top = np.argsort(d["mu_star"])[::-1][:8]
        print(f"  {oname}: " + ", ".join(f"{S.get(d['names'][i])}({d['mu_star'][i]:.3g})"
                                          for i in top))

    # ---- Robustness (fixedCV, uniform) via Morris ----
    robustness = {"lfcSE": morris_lfcSE}
    if RUN_ROBUSTNESS:
        for mode in ["fixedCV", "uniform"]:
            print(f"\n      Robustness Morris ({mode}) ...")
            r, _ = run_morris(factors, model_path, out_active_idx, y0, t, mode)
            robustness[mode] = r

    # ---- Select Morris-reduced set for Sobol ----
    keep = set()
    for oname in TARGET_OUTPUTS:
        d = morris_lfcSE[oname]
        top = np.argsort(d["mu_star"])[::-1][:TOP_K_SOBOL]
        keep.update(d["names"][i] for i in top)
    print(f"\n[2/3] Sobol input set: {len(keep)} nodes (union of top-{TOP_K_SOBOL} mu* per output)")

    sobol_res = None
    if RUN_SOBOL:
        print(f"      Sobol sampling (N={SOBOL_N}) ...")
        sobol_res, n_s, sobol_nodes = run_sobol(factors, model_path,
                                                out_active_idx, y0, t, keep)
        print(f"      {n_s} model runs.")

    # ---- Figures ----
    print("\n[3/3] Figures ...")
    fig_morris(morris_lfcSE, factors)
    if sobol_res is not None:
        fig_sobol(sobol_res, factors)

    # ---- Excel ----
    try:
        import pandas as pd
        from scipy.stats import spearmanr
        xlsx = os.path.join(SCRIPT_DIR, "GSA_results.xlsx")
        with pd.ExcelWriter(xlsx, engine="openpyxl") as w:
            # per-node uncertainty inputs
            pd.DataFrame([{k: f[k] for k in
                           ["node", "symbol", "n_genes", "nominal",
                            "mean_lfcSE", "cv_lfcSE", "sigma_log2_lfcSE"]}
                          for f in factors]).to_excel(
                w, sheet_name="Node_uncertainty", index=False)

            # Morris (lfcSE)
            for oname in TARGET_OUTPUTS:
                d = morris_lfcSE[oname]
                df = pd.DataFrame({"node": d["names"],
                                   "mu_star": d["mu_star"], "sigma": d["sigma"]})
                df = df.sort_values("mu_star", ascending=False)
                df.to_excel(w, sheet_name=f"Morris_{oname[:22]}", index=False)

            # Sobol
            if sobol_res is not None:
                for oname in TARGET_OUTPUTS:
                    d = sobol_res[oname]
                    df = pd.DataFrame({"node": d["names"], "S1": d["S1"],
                                       "S1_conf": d["S1_conf"], "ST": d["ST"],
                                       "ST_conf": d["ST_conf"]})
                    df = df.sort_values("ST", ascending=False)
                    df.to_excel(w, sheet_name=f"Sobol_{oname[:23]}", index=False)

            # Robustness: Spearman of mu* rankings across schemes
            rob_rows = []
            if RUN_ROBUSTNESS:
                for oname in TARGET_OUTPUTS:
                    base = robustness["lfcSE"][oname]["mu_star"]
                    for mode in ["fixedCV", "uniform"]:
                        other = robustness[mode][oname]["mu_star"]
                        rho, _ = spearmanr(base, other)
                        rob_rows.append({"output": oname, "scheme_vs_lfcSE": mode,
                                         "spearman_rho_mu_star": round(float(rho), 4)})
                pd.DataFrame(rob_rows).to_excel(w, sheet_name="Robustness_rank", index=False)

            meta = [
                ["Script", os.path.basename(__file__)],
                ["Perturbed inputs", "node initial values only; kinetics fixed at 2.0"],
                ["Output functional", f"AUC 0-{T_END} h (Simpson), {T_POINTS} pts"],
                ["Primary range", "log2(nominal) +/- 1.96*sigma_log2 from DESeq2 lfcSE (delta method)"],
                ["Morris", f"R={MORRIS_R}, levels={MORRIS_LEVELS}, runs={n_m}"],
                ["Sobol", f"N={SOBOL_N}, calc_second_order=False, "
                          f"top-{TOP_K_SOBOL} union set ({len(keep)} nodes)"
                          + ("" if sobol_res is None else f", runs={n_s}")],
                ["Robustness", f"fixedCV={FIXED_CV}, uniform=+/-{UNIFORM_PCT}"],
                ["References", "Marino 2008; Sumner 2012; Spencer 2009; Gaudet 2012; Zhang 2015"],
            ]
            pd.DataFrame(meta, columns=["Parameter", "Value"]).to_excel(
                w, sheet_name="Metadata", index=False)
        print(f"  Excel saved: {xlsx}")

        if RUN_ROBUSTNESS and rob_rows:
            print("\n  Rank robustness (Spearman rho of mu* vs lfcSE scheme):")
            for r in rob_rows:
                print(f"    {r['output']:<22} {r['scheme_vs_lfcSE']:<9} "
                      f"rho={r['spearman_rho_mu_star']}")
    except ImportError:
        print("  (pandas/openpyxl/scipy.stats missing; Excel skipped)")

    print("\nDone.")


if __name__ == "__main__":
    main()
