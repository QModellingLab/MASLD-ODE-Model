#!/usr/bin/env python3
"""
MASLD Fig 5 (v10 Mean, 2×3) — Upstream molecular vs Downstream pathway
========================================================================
Direct topological upstream→downstream pairs in the ODE model:
  CASP7  → P_Hepatocyte_injury  (via activation edge)
  TNF    → P_Cell_death          (via activation edge)
  CXCL8  → P_Inflammation        (via activation edge)

Panels (2×3):
  Row 1 — Upstream molecular nodes (single-gene nodes, no aggregation)
    (a) CASP7   active   (NS/N AUC = 3.50;  Normal<Obese<NAFL<NASH)
    (b) TNF-α   active   (NS/N AUC = 1.62;  NAFL~NASH, both >> Normal)
    (c) CXCL8   active   (NS/N AUC = 10.49; dramatic NASH rise)
  Row 2 — Downstream pathway outputs (AUC-based selection)
    (d) P_Hepatocyte_injury   (NS/N AUC = 24.85 ↑↑)
    (e) P_Cell_death          (NS/N AUC = 1.45  ↑)
    (f) P_Inflammation        (NS/N AUC = 1.87  ↑)

Evaluation: AUC over 100 h as the primary metric for pathway activation,
consistent with the interpretation that MASLD severity reflects
cumulative pathway activity over time.

Method: pyDESeq2 ratios, v10 simple-mean multi-gene aggregation
        (standard SOP, consistent with Tseng 2024 N&M).

Author: Yu-Yao Tseng  |  Date: 2026-04-18
"""

import os, sys, importlib.util
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.ticker import MultipleLocator
from scipy.integrate import odeint

plt.rcParams.update({
    'font.family':       'sans-serif',
    'font.sans-serif':   ['Calibri', 'Arial', 'DejaVu Sans'],
    'font.size':         7,
    'axes.labelsize':    7.5,
    'axes.titlesize':    8,
    'axes.linewidth':    0.5,
    'xtick.labelsize':   6.5,
    'ytick.labelsize':   6.5,
    'xtick.direction':   'out',
    'ytick.direction':   'out',
    'xtick.major.width': 0.4,
    'ytick.major.width': 0.4,
    'xtick.major.size':  2.5,
    'ytick.major.size':  2.5,
    'lines.linewidth':   1.0,
    'grid.alpha':        0.15,
    'grid.linewidth':    0.3,
    'savefig.dpi':       300,
})

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
REPO_ROOT = os.path.dirname(os.path.dirname(SCRIPT_DIR))
# Model files live in <repo>/data/models (repo standard; same as the
# 05_cross_cohort_validation scripts). Fall back to a local
# models_pydeseq2_mean/ folder next to this script if present.
_MODEL_SEARCH_DIRS = [
    os.path.join(REPO_ROOT, 'data', 'models'),
    os.path.join(SCRIPT_DIR, 'models_pydeseq2_mean'),
]


def _resolve_model(fname):
    """Return the full path to a model file, searched by basename."""
    base = os.path.basename(fname)
    for d in _MODEL_SEARCH_DIRS:
        p = os.path.join(d, base)
        if os.path.exists(p):
            return p
    return None

# ============================================================
# CONFIG
# ============================================================
MODEL_FILES = {
    'Obese': 'models_pydeseq2_mean/ode_model_pydeseq2_Obese_vs_Normal_v10_mean.py',
    'NAFL':  'models_pydeseq2_mean/ode_model_pydeseq2_NAFL_vs_Normal_v10_mean.py',
    'NASH':  'models_pydeseq2_mean/ode_model_pydeseq2_NASH_vs_Normal_v10_mean.py',
}

T_END = 300
T_POINTS = 7501
COND_ORDER = ['Normal', 'Obese', 'NAFL', 'NASH']

STYLE = {
    'Normal': {'color': '#2C3E50', 'ls': '-',  'lw': 1.0},
    'Obese':  {'color': '#3498DB', 'ls': '--', 'lw': 0.9},
    'NAFL':   {'color': '#E67E22', 'ls': '-.', 'lw': 1.0},
    'NASH':   {'color': '#C0392B', 'ls': '-',  'lw': 1.3},
}

# 2×3 panel layout
# ode_var = ODE variable name (without '_active' suffix)
# display = label shown in panel title
PANELS = [
    # Row 1 — upstream molecular nodes
    {'label': '(a)', 'display': 'CASP7',       'ode_var': 'CASP7',
     'ylabel': 'Relative level (a.u.)'},
    {'label': '(b)', 'display': r'TNF-$\alpha$', 'ode_var': 'TNFa',
     'ylabel': None},
    {'label': '(c)', 'display': 'IL-8',        'ode_var': 'IL_8',
     'ylabel': None},
    # Row 2 — downstream pathway outputs
    {'label': '(d)', 'display': 'P_Hepatocyte injury',
     'ode_var': 'P_Hepatocyte_injury',
     'ylabel': r'P$_{\mathrm{active}}$ (a.u.)'},
    {'label': '(e)', 'display': 'P_Cell death',
     'ode_var': 'P_Cell_death',
     'ylabel': None},
    {'label': '(f)', 'display': 'P_Inflammation',
     'ode_var': 'P_Inflammation',
     'ylabel': None},
]


# ============================================================
# HELPERS
# ============================================================
def load_model(filepath):
    name = os.path.splitext(os.path.basename(filepath))[0]
    spec = importlib.util.spec_from_file_location(name, filepath)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def build_normal_y0(SV):
    n = len(SV)
    y0 = np.zeros(n)
    for i, s in enumerate(SV):
        if s.startswith('P_') and s.endswith('_inactive'):
            y0[i] = 100.0
        elif s.endswith('_inactive'):
            y0[i] = 1.0
    return y0


# ============================================================
# MAIN
# ============================================================
def main():
    print("Loading models...")
    models = {}
    for label, fname in MODEL_FILES.items():
        fpath = _resolve_model(fname)
        if fpath is None:
            print(f"ERROR: model '{os.path.basename(fname)}' not found in {_MODEL_SEARCH_DIRS}"); sys.exit(1)
        models[label] = load_model(fpath)
        print(f"  {label}: {fname}")

    ref = models['NASH']
    SV = ref.STATE_VARS

    # Normal y0
    y0_normal = build_normal_y0(SV)

    # Build active-index map
    active_map = {s[:-7]: i for i, s in enumerate(SV)
                  if s.endswith('_active')}

    # Simulate
    print("\nRunning simulations...")
    t = np.linspace(0, T_END, T_POINTS)
    results = {}
    results['Normal'] = odeint(ref.ode_system, y0_normal, t,
                               args=(None,), mxstep=10000)
    for label in ['Obese', 'NAFL', 'NASH']:
        mod = models[label]
        results[label] = odeint(mod.ode_system, mod.Y0, t,
                                args=(None,), mxstep=10000)
    print("  Done.")

    # Summary — AUC + steady state
    last_n = max(1, int(T_POINTS * 0.1))
    summary_rows = []
    print(f"\n{'Panel':<5} {'Node':<28} {'NS/N AUC':>10} {'NS/N SS':>10}")
    print("-" * 65)
    for pdef in PANELS:
        ode_var = pdef['ode_var']
        pi = active_map.get(ode_var)
        if pi is None:
            print(f"{pdef['label']:<5} {ode_var:<28}  NOT FOUND")
            continue
        aucs = {c: float(np.trapezoid(results[c][:, pi], t))
                for c in COND_ORDER}
        sss = {c: float(np.mean(results[c][-last_n:, pi]))
               for c in COND_ORDER}
        ns_auc = aucs['NASH'] / aucs['Normal'] if aucs['Normal'] > 0 else 0
        ns_ss = sss['NASH'] / sss['Normal'] if sss['Normal'] > 0 else 0
        print(f"{pdef['label']:<5} {ode_var:<28} {ns_auc:>10.2f} {ns_ss:>10.2f}")
        summary_rows.append({
            'Panel': pdef['label'],
            'Node': ode_var,
            'Display': pdef['display'],
            'Normal_AUC': round(aucs['Normal'], 2),
            'Obese_AUC':  round(aucs['Obese'], 2),
            'NAFL_AUC':   round(aucs['NAFL'], 2),
            'NASH_AUC':   round(aucs['NASH'], 2),
            'AUC_NS/N':   round(ns_auc, 3),
            'Normal_SS':  round(sss['Normal'], 3),
            'NASH_SS':    round(sss['NASH'], 3),
            'SS_NS/N':    round(ns_ss, 3),
        })

    # Export Excel
    try:
        import pandas as pd
        xlsx_path = os.path.join(SCRIPT_DIR, 'Fig5_summary.xlsx')
        with pd.ExcelWriter(xlsx_path, engine='openpyxl') as writer:
            pd.DataFrame(summary_rows).to_excel(writer,
                sheet_name='Fig5_Summary', index=False)
            meta = pd.DataFrame([
                ['Script',       os.path.basename(__file__)],
                ['Date',         '2026-04-18'],
                ['Method',       'pyDESeq2 v10 Mean multi-gene aggregation'],
                ['Evaluation',   'AUC (trapezoidal integration) and steady state (last 10% mean)'],
                ['Simulation',   f't=0..{T_END}h, {T_POINTS} points'],
                ['Layout',       '2×3 (Row1: CASP7/TNF-α/CXCL8 ; Row2: P_Injury/P_Cell_death/P_Inflammation)'],
                ['Topology',     'CASP7→P_Hepatocyte_injury ; TNF→P_Cell_death ; CXCL8→P_Inflammation'],
            ], columns=['Parameter', 'Value'])
            meta.to_excel(writer, sheet_name='Metadata', index=False)
        print(f"\n  Excel saved: {xlsx_path}")
    except ImportError:
        pass

    # Plot 2×3
    print("\nGenerating figure...")
    fig, axes = plt.subplots(2, 3, figsize=(6.5, 4.3))
    fig.subplots_adjust(top=0.94, bottom=0.14, left=0.085, right=0.985,
                        hspace=0.55, wspace=0.35)

    for idx, pdef in enumerate(PANELS):
        row, col = idx // 3, idx % 3
        ax = axes[row][col]
        ode_var = pdef['ode_var']
        pi = active_map.get(ode_var)

        if pi is None:
            ax.text(0.5, 0.5, f"{ode_var}\nnot found",
                    transform=ax.transAxes, ha='center', va='center',
                    fontsize=8, color='red')
            continue

        for c in COND_ORDER:
            ax.plot(t, results[c][:, pi],
                    color=STYLE[c]['color'],
                    linestyle=STYLE[c]['ls'],
                    linewidth=STYLE[c]['lw'])

        ttl = ax.set_title(f"{pdef['label']} {pdef['display']}",
                           fontsize=8, fontweight='bold',
                           loc='left', pad=4)
        ttl.set_position((-0.05, 1.0))

        ax.set_xlabel('Time (h)')
        if pdef['ylabel'] is not None:
            ax.set_ylabel(pdef['ylabel'])
        ax.set_xlim(0, T_END)
        ax.xaxis.set_major_locator(MultipleLocator(75))
        ax.grid(True, linestyle=':', color='0.85', linewidth=0.25)
        for s in ax.spines.values():
            s.set_linewidth(0.4)

    # Single bottom legend
    leg_handles = [Line2D([0], [0],
                          color=STYLE[c]['color'],
                          ls=STYLE[c]['ls'],
                          lw=1.3, label=c)
                   for c in COND_ORDER]
    fig.legend(handles=leg_handles, loc='lower center',
               ncol=4, frameon=False, fontsize=8.5,
               bbox_to_anchor=(0.5, 0.01),
               handlelength=2.2, columnspacing=2.0)

    for ext in ['png', 'pdf']:
        out = os.path.join(SCRIPT_DIR, f'Fig5_dynamics_2x3.{ext}')
        fig.savefig(out, dpi=300 if ext == 'png' else None,
                    bbox_inches='tight', pad_inches=0.10)
        print(f"  Saved: {out}")
    plt.close(fig)
    print("\nDone!")


if __name__ == '__main__':
    main()
