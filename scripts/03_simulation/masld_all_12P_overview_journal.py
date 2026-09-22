#!/usr/bin/env python3
"""
All 12 P_* Pathway Outputs — 3×4 overview (v10)
================================================
Quick diagnostic: view all pathway outputs across 4 conditions.
Usage: python masld_all_12P_overview.py
"""

import os, sys, importlib.util
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from scipy.integrate import odeint

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))

MODEL_FILES = {
    'Obese': 'models_pydeseq2_mean/ode_model_pydeseq2_Obese_vs_Normal_v10_mean.py',
    'NAFL':  'models_pydeseq2_mean/ode_model_pydeseq2_NAFL_vs_Normal_v10_mean.py',
    'NASH':  'models_pydeseq2_mean/ode_model_pydeseq2_NASH_vs_Normal_v10_mean.py',
}

T_END = 200
T_POINTS = 5001
COND_ORDER = ['Normal', 'Obese', 'NAFL', 'NASH']

STYLE = {
    'Normal': {'color': '#2C3E50', 'ls': '-',  'lw': 1.0},
    'Obese':  {'color': '#3498DB', 'ls': '--', 'lw': 0.9},
    'NAFL':   {'color': '#E67E22', 'ls': '-.', 'lw': 1.0},
    'NASH':   {'color': '#C0392B', 'ls': '-',  'lw': 1.3},
}


def load_model(filepath):
    name = os.path.splitext(os.path.basename(filepath))[0]
    spec = importlib.util.spec_from_file_location(name, filepath)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def main():
    plt.rcParams.update({
        'font.family': 'sans-serif',
        'font.sans-serif': ['Calibri', 'Arial', 'DejaVu Sans'],
        'font.size': 7, 'axes.linewidth': 0.5,
    })

    # Load models
    print("Loading models...")
    models = {}
    for label, fname in MODEL_FILES.items():
        fpath = os.path.join(SCRIPT_DIR, fname)
        if not os.path.exists(fpath):
            print(f"  ERROR: {fpath} not found"); sys.exit(1)
        models[label] = load_model(fpath)
        print(f"  {label}: {fname}")

    ref = models['NASH']
    SV = ref.STATE_VARS
    n = len(SV)

    # Normal y0
    y0_normal = np.zeros(n)
    for i, s in enumerate(SV):
        if s.startswith('P_') and s.endswith('_inactive'):
            y0_normal[i] = 100.0
        elif s.endswith('_inactive'):
            y0_normal[i] = 1.0

    # Simulate
    print("Running simulations...")
    t = np.linspace(0, T_END, T_POINTS)
    results = {}
    results['Normal'] = odeint(ref.ode_system, y0_normal, t, args=(None,), mxstep=10000)
    for label in ['Obese', 'NAFL', 'NASH']:
        mod = models[label]
        results[label] = odeint(mod.ode_system, mod.Y0, t, args=(None,), mxstep=10000)
    print("  Done.")

    # P_* panels
    p_panels = [(i, s) for i, s in enumerate(SV)
                if s.startswith('P_') and s.endswith('_active')]
    last_n = max(1, int(T_POINTS * 0.1))

    # Steady-state summary
    print(f"\n{'#':<3} {'P_* Output':<30} {'Normal':>8} {'Obese':>8} "
          f"{'NAFL':>8} {'NASH':>8} | {'Ob/N':>6} {'NF/N':>6} {'NS/N':>6}")
    print("-" * 105)
    for idx, (pi, pname) in enumerate(p_panels):
        vals = {c: float(np.mean(results[c][-last_n:, pi])) for c in COND_ORDER}
        nv = vals['Normal']
        ratios = {c: vals[c]/nv if nv > 1e-10 else 0 for c in ['Obese','NAFL','NASH']}
        short = pname.replace('_active','').replace('P_','')
        print(f"{idx+1:<3} {short:<30} {vals['Normal']:8.2f} {vals['Obese']:8.2f} "
              f"{vals['NAFL']:8.2f} {vals['NASH']:8.2f} | "
              f"{ratios['Obese']:6.2f} {ratios['NAFL']:6.2f} {ratios['NASH']:6.2f}")

    # Plot 3×4
    fig, axes = plt.subplots(3, 4, figsize=(7.205, 6.1))   # npj SBA 雙欄寬度 183mm，加高留空間給 legend
    fig.subplots_adjust(hspace=0.65, wspace=0.35,
                        left=0.05, right=0.93, top=0.93, bottom=0.10)

    for idx, (pi, pname) in enumerate(p_panels):
        row, col = idx // 4, idx % 4
        ax = axes[row][col]
        short = pname.replace('_active','').replace('P_','').replace('_', ' ')

        for c in COND_ORDER:
            ax.plot(t, results[c][:, pi],
                    color=STYLE[c]['color'], linestyle=STYLE[c]['ls'],
                    linewidth=STYLE[c]['lw'])

        nv = float(np.mean(results['Normal'][-last_n:, pi]))
        nsv = float(np.mean(results['NASH'][-last_n:, pi]))
        ratio = nsv/nv if nv > 1e-10 else 0

        ax.set_title(f"P_{short}\nNS/N={ratio:.2f}", fontsize=8, fontweight='bold')
        ax.set_xlabel('Time (h)', fontsize=6.5)
        ax.set_xlim(0, T_END)
        ax.grid(True, alpha=0.15)
        ax.tick_params(labelsize=6)

        # Annotate steady-state values（固定用 axes 座標堆疊，跟資料數值
        # 接近程度無關，保證 4 個標籤永遠等距、不會重疊）
        vals = {c: float(np.mean(results[c][-last_n:, pi])) for c in COND_ORDER}
        order = sorted(COND_ORDER, key=lambda c: vals[c])  # 由小到大排序
        n_c = len(order)
        for k, c in enumerate(order):
            y_frac = 0.5 + (k - (n_c - 1) / 2) * 0.13  # 在 axes 中央上下對稱排開
            ax.annotate(f'{vals[c]:.1f}', xy=(1.0, y_frac), xycoords='axes fraction',
                        fontsize=4.5, color=STYLE[c]['color'], va='center', ha='left',
                        xytext=(3, 0), textcoords='offset points')

    # Legend
    leg_handles = [Line2D([0],[0], color=STYLE[c]['color'], ls=STYLE[c]['ls'],
                          lw=1.2, label=c) for c in COND_ORDER]
    fig.legend(handles=leg_handles, loc='lower center', ncol=4,
               frameon=False, fontsize=9, bbox_to_anchor=(0.5, -0.01))

    out = os.path.join(SCRIPT_DIR, 'All_12_P_outputs_journal.png')
    fig.savefig(out, dpi=300, bbox_inches='tight', facecolor='white')
    print(f"\nSaved: {out}")
    plt.close(fig)

    # Also save PDF
    out_pdf = os.path.join(SCRIPT_DIR, 'All_12_P_outputs_journal.pdf')
    fig2, axes2 = plt.subplots(3, 4, figsize=(7.205, 6.1))   # npj SBA 雙欄寬度 183mm
    fig2.subplots_adjust(hspace=0.65, wspace=0.35,
                         left=0.05, right=0.93, top=0.93, bottom=0.10)
    for idx, (pi, pname) in enumerate(p_panels):
        row, col = idx // 4, idx % 4
        ax = axes2[row][col]
        short = pname.replace('_active','').replace('P_','').replace('_', ' ')
        for c in COND_ORDER:
            ax.plot(t, results[c][:, pi],
                    color=STYLE[c]['color'], linestyle=STYLE[c]['ls'],
                    linewidth=STYLE[c]['lw'])
        nv = float(np.mean(results['Normal'][-last_n:, pi]))
        nsv = float(np.mean(results['NASH'][-last_n:, pi]))
        ratio = nsv/nv if nv > 1e-10 else 0
        ax.set_title(f"P_{short}\nNS/N={ratio:.2f}", fontsize=8, fontweight='bold')
        ax.set_xlabel('Time (h)', fontsize=6.5)
        ax.set_xlim(0, T_END)
        ax.grid(True, alpha=0.15)
        ax.tick_params(labelsize=6)
        vals = {c: float(np.mean(results[c][-last_n:, pi])) for c in COND_ORDER}
        order = sorted(COND_ORDER, key=lambda c: vals[c])
        n_c = len(order)
        for k, c in enumerate(order):
            y_frac = 0.5 + (k - (n_c - 1) / 2) * 0.13
            ax.annotate(f'{vals[c]:.1f}', xy=(1.0, y_frac), xycoords='axes fraction',
                        fontsize=4.5, color=STYLE[c]['color'], va='center', ha='left',
                        xytext=(3, 0), textcoords='offset points')
    fig2.legend(handles=leg_handles, loc='lower center', ncol=4,
                frameon=False, fontsize=9, bbox_to_anchor=(0.5, -0.01))
    fig2.savefig(out_pdf, bbox_inches='tight', facecolor='white')
    print(f"Saved: {out_pdf}")
    plt.close(fig2)
    print("\nDone!")


if __name__ == '__main__':
    main()
