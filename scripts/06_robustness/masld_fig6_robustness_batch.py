#!/usr/bin/env python3
"""
MASLD Fig 6 (v10 Mean, 2×3) — Silymarin intervention in NASH vs NAFL
======================================================================
In silico simulation of silymarin, a multi-target nutraceutical with strong
clinical evidence in NAFLD (Loguercio 2012; Wah-Kheong 2017; Chan 2017).

Silymarin is modeled as a simultaneous multi-node inhibitor of eight direct
molecular targets supported by published evidence:
  - CASP3, CASP7, CASP8  (caspase cascade; Amini 2011; Polyak 2007)
  - CYP2E1               (oxidative stress hub; Zhu 2014; Loguercio 2012)
  - IL-8 (CXCL8)         (Morishima 2010; Polyak 2007)
  - TNF-a                (Kang 2003; Manna 1999)
  - NF-kB                (Polyak 2007; Trappoliere 2009)
  - TGF-b1               (Trappoliere 2009)

These nodes regulate, directly or indirectly, the three validated pathway
outputs identified in Fig 5:
  - P_Hepatocyte_injury  (via caspases and oxidative stress)
  - P_Cell_death         (via TNF/FasL axis)
  - P_Inflammation       (via TGF-b1 and IL-8)

Intervention model (Tseng 2024, N&M 21:65, Eq. 4-5):
  activation_rate_intervention = activation_rate * ki
  ki = 1.0 (no drug), 0.7 (low), 0.5 (med), 0.3 (high)

Panels (2×3 layout):
  Row 1 (NASH intervention):
    (a) P_Hepatocyte_injury
    (b) P_Cell_death
    (c) P_Inflammation
  Row 2 (NAFL intervention — early-stage):
    (d) P_Hepatocyte_injury
    (e) P_Cell_death
    (f) P_Inflammation

Narrative: Compares silymarin efficacy in advanced NASH vs early NAFL.
Early-stage intervention is expected to be more effective, supporting
nutraceutical use as preventive therapy in early MASLD where no FDA-
approved pharmacotherapy is currently available.

Author: Yu-Yao Tseng  |  Date: 2026-04-18
"""

import os, sys, importlib.util, types
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
# ------------------------------------------------------------------
# ROBUSTNESS SWITCH: set PARAM_TAG to a variant suffix to load the
# parameter-tagged models produced by generate_robustness_models.py,
# e.g. PARAM_TAG = ''.  '' = original primary model.
# Output figure/Excel filenames inherit the same tag.
PARAM_TAG = ''   # overridden by batch loop below

# ==================================================================
# BATCH ROBUSTNESS: generate one figure set per parameter variant.
# Set PARAM_TAGS to the variants you want (must match files produced
# by generate_robustness_models.py). '' = primary (n=2, ksp=2) model.
# Each run writes files whose names encode the parameter combination,
# e.g. Fig6_silymarin_2x3_n1.0_ksp2.0.png .
# ==================================================================
PARAM_TAGS = [
    '',              # primary  (n=2.0, ksp=2.0)
    '_n1.0_ksp2.0',
    '_n4.0_ksp2.0',
    '_n2.0_ksp1.0',
    '_n2.0_ksp4.0',
]

_MODEL_SEARCH_DIRS = [
    os.path.join(REPO_ROOT, 'data', 'models', 'robustness_models'),
    os.path.join(SCRIPT_DIR, 'robustness_models'),
    os.path.join(REPO_ROOT, 'data', 'models'),
    os.path.join(SCRIPT_DIR, 'models_pydeseq2_mean'),
]


def _resolve_model(fname):
    """Return the full path to a model file, searched by basename."""
    base = os.path.basename(fname)
    if PARAM_TAG:
        base = base[:-3] + PARAM_TAG + '.py'
    for d in _MODEL_SEARCH_DIRS:
        p = os.path.join(d, base)
        if os.path.exists(p):
            return p
    return None

# ============================================================
# CONFIG
# ============================================================
BASE_MODELS = {
    'NASH': 'models_pydeseq2_mean/ode_model_pydeseq2_NASH_vs_Normal_v10_mean.py',
    'NAFL': 'models_pydeseq2_mean/ode_model_pydeseq2_NAFL_vs_Normal_v10_mean.py',
}

T_END    = 300
T_POINTS = 7501

SILYMARIN_DOSES = [
    {'ki': 1.0, 'label': 'Disease (no drug)', 'color': '#C0392B', 'ls': '-',  'lw': 1.5},
    {'ki': 0.7, 'label': 'Silymarin low',     'color': '#F39C12', 'ls': '--', 'lw': 1.1},
    {'ki': 0.5, 'label': 'Silymarin medium',  'color': '#16A085', 'ls': '-.', 'lw': 1.1},
    {'ki': 0.3, 'label': 'Silymarin high',    'color': '#2980B9', 'ls': '-',  'lw': 1.3},
]
NORMAL_STYLE = {'color': '#2C3E50', 'ls': ':', 'lw': 1.0, 'label': 'Normal (reference)'}

SILYMARIN_TARGETS = [
    ('CASP3',   18, 19),
    ('CASP7',   20, 21),
    ('CASP8',   22, 23),
    ('CYP2E1',  26, 27),
    ('IL_8',    58, 59),
    ('TNFa',   104, 105),
    ('NF_kB',   80, 81),
    ('TGF_b1', 100, 101),
]

P_OUTPUTS = [
    {'node': 'P_Hepatocyte_injury', 'title': 'P_Hepatocyte injury'},
    {'node': 'P_Cell_death',        'title': 'P_Cell death'},
    {'node': 'P_Inflammation',      'title': 'P_Inflammation'},
]
PANEL_LABELS = [['(a)', '(b)', '(c)'],
                ['(d)', '(e)', '(f)']]


# ============================================================
# HELPERS
# ============================================================
def build_silymarin_ode(base_source, ki):
    lines = base_source.split('\n')
    for _, i_idx, a_idx in SILYMARIN_TARGETS:
        for idx in (i_idx, a_idx):
            for li, line in enumerate(lines):
                prefix = f'dydt[{idx}]'
                if line.lstrip().startswith(prefix):
                    eq_pos = line.find('=')
                    rhs = line[eq_pos + 1:].strip()
                    if rhs.startswith('-('):
                        sign, body_start = '-', 1
                    elif rhs.startswith('('):
                        sign, body_start = '', 0
                    else:
                        continue
                    depth = 0
                    end_pos = None
                    for p, ch in enumerate(rhs[body_start:], start=body_start):
                        if ch == '(':
                            depth += 1
                        elif ch == ')':
                            depth -= 1
                            if depth == 0:
                                end_pos = p
                                break
                    if end_pos is None:
                        continue
                    act_block = rhs[body_start:end_pos + 1]
                    remainder = rhs[end_pos + 1:]
                    new_rhs = f"{sign}(({act_block}) * {ki:.4f}){remainder}"
                    indent = line[: len(line) - len(line.lstrip())]
                    lines[li] = f"{indent}dydt[{idx}] = {new_rhs}"
                    break
    return '\n'.join(lines)


def load_module_from_source(source, name):
    mod = types.ModuleType(name)
    mod.__file__ = name
    exec(compile(source, name, 'exec'), mod.__dict__)
    return mod


def load_base_model(filepath):
    spec = importlib.util.spec_from_file_location('base', filepath)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


# ============================================================
# MAIN
# ============================================================
def main():
    print(f"PARAM_TAG = {PARAM_TAG!r}  (empty = primary model)")
    base_mods, base_sources = {}, {}
    for cond, relpath in BASE_MODELS.items():
        fpath = _resolve_model(relpath)
        if fpath is None:
            print(f"ERROR: model '{os.path.basename(relpath)}' not found in {_MODEL_SEARCH_DIRS}"); sys.exit(1)
        base_mods[cond] = load_base_model(fpath)
        with open(fpath, 'r', encoding='utf-8') as f:
            base_sources[cond] = f.read()
        print(f"Loaded {cond}: {os.path.basename(fpath)}")

    SV = base_mods['NASH'].STATE_VARS
    y0_normal = np.zeros(len(SV))
    for i, s in enumerate(SV):
        if s.startswith('P_') and s.endswith('_inactive'):
            y0_normal[i] = 100.0
        elif s.endswith('_inactive'):
            y0_normal[i] = 1.0

    active_map = {s[:-7]: i for i, s in enumerate(SV) if s.endswith('_active')}
    t = np.linspace(0, T_END, T_POINTS)

    print(f"\nSimulating over t=0..{T_END}h ({T_POINTS} points)...")
    results = {}
    results['Normal'] = odeint(base_mods['NASH'].ode_system, y0_normal, t,
                               args=(None,), mxstep=10000)
    print("  Normal: done")

    for cond in ['NASH', 'NAFL']:
        base_mod = base_mods[cond]
        base_src = base_sources[cond]
        Y0 = base_mod.Y0
        for d in SILYMARIN_DOSES:
            ki = d['ki']
            key = f"{cond}_ki{ki}"
            if ki == 1.0:
                mod = base_mod
            else:
                mod_src = build_silymarin_ode(base_src, ki)
                mod = load_module_from_source(mod_src, f'{cond}_silymarin_ki{ki}')
            results[key] = odeint(mod.ode_system, Y0, t, args=(None,), mxstep=10000)
            print(f"  {key}: done")

    # Summary
    last_n = max(1, int(T_POINTS * 0.1))
    normal_auc = {p['node']: float(np.trapezoid(results['Normal'][:, active_map[p['node']]], t))
                  for p in P_OUTPUTS}
    summary_rows = []
    print(f"\n{'Cond':<5} {'P_*':<24} {'ki':>4} {'AUC':>10} {'Steady':>8} {'AUC/Norm':>9}")
    print("-" * 70)
    for cond in ['NASH', 'NAFL']:
        for pdef in P_OUTPUTS:
            pi = active_map[pdef['node']]
            for d in SILYMARIN_DOSES:
                ki = d['ki']
                y = results[f"{cond}_ki{ki}"][:, pi]
                auc = float(np.trapezoid(y, t))
                ss  = float(np.mean(y[-last_n:]))
                auc_n = auc / normal_auc[pdef['node']]
                print(f"{cond:<5} {pdef['node']:<24} {ki:>4.1f} "
                      f"{auc:>10.1f} {ss:>8.2f} {auc_n:>9.2f}")
                summary_rows.append({
                    'Condition': cond, 'P_Output': pdef['node'],
                    'Silymarin_ki': ki, 'AUC': round(auc, 2),
                    'Steady_state': round(ss, 3),
                    'AUC_vs_Normal': round(auc_n, 3),
                })
            print()

    # Excel
    try:
        import pandas as pd
        xlsx = os.path.join(SCRIPT_DIR, f'Fig6_summary{PARAM_TAG}.xlsx')
        with pd.ExcelWriter(xlsx, engine='openpyxl') as w:
            pd.DataFrame(summary_rows).to_excel(w, sheet_name='Fig6_Summary', index=False)
            meta = pd.DataFrame([
                ['Script',       os.path.basename(__file__)],
                ['Date',         '2026-04-18'],
                ['Intervention', 'Silymarin (multi-target nutraceutical, 8 direct targets)'],
                ['Targets',      'CASP3, CASP7, CASP8, CYP2E1, IL-8, TNFa, NF-kB, TGF-b1'],
                ['Method',       'Tseng 2024 N&M ki multiplier (Eq. 4-5)'],
                ['Dose levels',  f"ki = {[d['ki'] for d in SILYMARIN_DOSES]}"],
                ['Conditions',   'NASH (advanced) + NAFL (early-stage)'],
                ['Simulation',   f"t=0..{T_END}h, {T_POINTS} points"],
                ['Clinical ref', 'Wah-Kheong 2017 CGH 15:1940; Loguercio 2012 FRBM 52:1658'],
            ], columns=['Parameter', 'Value'])
            meta.to_excel(w, sheet_name='Metadata', index=False)
        print(f"  Excel saved: {xlsx}")
    except ImportError:
        pass

    # Plot 2×3
    print("\nGenerating figure...")
    fig, axes = plt.subplots(2, 3, figsize=(6.8, 4.8))
    fig.subplots_adjust(top=0.93, bottom=0.13, left=0.09, right=0.985,
                        hspace=0.55, wspace=0.33)

    for row, cond in enumerate(['NASH', 'NAFL']):
        for col, pdef in enumerate(P_OUTPUTS):
            ax = axes[row][col]
            pi = active_map[pdef['node']]
            label = PANEL_LABELS[row][col]

            ax.plot(t, results['Normal'][:, pi],
                    color=NORMAL_STYLE['color'], linestyle=NORMAL_STYLE['ls'],
                    linewidth=NORMAL_STYLE['lw'])
            for d in SILYMARIN_DOSES:
                ax.plot(t, results[f"{cond}_ki{d['ki']}"][:, pi],
                        color=d['color'], linestyle=d['ls'],
                        linewidth=d['lw'])

            title = f"{label} {cond}: {pdef['title']}"
            ttl = ax.set_title(title, fontsize=8, fontweight='bold',
                               loc='left', pad=4)
            ttl.set_position((-0.05, 1.0))

            if row == 1:
                ax.set_xlabel('Time (h)')
            if col == 0:
                ax.set_ylabel(r'P$_{\mathrm{active}}$ (a.u.)')
            ax.set_xlim(0, T_END)
            ax.xaxis.set_major_locator(MultipleLocator(75))
            ax.grid(True, linestyle=':', color='0.85', linewidth=0.25)
            for s in ax.spines.values():
                s.set_linewidth(0.4)

    leg_handles = [
        Line2D([0], [0], color=NORMAL_STYLE['color'],
               ls=NORMAL_STYLE['ls'], lw=1.3,
               label=NORMAL_STYLE['label']),
    ] + [
        Line2D([0], [0], color=d['color'], ls=d['ls'], lw=1.3,
               label=d['label'])
        for d in SILYMARIN_DOSES
    ]
    fig.legend(handles=leg_handles, loc='lower center',
               ncol=5, frameon=False, fontsize=7.5,
               bbox_to_anchor=(0.5, 0.005),
               handlelength=2.2, columnspacing=1.2)

    for ext in ['png', 'pdf']:
        out = os.path.join(SCRIPT_DIR, f'Fig6_silymarin_2x3{PARAM_TAG}.{ext}')
        fig.savefig(out, dpi=300 if ext == 'png' else None,
                    bbox_inches='tight', pad_inches=0.10)
        print(f"  Saved: {out}")
    plt.close(fig)
    print("\nDone!")


if __name__ == '__main__':
    _tags = PARAM_TAGS if 'PARAM_TAGS' in globals() and PARAM_TAGS else ['']
    for _tag in _tags:
        PARAM_TAG = _tag                      # rebind module global
        globals()['PARAM_TAG'] = _tag
        print("\n" + "=" * 60)
        print(f"RUN  PARAM_TAG = {_tag!r}")
        print("=" * 60)
        try:
            main()
        except SystemExit as e:
            print(f"  skipped {_tag!r}: {e}")
