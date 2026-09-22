#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
fig_ValidationCrossCohort_6datasets.py
================================================================
6-dataset cross-cohort validation figures for P_Cell_death,
P_Hepatocyte_injury, P_Inflammation.

Datasets (in row order):
  1. GSE126848 (RNA-seq, main model)         — Suppli 2019
  2. GSE48452  (Affymetrix, cohort 1)         — Ahrens 2013
  3. GSE89632  (Illumina, cohort 2)           — Arendt 2015
  4. GSE130970 (RNA-seq, cohort 3)            — Hoang 2019
  5. GSE162694 (RNA-seq, cohort 4)            — Pantano 2021
  6. GSE213621 (RNA-seq, cohort 5)            — Chen 2023

Figure A (6x3): Disease progression dynamics
  Rows = 6 datasets | Cols = 3 P_* outputs
  Lines = Normal / Obese (N/A if unavailable) / NAFL / NASH
  Each row label includes sample n (Normal/NAFL/NASH, +Obese if available)

Figure B (12x3): Silymarin intervention
  Rows = NASH/NAFL x 6 datasets = 12 rows | Cols = 3 P_* outputs
  Lines = Normal(ref) / ki=1.0(disease) / 0.7 / 0.5 / 0.3

FONT SIZES: enlarged throughout for A4 print legibility (titles 11-13pt,
row/col labels 10pt, ticks 9pt, legend 10-11pt) compared to the original
3-dataset version (7-8.5pt).

DIRECTORY STRUCTURE EXPECTED (fully self-contained within this repository):
  <REPO_ROOT>/
    data/models/
      ode_model_pydeseq2_{NASH,NAFL,Obese}_vs_Normal_v10_mean.py
    data/ratios/
      GSE48452_node_initial_ratios.xlsx
      GSE89632_node_initial_ratios.xlsx
      GSE130970_node_initial_ratios.xlsx
      GSE162694_node_initial_ratios.xlsx
      GSE213621_node_initial_ratios.xlsx

No external paths need to be set; everything resolves relative to this
script's location within the repository.

Author: Yu-Yao Tseng  |  2026
================================================================
"""
import os, importlib.util, types
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.integrate import odeint

# ================================================================
# Fully self-contained: all data paths resolve relative to this repository
# ================================================================
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
REPO_ROOT = os.path.dirname(os.path.dirname(SCRIPT_DIR))
MODEL_DIR = os.path.join(REPO_ROOT, 'data', 'models')
RATIOS_DIR = os.path.join(REPO_ROOT, 'data', 'ratios')
# ================================================================

MODEL_NASH  = os.path.join(MODEL_DIR, 'ode_model_pydeseq2_NASH_vs_Normal_v10_mean.py')
MODEL_NAFL  = os.path.join(MODEL_DIR, 'ode_model_pydeseq2_NAFL_vs_Normal_v10_mean.py')
MODEL_OBESE = os.path.join(MODEL_DIR, 'ode_model_pydeseq2_Obese_vs_Normal_v10_mean.py')

RATIOS_48452  = os.path.join(RATIOS_DIR, 'GSE48452_node_initial_ratios.xlsx')
RATIOS_89632  = os.path.join(RATIOS_DIR, 'GSE89632_node_initial_ratios.xlsx')
RATIOS_130970 = os.path.join(RATIOS_DIR, 'GSE130970_node_initial_ratios.xlsx')
RATIOS_162694 = os.path.join(RATIOS_DIR, 'GSE162694_node_initial_ratios.xlsx')
RATIOS_213621 = os.path.join(RATIOS_DIR, 'GSE213621_node_initial_ratios.xlsx')

OUT_FIG_A = os.path.join(SCRIPT_DIR,
    'Fig_ValidationCrossCohort_DiseaseProgression_6datasets_3outputs_6x3.png')
OUT_FIG_B = os.path.join(SCRIPT_DIR,
    'Fig_ValidationCrossCohort_Silymarin_NAFL_NASH_6datasets_3outputs_6x6.png')
OUT_AUC_XLSX = os.path.join(SCRIPT_DIR,
    'Table_Silymarin_AUCreduction_NAFL_vs_NASH_6datasets.xlsx')

# ----------------------------------------------------------------
T_MAX, N_POINTS = 300.0, 3001
t = np.linspace(0, T_MAX, N_POINTS)
HALF_TARGET = 50.0
OUTPUTS = ['P_Cell_death', 'P_Hepatocyte_injury', 'P_Inflammation']
OUTPUT_LABELS = ['P_Cell death', 'P_Hepatocyte injury', 'P_Inflammation']

SILYMARIN_TARGETS = [
    ('CASP3',18,19),('CASP7',20,21),('CASP8',22,23),
    ('CYP2E1',26,27),('IL_8',58,59),('TNFa',104,105),
    ('NF_kB',80,81),('TGF_b1',100,101),
]

DOSE_STYLES = {
    1.0: {'label': 'Disease (no drug)', 'color': '#C0392B', 'ls': '-',  'lw': 1.8},
    0.7: {'label': 'Silymarin low',     'color': '#E67E22', 'ls': '--', 'lw': 1.3},
    0.5: {'label': 'Silymarin medium',  'color': '#16A085', 'ls': '-.', 'lw': 1.3},
    0.3: {'label': 'Silymarin high',    'color': '#2980B9', 'ls': '-',  'lw': 1.8},
}
NORMAL_STYLE = {'label': 'Normal (ref)', 'color': '#2C3E50', 'ls': ':', 'lw': 1.3}

STAGE_STYLES = {
    'Normal': {'color': '#2C3E50', 'ls': '-',  'lw': 1.3},
    'Obese':  {'color': '#3498DB', 'ls': '--', 'lw': 1.2},
    'NAFL':   {'color': '#E67E22', 'ls': '-.', 'lw': 1.3},
    'NASH':   {'color': '#C0392B', 'ls': '-',  'lw': 1.6},
}

# Font settings for A3 printing (substantially enlarged; title/axis labels/ticks/legend are all larger than the A4 version)
plt.rcParams.update({
    'font.family': 'sans-serif',
    'font.sans-serif': ['Arial', 'DejaVu Sans'],
    'font.size': 13,
    'axes.labelsize': 13,
    'axes.titlesize': 16,
    'axes.linewidth': 1.0,
    'xtick.labelsize': 11,
    'ytick.labelsize': 11,
    'lines.linewidth': 1.6,
    'savefig.dpi': 300,
})

# A3 paper size (inch): portrait 11.7x16.5, landscape 16.5x11.7
A3_PORTRAIT  = (11.7, 16.5)
A3_LANDSCAPE = (16.5, 11.7)

# ================================================================
# Sample sizes for the 6 datasets (all verified one by one against the
# official GEO series_matrix metadata, not estimated from literature citations):
#   GSE126848: paper Table 1 (Suppli 2019)
#   GSE48452 : measured from the series_matrix char_2 "group" field (Control/Healthy
#              obese/Steatosis/Nash = 14/27/14/18, consistent with Ahrens 2013)
#   GSE89632 : measured from the series_matrix char_2 "diagnosis" field (HC/SS/NASH
#              = 24/20/19, consistent with Arendt 2015)
#   GSE130970/162694/213621: sample sizes after grouping, as actually
#              produced by running compute_node_ratios_*.py in this project
# ================================================================
SAMPLE_N = {
    'GSE126848': {'Normal': 14, 'Obese': 12, 'NAFL': 15, 'NASH': 16},
    'GSE48452':  {'Normal': 14, 'Obese': 27, 'NAFL': 14, 'NASH': 18},
    'GSE89632':  {'Normal': 24, 'NAFL': 20, 'NASH': 19},
    'GSE130970': {'Normal': 4,  'NAFL': 10, 'NASH': 26},
    'GSE162694': {'Normal': 31, 'NAFL': 65, 'NASH': 47},
    'GSE213621': {'Normal': 69, 'NAFL': 97, 'NASH': 202},
}


def n_label(ds, stages):
    """Builds a sample-size label string like '(n=14/15/16)', in stage order."""
    d = SAMPLE_N.get(ds, {})
    parts = [str(d[s]) if s in d else 'NA' for s in stages]
    return '(n=' + '/'.join(parts) + ')'


# ================================================================
# HELPERS
# ================================================================
def load_model(path):
    spec = importlib.util.spec_from_file_location(
        os.path.splitext(os.path.basename(path))[0], path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def make_normal_y0(model):
    y0 = np.zeros(len(model.STATE_VARS))
    for i, sv in enumerate(model.STATE_VARS):
        if sv.endswith('_inactive'):
            node = sv[:-len('_inactive')]
            if node in model.PATHWAY_METADATA:
                y0[i] = model.PATHWAY_METADATA[node].get('initial_level', 100.)
            else:
                y0[i] = 1.0
    return y0


def make_y0_from_ratios(model, ratio_map):
    y0 = np.zeros(len(model.STATE_VARS))
    for i, sv in enumerate(model.STATE_VARS):
        if sv.endswith('_inactive'):
            node = sv[:-len('_inactive')]
            if node in ratio_map:
                y0[i] = ratio_map[node]
            elif node in model.PATHWAY_METADATA:
                y0[i] = model.PATHWAY_METADATA[node].get('initial_level', 100.)
            else:
                y0[i] = 1.0
    return y0


def load_ratios(xlsx_path, stage):
    if not os.path.exists(xlsx_path):
        return None
    xls = pd.ExcelFile(xlsx_path)
    sheet = f'{stage}_ratios'
    if sheet not in xls.sheet_names:
        return None
    df = pd.read_excel(xls, sheet)
    return dict(zip(df['node'], df['ratio']))


def build_sily_source(src, ki):
    lines = src.split('\n')
    for _, ii, ai in SILYMARIN_TARGETS:
        for idx in (ii, ai):
            for li, line in enumerate(lines):
                if line.lstrip().startswith(f'dydt[{idx}]'):
                    eq  = line.find('='); rhs = line[eq+1:].strip()
                    if   rhs.startswith('-('):  sign, bs = '-', 1
                    elif rhs.startswith('('):   sign, bs = '',  0
                    else: continue
                    depth = 0; ep = None
                    for p, ch in enumerate(rhs[bs:], start=bs):
                        if ch == '(':  depth += 1
                        elif ch == ')':
                            depth -= 1
                            if depth == 0: ep = p; break
                    if ep is None: continue
                    act = rhs[bs:ep+1]; rem = rhs[ep+1:]
                    ind = line[:len(line)-len(line.lstrip())]
                    lines[li] = f"{ind}dydt[{idx}] = {sign}(({act})*{ki:.4f}){rem}"
                    break
    return '\n'.join(lines)


def make_mod_from_source(src, name):
    mod = types.ModuleType(name); mod.__file__ = name
    exec(compile(src, name, 'exec'), mod.__dict__)
    return mod


def run_sim(model, y0):
    return odeint(model.ode_system, y0, t, args=(None,), mxstep=10000)


def annotate_na(ax, text='Obese: N/A'):
    ax.text(0.97, 0.97, text, transform=ax.transAxes,
            ha='right', va='top', fontsize=10, color='#3498DB',
            fontstyle='italic',
            bbox=dict(boxstyle='round,pad=0.3', fc='white', ec='#3498DB',
                      alpha=0.75, linewidth=0.7))


def t_half_val(traj):
    above = np.where(traj >= HALF_TARGET)[0]
    return float(t[above[0]]) if len(above) else np.inf


def auc_pct_reduction(traj_drug, traj_nodrug):
    a0 = np.trapezoid(traj_nodrug, t)
    a1 = np.trapezoid(traj_drug, t)
    return (a0 - a1) / a0 * 100 if a0 > 0 else np.nan


# ================================================================
# MAIN
# ================================================================
def main():
    print('[File check]')
    required = [(MODEL_NASH,'MODEL_NASH'), (MODEL_NAFL,'MODEL_NAFL'),
                (RATIOS_48452,'RATIOS_48452'), (RATIOS_89632,'RATIOS_89632'),
                (RATIOS_130970,'RATIOS_130970'), (RATIOS_162694,'RATIOS_162694'),
                (RATIOS_213621,'RATIOS_213621')]
    for f, lb in required + [(MODEL_OBESE,'MODEL_OBESE')]:
        ok = os.path.exists(f)
        print(f'  {"OK  " if ok else "MISS"} {lb}: {f}')
        if not ok and lb not in ('MODEL_OBESE',):
            raise FileNotFoundError(f'Required: {f}')

    print('\n[Loading main models]')
    m_nash  = load_model(MODEL_NASH)
    m_nafl  = load_model(MODEL_NAFL)
    m_obese = load_model(MODEL_OBESE) if os.path.exists(MODEL_OBESE) else None

    SV = m_nash.STATE_VARS
    p_idx = {pw: SV.index(pw+'_active') for pw in OUTPUTS}

    with open(MODEL_NASH, encoding='utf-8') as fh:
        base_src = fh.read()

    # ---------- SIMULATE: GSE126848 ----------
    print('[Simulating GSE126848]')
    y0_norm  = make_normal_y0(m_nash)
    y0_nafl  = make_y0_from_ratios(m_nash, {n: m_nafl.GENE_METADATA[n]['initial_ratio']
                                             for n in m_nafl.GENE_METADATA})
    y0_nash  = make_y0_from_ratios(m_nash, {n: m_nash.GENE_METADATA[n]['initial_ratio']
                                             for n in m_nash.GENE_METADATA})

    traj_126 = {'Normal': run_sim(m_nash, y0_norm),
                'NAFL':   run_sim(m_nafl, y0_nafl),
                'NASH':   run_sim(m_nash, y0_nash)}
    if m_obese:
        y0_obese = make_y0_from_ratios(m_nash,
                        {n: m_obese.GENE_METADATA[n]['initial_ratio']
                         for n in m_obese.GENE_METADATA})
        traj_126['Obese'] = run_sim(m_obese, y0_obese)
    else:
        traj_126['Obese'] = None

    # ---------- SIMULATE: 5 external validation cohorts (sharing the same Normal trajectory) ----------
    def sim_external(ratios_path, stages):
        traj = {'Normal': traj_126['Normal']}
        ratio_dict = {}
        for stage in stages:
            r = load_ratios(ratios_path, stage)
            if r:
                ratio_dict[stage] = r
                y0 = make_y0_from_ratios(m_nash, r)
                traj[stage] = run_sim(m_nash, y0)
            else:
                traj[stage] = None
        if 'Obese' not in traj:
            traj['Obese'] = None
        return traj, ratio_dict

    print('[Simulating GSE48452]')
    traj_48, ratios48 = sim_external(RATIOS_48452, ['Obese', 'NAFL', 'NASH'])
    print('[Simulating GSE89632]')
    traj_89, ratios89 = sim_external(RATIOS_89632, ['NAFL', 'NASH'])
    print('[Simulating GSE130970]')
    traj_130970, ratios130970 = sim_external(RATIOS_130970, ['NAFL', 'NASH'])
    print('[Simulating GSE162694]')
    traj_162694, ratios162694 = sim_external(RATIOS_162694, ['NAFL', 'NASH'])
    print('[Simulating GSE213621]')
    traj_213621, ratios213621 = sim_external(RATIOS_213621, ['NAFL', 'NASH'])

    # ============================================================
    # FIGURE A: Disease Progression 6×3
    # ============================================================
    print('\n[Figure A: Disease progression 6x3]')

    dataset_trajs = [
        ('GSE126848', 'GSE126848 (RNA-seq,\nmain model)', traj_126,
         ['Normal','Obese','NAFL','NASH']),
        ('GSE48452',  'GSE48452 (Affymetrix,\ncohort 1)',  traj_48,
         ['Normal','Obese','NAFL','NASH']),
        ('GSE89632',  'GSE89632 (Illumina,\ncohort 2)',    traj_89,
         ['Normal','NAFL','NASH']),
        ('GSE130970', 'GSE130970 (RNA-seq,\ncohort 3)',    traj_130970,
         ['Normal','NAFL','NASH']),
        ('GSE162694', 'GSE162694 (RNA-seq,\ncohort 4)',    traj_162694,
         ['Normal','NAFL','NASH']),
        ('GSE213621', 'GSE213621 (RNA-seq,\ncohort 5)',    traj_213621,
         ['Normal','NAFL','NASH']),
    ]

    fig_a, axes_a = plt.subplots(6, 3, figsize=A3_PORTRAIT, constrained_layout=True)

    for row, (ds_id, ds_label, traj_dict, stage_order) in enumerate(dataset_trajs):
        for col, (pw, pw_label) in enumerate(zip(OUTPUTS, OUTPUT_LABELS)):
            ax = axes_a[row][col]
            i  = p_idx[pw]

            for stage in ['Normal','Obese','NAFL','NASH']:
                tr = traj_dict.get(stage)
                if tr is None:
                    continue
                st = STAGE_STYLES[stage]
                ax.plot(t, tr[:, i], color=st['color'],
                        ls=st['ls'], lw=st['lw'], label=stage)

            if traj_dict.get('Obese') is None and 'Obese' not in stage_order:
                annotate_na(ax)

            if col == 0:
                ax.set_ylabel('P_active (a.u.)', fontsize=13)
                ax.text(-0.38, 0.5, f'{ds_label}\n{n_label(ds_id, stage_order)}',
                        transform=ax.transAxes, rotation=90,
                        va='center', ha='center', fontsize=12.5, fontweight='bold')
            if row == 5:
                ax.set_xlabel('Time (h)', fontsize=13)
            if row == 0:
                ax.set_title(pw_label, fontsize=16, fontweight='bold')

            ax.grid(alpha=0.2)
            ax.spines['top'].set_visible(False)
            ax.spines['right'].set_visible(False)

    legend_elements = [
        plt.Line2D([0],[0], color=STAGE_STYLES[s]['color'], ls=STAGE_STYLES[s]['ls'],
                   lw=2.4, label=s) for s in ['Normal','Obese','NAFL','NASH']
    ]
    fig_a.legend(handles=legend_elements, loc='lower center', ncol=4, fontsize=14,
                 bbox_to_anchor=(0.5, -0.01), framealpha=0.9)
    fig_a.suptitle(
        'Disease-stage P_* dynamics: cross-cohort comparison across 6 datasets (t = 300 h)',
        fontsize=18, fontweight='bold')

    fig_a.savefig(OUT_FIG_A, dpi=300, bbox_inches='tight')
    plt.close(fig_a)
    print(f'  Saved: {OUT_FIG_A}')

    # ============================================================
    # FIGURE B: Silymarin 6x6 (rows=6 datasets, cols=NASH x3 outputs + NAFL x3 outputs)
    # Uses 6x6 instead of 12x3, paired with A3 landscape printing, to avoid excessive blank space from an overly long vertical layout
    # ============================================================
    print('\n[Figure B: Silymarin 6x6 (A3 landscape)]')

    sily_mods = {1.0: m_nash}
    for ki in [0.7, 0.5, 0.3]:
        src = build_sily_source(base_src, ki)
        sily_mods[ki] = make_mod_from_source(src, f'sily_{ki}')
    print('  Silymarin modules compiled.')

    stage_y0s = {
        'GSE126848': {'NAFL': y0_nafl, 'NASH': y0_nash},
        'GSE48452':  {s: (make_y0_from_ratios(m_nash, ratios48[s]) if s in ratios48 else None)
                      for s in ['NAFL','NASH']},
        'GSE89632':  {s: (make_y0_from_ratios(m_nash, ratios89[s]) if s in ratios89 else None)
                      for s in ['NAFL','NASH']},
        'GSE130970': {s: (make_y0_from_ratios(m_nash, ratios130970[s]) if s in ratios130970 else None)
                      for s in ['NAFL','NASH']},
        'GSE162694': {s: (make_y0_from_ratios(m_nash, ratios162694[s]) if s in ratios162694 else None)
                      for s in ['NAFL','NASH']},
        'GSE213621': {s: (make_y0_from_ratios(m_nash, ratios213621[s]) if s in ratios213621 else None)
                      for s in ['NAFL','NASH']},
    }

    ds_short_labels = {
        'GSE126848': 'GSE126848 (RNA-seq)',
        'GSE48452':  'GSE48452 (Affymetrix)',
        'GSE89632':  'GSE89632 (Illumina)',
        'GSE130970': 'GSE130970 (RNA-seq)',
        'GSE162694': 'GSE162694 (RNA-seq)',
        'GSE213621': 'GSE213621 (RNA-seq)',
    }
    ds_order = ['GSE126848','GSE48452','GSE89632','GSE130970','GSE162694','GSE213621']

    # 6 cols: (stage, output) — the three NASH outputs on the left half, the three NAFL outputs on the right half
    col_order = [('NASH', pw) for pw in OUTPUTS] + [('NAFL', pw) for pw in OUTPUTS]
    col_labels = [f'NASH\n{lb}' for lb in OUTPUT_LABELS] + [f'NAFL\n{lb}' for lb in OUTPUT_LABELS]

    fig_b, axes_b = plt.subplots(6, 6, figsize=A3_LANDSCAPE, constrained_layout=True)

    y0_norm_sily = make_normal_y0(m_nash)
    norm_traj    = run_sim(m_nash, y0_norm_sily)

    print('\n  AUC reduction % at ki=0.3 (NAFL vs NASH):')
    print(f'  {"Dataset":<12} {"Output":<28} {"NAFL":>8} {"NASH":>8} {"ratio":>7}')
    print('  ' + '-'*65)

    auc_summary = []

    for row, ds in enumerate(ds_order):
        ds_short = ds_short_labels[ds]
        for col, (stage, pw) in enumerate(col_order):
            ax = axes_b[row][col]
            y0_stage = stage_y0s[ds][stage]
            i = p_idx[pw]

            if y0_stage is None:
                ax.text(0.5, 0.5, 'N/A', transform=ax.transAxes,
                        ha='center', va='center', fontsize=13, color='gray')
            else:
                ax.plot(t, norm_traj[:, i],
                        color=NORMAL_STYLE['color'], ls=NORMAL_STYLE['ls'],
                        lw=NORMAL_STYLE['lw'], label=NORMAL_STYLE['label'])
                for ki in [1.0, 0.7, 0.5, 0.3]:
                    ds_ki = DOSE_STYLES[ki]
                    drug_traj = run_sim(sily_mods[ki], y0_stage)
                    ax.plot(t, drug_traj[:, i],
                            color=ds_ki['color'], ls=ds_ki['ls'],
                            lw=ds_ki['lw'], label=ds_ki['label'])

            ax.grid(alpha=0.2)
            ax.spines['top'].set_visible(False)
            ax.spines['right'].set_visible(False)
            ax.tick_params(labelsize=9.5)

            n_str = SAMPLE_N.get(ds, {}).get(stage, 'NA')
            if col == 0:
                ax.set_ylabel('P_active (a.u.)', fontsize=11)
                ax.text(-0.42, 0.5, f'{ds_short}',
                        transform=ax.transAxes, rotation=90,
                        va='center', ha='center', fontsize=11.5, fontweight='bold')
            if row == 0:
                ax.set_title(col_labels[col], fontsize=12.5, fontweight='bold')
            # small n label in the upper-right corner of each panel
            ax.text(0.97, 0.06, f'n={n_str}', transform=ax.transAxes,
                    ha='right', va='bottom', fontsize=8, color='#555',
                    fontstyle='italic')
            if row == 5:
                ax.set_xlabel('Time (h)', fontsize=11)

    # AUC summary
    for ds in ds_order:
        for pw in OUTPUTS:
            i = p_idx[pw]
            vals = {}
            for stage in ['NAFL','NASH']:
                y0_s = stage_y0s[ds][stage]
                if y0_s is None:
                    vals[stage] = np.nan; continue
                tr_nodrug = run_sim(m_nash, y0_s)
                tr_drug   = run_sim(sily_mods[0.3], y0_s)
                vals[stage] = auc_pct_reduction(tr_drug[:,i], tr_nodrug[:,i])
            nafl_v = vals.get('NAFL', np.nan)
            nash_v = vals.get('NASH', np.nan)
            ratio  = nafl_v/nash_v if (np.isfinite(nafl_v) and np.isfinite(nash_v)
                                        and nash_v > 0) else np.nan
            auc_summary.append({
                'Dataset': ds, 'P_output': pw,
                'NAFL_AUC_red%': round(nafl_v,2) if np.isfinite(nafl_v) else np.nan,
                'NASH_AUC_red%': round(nash_v,2) if np.isfinite(nash_v) else np.nan,
                'NAFL/NASH_ratio': round(ratio,2) if np.isfinite(ratio) else np.nan,
            })
            nafl_s = f'{nafl_v:.1f}%' if np.isfinite(nafl_v) else 'N/A'
            nash_s = f'{nash_v:.1f}%' if np.isfinite(nash_v) else 'N/A'
            ratio_s = f'{ratio:.2f}x' if np.isfinite(ratio) else 'N/A'
            print(f'  {ds:<12} {pw:<28} {nafl_s:>8} {nash_s:>8} {ratio_s:>7}')

    auc_df = pd.DataFrame(auc_summary)
    auc_df.to_excel(OUT_AUC_XLSX, index=False)
    print(f'\n  AUC summary saved: {OUT_AUC_XLSX}')

    legend_b = [
        plt.Line2D([0],[0], color=NORMAL_STYLE['color'], ls=NORMAL_STYLE['ls'],
                   lw=2.0, label='Normal (reference)'),
    ] + [
        plt.Line2D([0],[0], color=DOSE_STYLES[ki]['color'],
                   ls=DOSE_STYLES[ki]['ls'], lw=2.0,
                   label=DOSE_STYLES[ki]['label'])
        for ki in [1.0, 0.7, 0.5, 0.3]
    ]
    fig_b.legend(handles=legend_b, loc='lower center', ncol=5, fontsize=13,
                 bbox_to_anchor=(0.5, -0.02), framealpha=0.9)
    fig_b.suptitle(
        'Silymarin intervention: cross-cohort validation across 6 datasets '
        '(t = 300 h, ki = 0.3 = 70% inhibition)',
        fontsize=17, fontweight='bold')

    fig_b.savefig(OUT_FIG_B, dpi=300, bbox_inches='tight')
    plt.close(fig_b)
    print(f'  Saved: {OUT_FIG_B}')

    print('\nDone.')


if __name__ == '__main__':
    main()
