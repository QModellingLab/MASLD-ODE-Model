#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
quantify_all_outputs_separation_vs_intervention.py
================================================================
Extends the validation of disease-normal separation vs the Silymarin early-intervention effect from
the 3 core P_* outputs (P_Cell_death/P_Hepatocyte_injury/P_Inflammation)
to all 12 P_* outputs, giving 6 datasets x 12 outputs = 72 data points.

For each (dataset, P_*) combination, independently computes:
  1. Separation Index (SI) = (AUC_NASH - AUC_Normal)/AUC_Normal*100, t=0-300h
  2. Silymarin AUC reduction ratio (NAFL / NASH), ki=0.3 (high dose)

Statistical tests:
  - Spearman correlation (continuous trend, secondary evidence)
  - Mann-Whitney U (SI > threshold "separated" vs <= threshold
    "not separated"; the primary test, corresponding to the original hypothesis stated in binary form)

Output:
  Table_AllOutputs_SeparationIndex_vs_AUCratio.xlsx
  Fig_AllOutputs_SeparationIndex_vs_AUCratio_scatter.png

Place in the same folder as fig_ValidationCrossCohort_6datasets.py.
All data paths resolve relative to this repository (no external
folders required).
================================================================
"""
import os, importlib.util, types
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.integrate import odeint
from scipy import stats

# ================================================================
# Fully self-contained: all data paths resolve relative to this repository
# ================================================================
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
REPO_ROOT = os.path.dirname(os.path.dirname(SCRIPT_DIR))
MODEL_DIR = os.path.join(REPO_ROOT, 'data', 'models')
RATIOS_DIR = os.path.join(REPO_ROOT, 'data', 'ratios')
# ================================================================

MODEL_NASH = os.path.join(MODEL_DIR, 'ode_model_pydeseq2_NASH_vs_Normal_v10_mean.py')
MODEL_NAFL = os.path.join(MODEL_DIR, 'ode_model_pydeseq2_NAFL_vs_Normal_v10_mean.py')

RATIOS = {
    'GSE48452':  os.path.join(RATIOS_DIR, 'GSE48452_node_initial_ratios.xlsx'),
    'GSE89632':  os.path.join(RATIOS_DIR, 'GSE89632_node_initial_ratios.xlsx'),
    'GSE130970': os.path.join(RATIOS_DIR, 'GSE130970_node_initial_ratios.xlsx'),
    'GSE162694': os.path.join(RATIOS_DIR, 'GSE162694_node_initial_ratios.xlsx'),
    'GSE213621': os.path.join(RATIOS_DIR, 'GSE213621_node_initial_ratios.xlsx'),
}
DATASETS = ['GSE126848', 'GSE48452', 'GSE89632', 'GSE130970', 'GSE162694', 'GSE213621']

OUT_XLSX = os.path.join(SCRIPT_DIR, 'Table_AllOutputs_SeparationIndex_vs_AUCratio.xlsx')
OUT_FIG  = os.path.join(SCRIPT_DIR, 'Fig_AllOutputs_SeparationIndex_vs_AUCratio_scatter.png')

T_MAX, N_POINTS = 300.0, 3001
t = np.linspace(0, T_MAX, N_POINTS)

# All 12 P_* outputs (consistent with node_name_table / the ODE model)
ALL_OUTPUTS = [
    'P_Hyperinsulinemia', 'P_De_novo_fatty_acid_synthesis', 'P_Improvement_of_NAFLD',
    'P_Development_of_NAFLD', 'P_Adipogenesis', 'P_HCC_proliferation',
    'P_Apoptosis', 'P_Development_of_steatohepatitis', 'P_Cell_death',
    'P_Inflammation', 'P_Fibrosis', 'P_Hepatocyte_injury',
]

# Disease-direction classification (verified from the ODE network topology's upstream driver nodes, not a guess):
#   P_Improvement_of_NAFLD has only two upstream edges, AMPK and PPAR_a
#   (("AMPK","P_Improvement_of_NAFLD"), ("PPAR_a","P_Improvement_of_NAFLD")),
#   both are protective/anti-inflammatory pathways reported in the literature to be downregulated in NASH (PPARA down-regulation is a highly consistent finding).
#   This node's expected direction is therefore opposite to the other 11 outputs: Normal should be > NASH,
#   so the sign must be flipped when computing SI, so that the SI of all 12 outputs consistently represents
#   "the degree of separation toward the direction of disease severity."
NEGATIVE_DIRECTION_OUTPUTS = {'P_Improvement_of_NAFLD'}

SILYMARIN_TARGETS = [
    ('CASP3', 18, 19), ('CASP7', 20, 21), ('CASP8', 22, 23),
    ('CYP2E1', 26, 27), ('IL_8', 58, 59), ('TNFa', 104, 105),
    ('NF_kB', 80, 81), ('TGF_b1', 100, 101),
]

KNOWN_EXCEPTIONS = {('GSE48452', 'P_Hepatocyte_injury'), ('GSE89632', 'P_Hepatocyte_injury')}
SEP_THRESHOLD = 10.0

plt.rcParams.update({
    'font.family': 'sans-serif', 'font.sans-serif': ['Arial', 'DejaVu Sans'],
    'font.size': 12, 'axes.labelsize': 13, 'axes.titlesize': 14,
    'xtick.labelsize': 11, 'ytick.labelsize': 11, 'savefig.dpi': 300,
})


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
            y0[i] = model.PATHWAY_METADATA.get(node, {}).get('initial_level', 100.) \
                if node in model.PATHWAY_METADATA else 1.0
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
                    eq = line.find('='); rhs = line[eq+1:].strip()
                    if rhs.startswith('-('): sign, bs = '-', 1
                    elif rhs.startswith('('): sign, bs = '', 0
                    else: continue
                    depth = 0; ep = None
                    for p, ch in enumerate(rhs[bs:], start=bs):
                        if ch == '(': depth += 1
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


def auc_of(traj):
    return float(np.trapezoid(traj, t))


def auc_pct_reduction(traj_drug, traj_nodrug):
    a0 = auc_of(traj_nodrug)
    a1 = auc_of(traj_drug)
    return (a0 - a1) / a0 * 100 if a0 > 0 else np.nan


def main():
    print('[File check]')
    for f, lb in [(MODEL_NASH, 'MODEL_NASH'), (MODEL_NAFL, 'MODEL_NAFL')]:
        ok = os.path.exists(f)
        print(f'  {"OK  " if ok else "MISS"} {lb}: {f}')
        if not ok:
            raise FileNotFoundError(f)

    m_nash = load_model(MODEL_NASH)
    m_nafl = load_model(MODEL_NAFL)
    SV = m_nash.STATE_VARS
    p_idx = {pw: SV.index(pw + '_active') for pw in ALL_OUTPUTS}

    with open(MODEL_NASH, encoding='utf-8') as fh:
        base_src = fh.read()

    sily_mods = {}
    for ki in [0.3]:
        src = build_sily_source(base_src, ki)
        sily_mods[ki] = make_mod_from_source(src, f'sily_{ki}')

    y0_norm = make_normal_y0(m_nash)
    traj_norm = run_sim(m_nash, y0_norm)

    print('\n[Building y0 for NAFL/NASH per dataset]')
    y0_map = {}  # (dataset, stage) -> y0
    y0_map[('GSE126848', 'NASH')] = make_y0_from_ratios(
        m_nash, {n: m_nash.GENE_METADATA[n]['initial_ratio'] for n in m_nash.GENE_METADATA})
    y0_map[('GSE126848', 'NAFL')] = make_y0_from_ratios(
        m_nash, {n: m_nafl.GENE_METADATA[n]['initial_ratio'] for n in m_nafl.GENE_METADATA})
    for ds in ['GSE48452', 'GSE89632', 'GSE130970', 'GSE162694', 'GSE213621']:
        for stage in ['NASH', 'NAFL']:
            r = load_ratios(RATIOS[ds], stage)
            y0_map[(ds, stage)] = make_y0_from_ratios(m_nash, r) if r else None

    print('\n[Simulating all dataset x stage trajectories (no-drug + ki=0.3)]')
    traj_nodrug = {}
    traj_drug = {}
    for (ds, stage), y0 in y0_map.items():
        if y0 is None:
            continue
        traj_nodrug[(ds, stage)] = run_sim(m_nash, y0)
        traj_drug[(ds, stage)] = run_sim(sily_mods[0.3], y0)
        print(f'  done: {ds} {stage}')

    print('\n[Computing SI and AUC reduction ratio for all 6 dataset x 12 output = 72 combos]')
    rows = []
    for ds in DATASETS:
        if (ds, 'NASH') not in traj_nodrug:
            print(f'  [WARN] {ds}: NASH missing, skip')
            continue
        for pw in ALL_OUTPUTS:
            i = p_idx[pw]
            auc_norm = auc_of(traj_norm[:, i])
            auc_nash_nodrug = auc_of(traj_nodrug[(ds, 'NASH')][:, i])
            si = (auc_nash_nodrug - auc_norm) / auc_norm * 100 if auc_norm > 0 else np.nan

            vals = {}
            for stage in ['NAFL', 'NASH']:
                if (ds, stage) not in traj_nodrug:
                    vals[stage] = np.nan
                    continue
                vals[stage] = auc_pct_reduction(traj_drug[(ds, stage)][:, i],
                                                 traj_nodrug[(ds, stage)][:, i])
            nafl_v, nash_v = vals.get('NAFL', np.nan), vals.get('NASH', np.nan)
            ratio = nafl_v / nash_v if (np.isfinite(nafl_v) and np.isfinite(nash_v) and nash_v > 0) else np.nan

            rows.append({
                'Dataset': ds, 'P_output': pw,
                'Separation_Index_%': round(si, 2) if np.isfinite(si) else np.nan,
                'NAFL_AUC_red%': round(nafl_v, 2) if np.isfinite(nafl_v) else np.nan,
                'NASH_AUC_red%': round(nash_v, 2) if np.isfinite(nash_v) else np.nan,
                'NAFL/NASH_ratio': round(ratio, 3) if np.isfinite(ratio) else np.nan,
            })
        print(f'  {ds}: done')

    df = pd.DataFrame(rows)
    df['is_known_exception'] = df.apply(
        lambda r: (r['Dataset'], r['P_output']) in KNOWN_EXCEPTIONS, axis=1)
    df['output_direction'] = df['P_output'].apply(
        lambda p: 'inverse (protective pathway)' if p in NEGATIVE_DIRECTION_OUTPUTS
        else 'disease-positive')

    print(f'\n[Total rows] {len(df)} (expect 6 x 12 = 72)')

    # P_Improvement_of_NAFLD's upstream (AMPK/PPAR_a) does not overlap at all
    # with the 8 Silymarin targets; the drug has no direct path to this node, and the sign of its SI also differs
    # in meaning from the other 11 outputs, so it is excluded from the main correlation/group test and listed separately for reference.
    df_main = df[~df['P_output'].isin(NEGATIVE_DIRECTION_OUTPUTS)].copy()
    df_excluded_output = df[df['P_output'].isin(NEGATIVE_DIRECTION_OUTPUTS)].copy()

    print(f'\n[Main analysis] after excluding {NEGATIVE_DIRECTION_OUTPUTS}: {len(df_main)} rows (expect 6 x 11 = 66)')
    if len(df_excluded_output):
        print(f'\n[P_Improvement_of_NAFLD listed separately, not included in the main test]')
        print(df_excluded_output[['Dataset', 'Separation_Index_%', 'NAFL/NASH_ratio']]
              .to_string(index=False))

    valid = df_main.dropna(subset=['Separation_Index_%', 'NAFL/NASH_ratio'])
    x = valid['Separation_Index_%'].values
    y = valid['NAFL/NASH_ratio'].values

    spearman_r, spearman_p = stats.spearmanr(x, y)
    pearson_r, pearson_p = stats.pearsonr(x, y)

    valid = valid.copy()
    valid['separated'] = valid['Separation_Index_%'] > SEP_THRESHOLD
    grp_sep = valid.loc[valid['separated'], 'NAFL/NASH_ratio']
    grp_nosep = valid.loc[~valid['separated'], 'NAFL/NASH_ratio']

    print('\n' + '=' * 65)
    print(f'[Main analysis: 11 disease-positive outputs] N = {len(valid)}')
    print('=' * 65)
    print(f'  Pearson  r = {pearson_r:.3f}, p = {pearson_p:.4f}')
    print(f'  Spearman r = {spearman_r:.3f}, p = {spearman_p:.4f}')
    print(f'\n  Separated     (N={len(grp_sep)}): median ratio = {grp_sep.median():.2f}, mean = {grp_sep.mean():.2f}')
    print(f'  Not separated (N={len(grp_nosep)}): median ratio = {grp_nosep.median():.2f}, mean = {grp_nosep.mean():.2f}')
    if len(grp_sep) >= 2 and len(grp_nosep) >= 2:
        mw_stat, mw_p = stats.mannwhitneyu(grp_sep, grp_nosep, alternative='greater')
        print(f'  Mann-Whitney U (one-sided): U={mw_stat:.1f}, p={mw_p:.4f}')
    else:
        mw_stat, mw_p = np.nan, np.nan

    # Excluding known exceptions
    valid_noexc = valid[~valid.apply(lambda r: (r['Dataset'], r['P_output']) in KNOWN_EXCEPTIONS, axis=1)]
    grp_sep2 = valid_noexc.loc[valid_noexc['separated'], 'NAFL/NASH_ratio']
    grp_nosep2 = valid_noexc.loc[~valid_noexc['separated'], 'NAFL/NASH_ratio']
    if len(grp_sep2) >= 2 and len(grp_nosep2) >= 2:
        mw_stat2, mw_p2 = stats.mannwhitneyu(grp_sep2, grp_nosep2, alternative='greater')
        print(f'\n  After excluding known exceptions (N={len(valid_noexc)}): Mann-Whitney p = {mw_p2:.4f}')
    else:
        mw_stat2, mw_p2 = np.nan, np.nan

    # ---------- Save Excel ----------
    with pd.ExcelWriter(OUT_XLSX, engine='openpyxl') as w:
        df.to_excel(w, sheet_name='AllOutputs_72combos', index=False)
        stat_rows = [
            ('Note', 'P_Improvement_of_NAFLD excluded from main test: '
                      'upstream(AMPK/PPAR_a) has no overlap with Silymarin targets, '
                      'inverse disease direction. See its own rows for reference.'),
            ('N (main, 11 disease-positive outputs x 6 datasets)', len(valid)),
            ('Pearson r', round(pearson_r, 3)), ('Pearson p', round(pearson_p, 4)),
            ('Spearman r', round(spearman_r, 3)), ('Spearman p', round(spearman_p, 4)),
            ('--- Binary group test (SI > 10%) ---', ''),
            ('N separated', len(grp_sep)), ('N not separated', len(grp_nosep)),
            ('Median ratio (separated)', round(grp_sep.median(), 3)),
            ('Median ratio (not separated)', round(grp_nosep.median(), 3)),
            ('Mann-Whitney p (one-sided)', round(mw_p, 5) if np.isfinite(mw_p) else np.nan),
            ('Mann-Whitney p excl. known exceptions', round(mw_p2, 5) if np.isfinite(mw_p2) else np.nan),
        ]
        pd.DataFrame(stat_rows, columns=['Statistic', 'Value']).to_excel(
            w, sheet_name='Summary', index=False)
    print(f'\n  Saved: {OUT_XLSX}')

    # ---------- Scatter plot ----------
    fig, ax = plt.subplots(figsize=(9, 7))
    cmap = plt.get_cmap('tab20')
    main_outputs = [pw for pw in ALL_OUTPUTS if pw not in NEGATIVE_DIRECTION_OUTPUTS]
    colors = {pw: cmap(i / len(main_outputs)) for i, pw in enumerate(main_outputs)}
    x_min_shift = x.min() - 1.0
    for pw in main_outputs:
        sub = valid[valid['P_output'] == pw]
        if len(sub) == 0:
            continue
        sub_x_log = np.log10(sub['Separation_Index_%'].values - x_min_shift)
        is_exc = sub['is_known_exception'].values
        ax.scatter(sub_x_log[~is_exc], sub.loc[~is_exc, 'NAFL/NASH_ratio'],
                   color=colors[pw], s=55, label=pw.replace('P_', ''),
                   edgecolor='white', linewidth=0.5, zorder=3)
        if is_exc.any():
            ax.scatter(sub_x_log[is_exc], sub.loc[is_exc, 'NAFL/NASH_ratio'],
                       facecolor='none', edgecolor=colors[pw], s=120, linewidth=2.0,
                       marker='D', zorder=4)

    thr_log = np.log10(SEP_THRESHOLD - x_min_shift)
    ax.axvline(thr_log, color='#888', lw=1.2, ls=':', zorder=1,
               label=f'SI = {SEP_THRESHOLD}% threshold')
    ax.axhline(1.0, color='black', lw=0.8, alpha=0.5)
    ax.set_xlabel('log10(Disease-Normal Separation Index, shifted to positive)')
    ax.set_ylabel('Silymarin AUC reduction ratio (NAFL / NASH)\n[>1 = NAFL-stage more effective]')
    ax.set_title('11 disease-positive P_* outputs x 6 datasets (N=66)\n'
                  'P_Improvement_of_NAFLD excluded (inverse direction; see table)\n'
                  f'Spearman r={spearman_r:.2f} p={spearman_p:.3f}; '
                  f'Mann-Whitney (separated vs not) p={mw_p:.4f}')
    ax.legend(fontsize=8, loc='upper left', ncol=2, bbox_to_anchor=(1.01, 1.0))
    ax.grid(alpha=0.25)
    fig.tight_layout()
    fig.savefig(OUT_FIG, dpi=300, bbox_inches='tight')
    plt.close(fig)
    print(f'  Saved: {OUT_FIG}')

    print('\nDone.')


if __name__ == '__main__':
    main()
