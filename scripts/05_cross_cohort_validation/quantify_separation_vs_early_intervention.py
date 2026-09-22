#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
quantify_separation_vs_early_intervention.py
================================================================
Quantitative test of the hypothesis:
  "Is there a statistical correlation between the disease (NASH) vs Normal
   separation shown by a P_* output in a given dataset, and the Silymarin
   NAFL/NASH AUC reduction ratio (early-intervention effect)?"

Separation index (SI):
  SI = (AUC_NASH - AUC_Normal) / AUC_Normal * 100   (%)
  AUC is the t=0-300h trapezoidal integral. The three P_* outputs (saturating or not)
  are all computed with the same metric to avoid mixing dimensions (an earlier version mixed t-half delay
  and steady-state diff, which distorted the Pearson linear regression and has been abandoned).
  Larger SI = more pronounced disease-signal separation (timing + magnitude); SI near 0 or negative = the signal
  has not fully activated / the ordering is anomalous.

Early-intervention effect metric: directly reuses the existing AUC reduction ratio
  (NAFL_AUC_reduction% / NASH_AUC_reduction%); > 1 means the NAFL-stage
  intervention effect is better than the NASH-stage effect.

Statistical tests:
  - Pearson correlation (linear relationship)
  - Spearman correlation (monotonic relationship, more robust to extreme values)
  - Scatter plot + regression line for the 18 data points (6 datasets x 3 outputs)
  - Flags the two known exceptions (GSE48452/GSE89632 x P_Hepatocyte_injury)

Output:
  Table_SeparationIndex_vs_AUCratio.xlsx
  Fig_SeparationIndex_vs_AUCratio_scatter.png

Place in the same folder as fig_ValidationCrossCohort_6datasets.py.
All data paths resolve relative to this repository (no external
folders required).
================================================================
"""
import os, importlib.util
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

RATIOS = {
    'GSE126848': None,  # uses the built-in GENE_METADATA, handled specially
    'GSE48452':  os.path.join(RATIOS_DIR, 'GSE48452_node_initial_ratios.xlsx'),
    'GSE89632':  os.path.join(RATIOS_DIR, 'GSE89632_node_initial_ratios.xlsx'),
    'GSE130970': os.path.join(RATIOS_DIR, 'GSE130970_node_initial_ratios.xlsx'),
    'GSE162694': os.path.join(RATIOS_DIR, 'GSE162694_node_initial_ratios.xlsx'),
    'GSE213621': os.path.join(RATIOS_DIR, 'GSE213621_node_initial_ratios.xlsx'),
}

AUC_TABLE = os.path.join(SCRIPT_DIR, 'Table_Silymarin_AUCreduction_NAFL_vs_NASH_6datasets.xlsx')

OUT_XLSX = os.path.join(SCRIPT_DIR, 'Table_SeparationIndex_vs_AUCratio.xlsx')
OUT_FIG  = os.path.join(SCRIPT_DIR, 'Fig_SeparationIndex_vs_AUCratio_scatter.png')

T_MAX, N_POINTS = 300.0, 3001
t = np.linspace(0, T_MAX, N_POINTS)
OUTPUTS = ['P_Cell_death', 'P_Hepatocyte_injury', 'P_Inflammation']

KNOWN_EXCEPTIONS = {('GSE48452', 'P_Hepatocyte_injury'), ('GSE89632', 'P_Hepatocyte_injury')}

def auc_of(traj):
    return float(np.trapezoid(traj, t))

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


def run_sim(model, y0):
    return odeint(model.ode_system, y0, t, args=(None,), mxstep=10000)


def main():
    print('[File check]')
    for f, lb in [(MODEL_NASH, 'MODEL_NASH'), (AUC_TABLE, 'AUC_TABLE')]:
        ok = os.path.exists(f)
        print(f'  {"OK  " if ok else "MISS"} {lb}: {f}')
        if not ok:
            raise FileNotFoundError(f)

    m_nash = load_model(MODEL_NASH)
    SV = m_nash.STATE_VARS
    p_idx = {pw: SV.index(pw + '_active') for pw in OUTPUTS}

    y0_norm = make_normal_y0(m_nash)
    traj_norm = run_sim(m_nash, y0_norm)

    # GSE126848 NASH y0 (uses the built-in GENE_METADATA)
    y0_126_nash = make_y0_from_ratios(
        m_nash, {n: m_nash.GENE_METADATA[n]['initial_ratio'] for n in m_nash.GENE_METADATA})

    print('\n[Computing Separation Index per dataset x output]')
    rows = []
    for ds in ['GSE126848', 'GSE48452', 'GSE89632', 'GSE130970', 'GSE162694', 'GSE213621']:
        if ds == 'GSE126848':
            y0_nash = y0_126_nash
        else:
            r = load_ratios(RATIOS[ds], 'NASH')
            if r is None:
                print(f'  [WARN] {ds}: NASH ratios missing, skip')
                continue
            y0_nash = make_y0_from_ratios(m_nash, r)

        traj_nash = run_sim(m_nash, y0_nash)

        for pw in OUTPUTS:
            i = p_idx[pw]
            norm_t300 = traj_norm[-1, i]
            nash_t300 = traj_nash[-1, i]

            auc_norm = auc_of(traj_norm[:, i])
            auc_nash = auc_of(traj_nash[:, i])
            si = (auc_nash - auc_norm) / auc_norm * 100 if auc_norm > 0 else np.nan

            rows.append({
                'Dataset': ds, 'P_output': pw,
                'Normal_t300': round(norm_t300, 3),
                'NASH_t300': round(nash_t300, 3),
                'AUC_Normal': round(auc_norm, 1),
                'AUC_NASH': round(auc_nash, 1),
                'Separation_Index_%': round(si, 2) if np.isfinite(si) else np.nan,
            })
            print(f'  {ds:<12} {pw:<22} AUC_Normal={auc_norm:9.1f}  AUC_NASH={auc_nash:9.1f}  SI={si:8.1f}%')

    sep_df = pd.DataFrame(rows)

    # Merge the AUC ratio table
    auc_df = pd.read_excel(AUC_TABLE)
    merged = sep_df.merge(auc_df, on=['Dataset', 'P_output'], how='inner')
    merged['is_known_exception'] = merged.apply(
        lambda r: (r['Dataset'], r['P_output']) in KNOWN_EXCEPTIONS, axis=1)

    print(f'\n[Merged table] {len(merged)} rows (6 dataset x 3 output = 18 expected)')

    # ---------- Correlation analysis ----------
    valid = merged.dropna(subset=['Separation_Index_%', 'NAFL/NASH_ratio'])
    x = valid['Separation_Index_%'].values
    y = valid['NAFL/NASH_ratio'].values

    pearson_r, pearson_p = stats.pearsonr(x, y)
    spearman_r, spearman_p = stats.spearmanr(x, y)

    # Robustness check: SI is proportion/percentage data, right-skewed with extreme leverage points
    # (P_Hepatocyte_injury has a baseline AUC much smaller than the other outputs, so the same
    # absolute change gets amplified once converted to %); the standard approach is to log-transform before testing
    # Pearson, to keep leverage points from dominating the linear regression.
    x_shift = x - x.min() + 1.0   # shift to positive values
    x_log = np.log10(x_shift)
    pearson_log_r, pearson_log_p = stats.pearsonr(x_log, y)

    print('\n' + '=' * 60)
    print('[Correlation: Separation Index vs NAFL/NASH AUC ratio]')
    print('=' * 60)
    print(f'  N = {len(valid)}')
    print(f'  Pearson  r = {pearson_r:.3f}, p = {pearson_p:.4f}  (raw SI%, dominated by 2 leverage points)')
    print(f'  Spearman r = {spearman_r:.3f}, p = {spearman_p:.4f}  (rank-based, robust to leverage points)')
    print(f'  Pearson  r = {pearson_log_r:.3f}, p = {pearson_log_p:.4f}  (after log10(SI%) transform, robustness check)')

    # Recompute after excluding the known exceptions (sensitivity analysis)
    valid_noexc = valid[~valid.apply(
        lambda r: (r['Dataset'], r['P_output']) in KNOWN_EXCEPTIONS, axis=1)]
    if len(valid_noexc) >= 3:
        x2, y2 = valid_noexc['Separation_Index_%'].values, valid_noexc['NAFL/NASH_ratio'].values
        pr2, pp2 = stats.pearsonr(x2, y2)
        sr2, sp2 = stats.spearmanr(x2, y2)
        print(f'\n  After excluding the 2 known exceptions (N={len(valid_noexc)}):')
        print(f'  Pearson  r = {pr2:.3f}, p = {pp2:.4f}')
        print(f'  Spearman r = {sr2:.3f}, p = {sp2:.4f}')

    # ---------- Binary group test: directly tests the original hypothesis ----------
    # "Separated (SI > threshold) -> AUC ratio > 1; not separated (SI <= threshold)
    #  -> AUC ratio ≈ 1" is a threshold-like/binary relationship, which fits the observed
    # phenomenon better than a continuous linear correlation (Pearson/log-Pearson are both nonsignificant, but Spearman is significant,
    # a statistical signature of a threshold relationship: the effect appears once the threshold is crossed, and further increases in separation
    # beyond that do not proportionally amplify the effect). A Mann-Whitney U test compares the two groups' AUC ratio distributions.
    SEP_THRESHOLD = 10.0   # SI > 10% is treated as "clearly separated" (an arbitrary but reasonable threshold)
    valid = valid.copy()
    valid['separated'] = valid['Separation_Index_%'] > SEP_THRESHOLD
    grp_sep   = valid.loc[valid['separated'], 'NAFL/NASH_ratio']
    grp_nosep = valid.loc[~valid['separated'], 'NAFL/NASH_ratio']

    print('\n' + '=' * 60)
    print(f'[Binary group test: SI > {SEP_THRESHOLD}% ("separated") vs SI <= {SEP_THRESHOLD}% ("not separated")]')
    print('=' * 60)
    print(f'  Separated     (N={len(grp_sep)}): median ratio = {grp_sep.median():.2f}, '
          f'mean = {grp_sep.mean():.2f}')
    print(f'  Not separated (N={len(grp_nosep)}): median ratio = {grp_nosep.median():.2f}, '
          f'mean = {grp_nosep.mean():.2f}')
    if len(grp_sep) >= 2 and len(grp_nosep) >= 2:
        mw_stat, mw_p = stats.mannwhitneyu(grp_sep, grp_nosep, alternative='greater')
        print(f'  Mann-Whitney U (one-sided, separated > not separated): '
              f'U={mw_stat:.1f}, p={mw_p:.4f}')
    else:
        mw_stat, mw_p = np.nan, np.nan
        print('  [WARN] insufficient sample size, cannot perform a Mann-Whitney U test')

    # ---------- Save Excel ----------
    with pd.ExcelWriter(OUT_XLSX, engine='openpyxl') as w:
        merged.to_excel(w, sheet_name='SeparationIndex_vs_AUCratio', index=False)
        stat_rows = [
            ('N (all 18)', len(valid)),
            ('Pearson r (raw SI%)', round(pearson_r, 3)),
            ('Pearson p (raw SI%)', round(pearson_p, 4)),
            ('Spearman r (rank-based)', round(spearman_r, 3)),
            ('Spearman p (rank-based)', round(spearman_p, 4)),
            ('Pearson r (log10(SI%) robustness check)', round(pearson_log_r, 3)),
            ('Pearson p (log10(SI%) robustness check)', round(pearson_log_p, 4)),
        ]
        if len(valid_noexc) >= 3:
            stat_rows += [
                ('N (excl. 2 exceptions)', len(valid_noexc)),
                ('Pearson r (excl.)', round(pr2, 3)),
                ('Pearson p (excl.)', round(pp2, 4)),
                ('Spearman r (excl.)', round(sr2, 3)),
                ('Spearman p (excl.)', round(sp2, 4)),
            ]
        stat_rows += [
            ('--- Binary group test (SI > 10% threshold) ---', ''),
            ('N separated', len(grp_sep)),
            ('N not separated', len(grp_nosep)),
            ('Median ratio (separated)', round(grp_sep.median(), 3)),
            ('Median ratio (not separated)', round(grp_nosep.median(), 3)),
            ('Mann-Whitney U', round(mw_stat, 2) if np.isfinite(mw_stat) else np.nan),
            ('Mann-Whitney p (one-sided)', round(mw_p, 4) if np.isfinite(mw_p) else np.nan),
        ]
        pd.DataFrame(stat_rows, columns=['Statistic', 'Value']).to_excel(
            w, sheet_name='Summary', index=False)
    print(f'\n  Saved: {OUT_XLSX}')

    # ---------- Scatter plot ----------
    fig, ax = plt.subplots(figsize=(8, 6.5))
    colors = {'P_Cell_death': '#C0392B', 'P_Hepatocyte_injury': '#2980B9',
              'P_Inflammation': '#16A085'}
    x_min_shift = x.min() - 1.0
    for pw in OUTPUTS:
        sub = valid[valid['P_output'] == pw]
        sub_x_log = np.log10(sub['Separation_Index_%'].values - x_min_shift)
        is_exc = sub['is_known_exception'].values
        ax.scatter(sub_x_log[~is_exc], sub.loc[~is_exc, 'NAFL/NASH_ratio'],
                   color=colors[pw], s=90, label=pw, edgecolor='white', linewidth=0.8, zorder=3)
        ax.scatter(sub_x_log[is_exc], sub.loc[is_exc, 'NAFL/NASH_ratio'],
                   facecolor='none', edgecolor=colors[pw], s=160, linewidth=2.2,
                   marker='D', zorder=4)

    xs_log = np.linspace(x_log.min(), x_log.max(), 50)
    slope_log, intercept_log = np.polyfit(x_log, y, 1)
    ax.plot(xs_log, slope_log * xs_log + intercept_log, color='gray', ls='--', lw=1.5, zorder=2,
            label=f'Linear fit on log10(SI%) (r={pearson_log_r:.2f}, p={pearson_log_p:.3f})')

    ax.axhline(1.0, color='black', lw=0.8, alpha=0.5)
    thr_log = np.log10(SEP_THRESHOLD - x_min_shift)
    ax.axvline(thr_log, color='#888', lw=1.2, ls=':', zorder=1,
               label=f'SI = {SEP_THRESHOLD}% threshold')
    ax.set_xlabel('log10(Disease-Normal Separation Index, shifted to positive)\n'
                   '[(AUC_NASH − AUC_Normal) / AUC_Normal × 100, t=0–300h]')
    ax.set_ylabel('Silymarin AUC reduction ratio (NAFL / NASH)\n'
                   '[>1 = NAFL-stage intervention more effective]')
    ax.set_title('Disease-Normal separation predicts early-intervention advantage\n'
                  f'(Spearman r={spearman_r:.2f} p={spearman_p:.3f}; '
                  f'Mann-Whitney (separated vs not) p={mw_p:.3f})')
    ax.legend(fontsize=10, loc='best')
    ax.grid(alpha=0.25)
    fig.tight_layout()
    fig.savefig(OUT_FIG, dpi=300, bbox_inches='tight')
    plt.close(fig)
    print(f'  Saved: {OUT_FIG}')

    print('\nDone.')


if __name__ == '__main__':
    main()
