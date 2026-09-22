#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
step4_concordance.py
================================================================
Step 4 of GSE48452 independent validation pipeline.

Compare log2 fold changes between GSE48452 (validation, microarray) and
GSE126848 (training, RNA-seq) per pairwise disease-stage comparison.

Validation strategy (per NSTC grant 2.7 and Wang et al. 2014 SEQC):
- Do NOT merge cross-platform expression matrices (technical confound).
- Run DEA independently in each dataset, then compare summary statistics
  (log2FC, padj) at the gene level.

Metrics reported per comparison (Obese/NAFL/NASH vs Normal):
- Spearman & Pearson correlation of log2FC vectors
- Directional agreement % for strong-effect or significant-in-training genes
- Top-100 |log2FC| overlap with hypergeometric p-value
- ODE-58-gene directional check
- Sensitivity-top-7 gene check (stage-dependent drivers from sensitivity
  analysis: PRKAG2, PPARA, ADIPOR1, MAPK14, ADIPOQ, LEPR, LEP)

Inputs:
  - GSE48452_DEG_gene_level.xlsx  (from step 3)
  - DEG_pydeseq2_3comparisons.xlsx (from earlier GSE126848 work)
        EDIT GSE126848_PATH below to point to where this file lives.

Outputs:
  - Validation_concordance_metrics.xlsx
  - Fig_validation_concordance_scatter.png
  - Validation_ODE_gene_check.csv
================================================================
"""
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy import stats

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
GSE48452_XLSX = os.path.join(SCRIPT_DIR, 'GSE48452_DEG_gene_level.xlsx')

# EDIT THIS PATH to where DEG_pydeseq2_3comparisons.xlsx lives on your system
GSE126848_PATH = os.path.join(SCRIPT_DIR, 'DEG_pydeseq2_3comparisons.xlsx')

OUT_METRICS = os.path.join(SCRIPT_DIR, 'Validation_concordance_metrics.xlsx')
OUT_SCATTER = os.path.join(SCRIPT_DIR, 'Fig_validation_concordance_scatter.png')
OUT_ODECHK  = os.path.join(SCRIPT_DIR, 'Validation_ODE_gene_check.csv')

# ODE model 58 genes (KEGG hsa04932)
ODE_GENES = [
    'ADIPOQ','ADIPOR1','AKT3','ATF4','BAX','BCL2L11','BID','CASP3','CASP7','CASP8',
    'CDC42','CEBPA','COX6B2','CXCL8','CYCS','CYP2E1','DDIT3','EIF2AK3','EIF2S1',
    'ERN1','FAS','FASLG','FOS','GSK3A','IKBKB','IL1A','IL6','IL6R','INS','INSR',
    'IRS1','ITCH','LEP','LEPR','MAP3K11','MAP3K5','MAPK14','MAPK8','MLXIP',
    'NDUFC2-KCTD14','NFKB1','NR1H3','P3R3URF-PIK3R3','PKLR','PPARA','PPARG',
    'PRKAG2','RAC1','RXRA','SDHA','SOCS3','SREBF1','TGFB1','TNF','TNFRSF1A',
    'TRAF2','UQCR11','XBP1',
]

# Sensitivity top-7 (stage-dependent drivers, mostly AMPK-PPAR-Adipocytokine axis)
SENS_TOP_GENES = ['PRKAG2', 'PPARA', 'ADIPOR1', 'MAPK14', 'ADIPOQ', 'LEPR', 'LEP']

# Comparison name in GSE48452 -> suffix in GSE126848 Merged_3grp columns
COMPARISONS = {
    'Obese_vs_Normal': 'Obese',
    'NAFL_vs_Normal':  'NAFL',
    'NASH_vs_Normal':  'NASH',
}


def concordance_metrics(merged):
    """Compute summary concordance metrics from joined DEG table.
    Columns expected: log2FC_v, log2FC_t, padj_v, padj_t."""
    df = merged.dropna(subset=['log2FC_v', 'log2FC_t']).copy()
    n_total = len(df)
    if n_total < 10:
        return None

    rho, p_rho = stats.spearmanr(df['log2FC_v'], df['log2FC_t'])
    r, p_r = stats.pearsonr(df['log2FC_v'], df['log2FC_t'])

    strong = df[(df['log2FC_v'].abs() > 0.5) | (df['log2FC_t'].abs() > 0.5)]
    same_strong = (np.sign(strong['log2FC_v'])
                    == np.sign(strong['log2FC_t'])).mean()

    sig_train = df[df['padj_t'] < 0.05]
    same_sig = (np.sign(sig_train['log2FC_v'])
                 == np.sign(sig_train['log2FC_t'])).mean() \
                 if len(sig_train) > 0 else np.nan

    df['_abs_v'] = df['log2FC_v'].abs()
    df['_abs_t'] = df['log2FC_t'].abs()
    top_v = set(df.nlargest(100, '_abs_v')['gene_symbol'])
    top_t = set(df.nlargest(100, '_abs_t')['gene_symbol'])
    overlap = len(top_v & top_t)
    p_hyper = stats.hypergeom.sf(overlap - 1, len(df), 100, 100) \
              if overlap > 0 else 1.0

    return {
        'N_genes_compared': n_total,
        'Spearman_rho': rho, 'Spearman_p': p_rho,
        'Pearson_r': r,      'Pearson_p': p_r,
        'Directional_agree_strong_pct': same_strong * 100,
        'N_strong_either': len(strong),
        'Directional_agree_sig_train_pct': same_sig * 100,
        'N_sig_train': len(sig_train),
        'Top100_overlap': overlap,
        'Top100_hyper_p': p_hyper,
    }


def main():
    if not os.path.exists(GSE48452_XLSX):
        raise FileNotFoundError(f'{GSE48452_XLSX} not found. Run step 3.')
    if not os.path.exists(GSE126848_PATH):
        raise FileNotFoundError(
            f'GSE126848 DEG file not found at: {GSE126848_PATH}\n'
            'Edit GSE126848_PATH at top of this script.')

    gse48 = {n: pd.read_excel(GSE48452_XLSX, n) for n in COMPARISONS}
    gse126 = pd.read_excel(GSE126848_PATH, 'Merged_3grp')
    print(f'GSE126848 merged: {gse126.shape}')

    metrics_all, merged_all = {}, {}
    for name, suffix in COMPARISONS.items():
        print(f'\n[{name}]')
        v = gse48[name][['gene_symbol', 'log2FC', 'padj']].rename(
            columns={'log2FC': 'log2FC_v', 'padj': 'padj_v'})
        t = gse126[['Gene', f'log2FC_{suffix}', f'padj_{suffix}']].rename(
            columns={'Gene': 'gene_symbol',
                      f'log2FC_{suffix}': 'log2FC_t',
                      f'padj_{suffix}': 'padj_t'})
        m = v.merge(t, on='gene_symbol', how='inner')
        print(f'  Common genes: {len(m):,}')
        if len(m) > 10:
            mtr = concordance_metrics(m)
            metrics_all[name] = mtr
            merged_all[name] = m
            print(f'  Spearman rho = {mtr["Spearman_rho"]:+.3f}  '
                  f'(p={mtr["Spearman_p"]:.2e})')
            print(f'  Directional agreement (strong, n={mtr["N_strong_either"]}): '
                  f'{mtr["Directional_agree_strong_pct"]:.1f}%')
            print(f'  Directional agreement (sig in training, '
                  f'n={mtr["N_sig_train"]}): '
                  f'{mtr["Directional_agree_sig_train_pct"]:.1f}%')
            print(f'  Top-100 |log2FC| overlap: {mtr["Top100_overlap"]} '
                  f'(hyper p={mtr["Top100_hyper_p"]:.2e})')

    # --- ODE 58 gene check ---
    print('\n' + '=' * 60)
    print('ODE 58-gene directional check')
    print('=' * 60)
    ode_rows = []
    for name in COMPARISONS:
        if name not in merged_all:
            continue
        m = merged_all[name]
        sub = m[m['gene_symbol'].isin(ODE_GENES)].copy()
        sub['sign_match'] = (np.sign(sub['log2FC_v'])
                              == np.sign(sub['log2FC_t']))
        pct = sub['sign_match'].mean() * 100 if len(sub) else np.nan
        print(f'  {name:<20}  {len(sub):>3}/{len(ODE_GENES)} ODE genes mapped   '
              f'directional agreement = {pct:.1f}%')
        sub['comparison'] = name
        ode_rows.append(sub)
    if ode_rows:
        ode_check = pd.concat(ode_rows, ignore_index=True)
        ode_check.to_csv(OUT_ODECHK, index=False)
        print(f'  Saved: {OUT_ODECHK}')

    # --- Sensitivity top-7 gene check ---
    print('\n' + '=' * 60)
    print('Sensitivity top-7 gene check (NAFL stage drivers)')
    print('=' * 60)
    for name in COMPARISONS:
        if name not in merged_all:
            continue
        m = merged_all[name]
        sub = m[m['gene_symbol'].isin(SENS_TOP_GENES)]
        print(f'\n  --- {name} ---')
        if len(sub) == 0:
            print('    (no top-7 genes mapped; check annotation)')
            continue
        print(f'  {"Gene":<10}{"FC_GSE48452":>13}{"FC_GSE126848":>15}'
              f'{"sign_match":>13}')
        for _, r in sub.iterrows():
            sm = 'OK' if np.sign(r['log2FC_v']) == np.sign(r['log2FC_t']) else 'X'
            print(f'  {r["gene_symbol"]:<10}{r["log2FC_v"]:>+13.3f}'
                  f'{r["log2FC_t"]:>+15.3f}{sm:>13}')

    # --- Scatter plot ---
    fig, axes = plt.subplots(1, 3, figsize=(15, 5), constrained_layout=True)
    for ax, (name, suffix) in zip(axes, COMPARISONS.items()):
        if name not in merged_all:
            ax.set_title(f'{name}\n(no data)'); continue
        m = merged_all[name]
        ax.scatter(m['log2FC_t'], m['log2FC_v'], s=3, alpha=0.15,
                    c='gray', label=f'all (n={len(m):,})')
        ode_m = m[m['gene_symbol'].isin(ODE_GENES)]
        ax.scatter(ode_m['log2FC_t'], ode_m['log2FC_v'], s=40, alpha=0.85,
                    c='steelblue', edgecolor='navy', linewidth=0.5,
                    label=f'ODE 58 (n={len(ode_m)})')
        sens_m = m[m['gene_symbol'].isin(SENS_TOP_GENES)]
        ax.scatter(sens_m['log2FC_t'], sens_m['log2FC_v'], s=80, alpha=1.0,
                    c='crimson', edgecolor='black', linewidth=0.8, marker='D',
                    label=f'Sens top-7 (n={len(sens_m)})')
        for _, r in sens_m.iterrows():
            ax.annotate(r['gene_symbol'], (r['log2FC_t'], r['log2FC_v']),
                         fontsize=8, xytext=(3, 3), textcoords='offset points')
        ax.axhline(0, c='k', lw=0.5); ax.axvline(0, c='k', lw=0.5)
        lim = max(abs(m['log2FC_t']).max(), abs(m['log2FC_v']).max()) * 1.05
        ax.plot([-lim, lim], [-lim, lim], 'k--', lw=0.5, alpha=0.3)
        ax.set_xlim(-lim, lim); ax.set_ylim(-lim, lim)
        mtr = metrics_all.get(name)
        if mtr:
            ax.text(0.05, 0.95,
                     f'Spearman $\\rho$ = {mtr["Spearman_rho"]:.3f}\n'
                     f'Directional agree = '
                     f'{mtr["Directional_agree_strong_pct"]:.1f}%',
                     transform=ax.transAxes, va='top', fontsize=9,
                     bbox=dict(boxstyle='round', fc='white', ec='gray',
                                alpha=0.9))
        ax.set_xlabel('log$_2$FC  GSE126848 (training)')
        ax.set_ylabel('log$_2$FC  GSE48452 (validation)')
        ax.set_title(name)
        ax.legend(loc='lower right', fontsize=8)
    plt.suptitle('Validation: log$_2$FC concordance, '
                  'GSE126848 (RNA-seq) vs GSE48452 (microarray)',
                  fontsize=12, fontweight='bold')
    plt.savefig(OUT_SCATTER, dpi=300, bbox_inches='tight')
    plt.close()
    print(f'\nScatter plot saved: {OUT_SCATTER}')

    # --- Save metrics table ---
    mt_df = pd.DataFrame(metrics_all).T
    with pd.ExcelWriter(OUT_METRICS, engine='openpyxl') as xw:
        mt_df.to_excel(xw, sheet_name='global_metrics')
        for name, m in merged_all.items():
            m.to_excel(xw, sheet_name=f'merged_{name}', index=False)
    print(f'Metrics saved: {OUT_METRICS}')
    print('\nDone.')


if __name__ == '__main__':
    main()
