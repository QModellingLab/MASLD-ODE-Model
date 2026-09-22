#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
step4_concordance.py  (GSE89632 version)
================================================================
Step 4 of GSE89632 independent validation pipeline.

Compare log2FC between GSE89632 (validation 2, Illumina microarray)
and GSE126848 (training, RNA-seq).  Two comparisons:
  - NAFL_vs_Normal  (SS in GSE89632 ↔ NAFL in GSE126848)
  - NASH_vs_Normal  (NASH in both)

Output now also includes:
  - Stratified concordance by |log2FC_train| threshold (anticipated
    reviewer request — shows that strong-signal genes agree better)
  - Canonical NAFLD biomarker check (CYP7A1, AKR1B10, COL1A1, etc.)

Inputs:
  - GSE89632_DEG_gene_level.xlsx (from step 3)
  - DEG_pydeseq2_3comparisons.xlsx (GSE126848 DEG)
        EDIT GSE126848_PATH below to point to your local file.

Outputs:
  - Validation_concordance_metrics.xlsx
  - Validation_stratified_concordance.csv
  - Validation_biomarker_check.csv
  - Validation_ODE_gene_check.csv
  - Fig_validation_concordance_scatter.png
================================================================
"""
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy import stats

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
GSE89632_XLSX  = os.path.join(SCRIPT_DIR, 'GSE89632_DEG_gene_level.xlsx')

# EDIT THIS PATH to where DEG_pydeseq2_3comparisons.xlsx lives on your system
GSE126848_PATH = os.path.join(SCRIPT_DIR, 'DEG_pydeseq2_3comparisons.xlsx')

OUT_METRICS      = os.path.join(SCRIPT_DIR, 'Validation_concordance_metrics.xlsx')
OUT_STRATIFIED   = os.path.join(SCRIPT_DIR, 'Validation_stratified_concordance.csv')
OUT_BIOMARKER    = os.path.join(SCRIPT_DIR, 'Validation_biomarker_check.csv')
OUT_ODECHK       = os.path.join(SCRIPT_DIR, 'Validation_ODE_gene_check.csv')
OUT_SCATTER      = os.path.join(SCRIPT_DIR, 'Fig_validation_concordance_scatter.png')

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

# Sensitivity top-7 (NAFL stage drivers from sensitivity analysis)
SENS_TOP_GENES = ['PRKAG2', 'PPARA', 'ADIPOR1', 'MAPK14', 'ADIPOQ', 'LEPR', 'LEP']

# Canonical NAFLD biomarkers (literature-confirmed)
NAFLD_BIOMARKERS = [
    'CYP7A1',  'AKR1B10', 'COL1A1', 'COL3A1', 'MMP9',
    'PNPLA3',  'FASN',    'SREBF1', 'CXCL10', 'LUM',
    'THBS2',   'TGFB1',
]

# Comparison name -> GSE126848 column suffix
COMPARISONS = {
    'NAFL_vs_Normal': 'NAFL',
    'NASH_vs_Normal': 'NASH',
}


def concordance_metrics(merged):
    df = merged.dropna(subset=['log2FC_v', 'log2FC_t']).copy()
    n = len(df)
    if n < 10:
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
        'N_genes_compared': n,
        'Spearman_rho': rho,  'Spearman_p': p_rho,
        'Pearson_r': r,       'Pearson_p': p_r,
        'Directional_agree_strong_pct': same_strong * 100,
        'N_strong_either': len(strong),
        'Directional_agree_sig_train_pct': same_sig * 100,
        'N_sig_train': len(sig_train),
        'Top100_overlap': overlap,
        'Top100_hyper_p': p_hyper,
    }


def stratified_table(merged, name):
    """Stratify by |log2FC_train| threshold. Reviewer-anticipated analysis."""
    rows = []
    df = merged.dropna(subset=['log2FC_v', 'log2FC_t'])
    for thr in [0.0, 0.5, 1.0, 1.5, 2.0, 3.0]:
        sub = df[df['log2FC_t'].abs() > thr]
        if len(sub) < 10:
            rows.append({'comparison': name, 'threshold': thr,
                          'n_genes': len(sub), 'Spearman_rho': np.nan,
                          'directional_pct': np.nan})
            continue
        rho, _ = stats.spearmanr(sub['log2FC_v'], sub['log2FC_t'])
        same = (np.sign(sub['log2FC_v'])
                 == np.sign(sub['log2FC_t'])).mean() * 100
        rows.append({'comparison': name, 'threshold': thr,
                      'n_genes': len(sub), 'Spearman_rho': rho,
                      'directional_pct': same})
    return rows


def main():
    if not os.path.exists(GSE89632_XLSX):
        raise FileNotFoundError(f'{GSE89632_XLSX} not found. Run step 3.')
    if not os.path.exists(GSE126848_PATH):
        raise FileNotFoundError(
            f'GSE126848 DEG file not found: {GSE126848_PATH}\n'
            'Edit GSE126848_PATH at top of this script.')

    gse89 = {n: pd.read_excel(GSE89632_XLSX, n) for n in COMPARISONS}
    gse126 = pd.read_excel(GSE126848_PATH, 'Merged_3grp')
    print(f'GSE126848 merged: {gse126.shape}')

    metrics_all, merged_all = {}, {}
    for name, suffix in COMPARISONS.items():
        print(f'\n[{name}]')
        v = gse89[name][['gene_symbol', 'log2FC', 'padj']].rename(
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

    # --- Stratified concordance (reviewer-anticipated) ---
    print('\n' + '=' * 60)
    print('Stratified concordance by |log2FC_train| threshold')
    print('=' * 60)
    strat_rows = []
    for name, m in merged_all.items():
        rows = stratified_table(m, name)
        strat_rows.extend(rows)
        print(f'\n  {name}:')
        for r in rows:
            print(f'    |log2FC_train|>{r["threshold"]:.1f}: '
                  f'n={r["n_genes"]:>5,}  '
                  f'ρ={r["Spearman_rho"]:>+.3f}  '
                  f'directional={r["directional_pct"]:>5.1f}%')
    pd.DataFrame(strat_rows).to_csv(OUT_STRATIFIED, index=False)
    print(f'  Saved: {OUT_STRATIFIED}')

    # --- ODE 58 gene check ---
    print('\n' + '=' * 60)
    print('ODE 58-gene directional check')
    print('=' * 60)
    ode_rows = []
    for name, m in merged_all.items():
        sub = m[m['gene_symbol'].isin(ODE_GENES)].copy()
        sub['sign_match'] = (np.sign(sub['log2FC_v'])
                              == np.sign(sub['log2FC_t']))
        pct = sub['sign_match'].mean() * 100 if len(sub) else np.nan
        print(f'  {name:<20}  {len(sub):>3}/{len(ODE_GENES)} ODE genes mapped  '
              f'directional = {pct:.1f}%')
        sub['comparison'] = name
        ode_rows.append(sub)
    if ode_rows:
        ode_check = pd.concat(ode_rows, ignore_index=True)
        ode_check.to_csv(OUT_ODECHK, index=False)
        print(f'  Saved: {OUT_ODECHK}')

    # --- Sensitivity top-7 ---
    print('\n' + '=' * 60)
    print('Sensitivity top-7 gene check')
    print('=' * 60)
    for name, m in merged_all.items():
        sub = m[m['gene_symbol'].isin(SENS_TOP_GENES)]
        print(f'\n  --- {name} ---')
        if len(sub) == 0:
            print('    (none mapped)')
            continue
        print(f'  {"Gene":<10}{"FC_GSE89632":>13}{"FC_GSE126848":>15}'
              f'{"sign_match":>13}')
        for _, r in sub.iterrows():
            sm = 'OK' if np.sign(r['log2FC_v']) == np.sign(r['log2FC_t']) else 'X'
            print(f'  {r["gene_symbol"]:<10}{r["log2FC_v"]:>+13.3f}'
                  f'{r["log2FC_t"]:>+15.3f}{sm:>13}')

    # --- Canonical NAFLD biomarker check ---
    print('\n' + '=' * 60)
    print('Canonical NAFLD biomarkers')
    print('=' * 60)
    bio_rows = []
    for name, m in merged_all.items():
        sub = m[m['gene_symbol'].isin(NAFLD_BIOMARKERS)].copy()
        sub['sign_match'] = (np.sign(sub['log2FC_v'])
                              == np.sign(sub['log2FC_t']))
        sub['comparison'] = name
        print(f'\n  --- {name} ---  ({sub["sign_match"].sum()}/{len(sub)} match)')
        for _, r in sub.iterrows():
            sm = 'OK' if r['sign_match'] else 'X'
            print(f'  {r["gene_symbol"]:<8}  '
                  f'GSE89632 = {r["log2FC_v"]:>+6.2f}   '
                  f'GSE126848 = {r["log2FC_t"]:>+6.2f}   {sm}')
        bio_rows.append(sub)
    if bio_rows:
        pd.concat(bio_rows, ignore_index=True).to_csv(OUT_BIOMARKER, index=False)
        print(f'\n  Saved: {OUT_BIOMARKER}')

    # --- Scatter plot ---
    fig, axes = plt.subplots(1, 2, figsize=(11, 5), constrained_layout=True)
    for ax, (name, suffix) in zip(axes, COMPARISONS.items()):
        if name not in merged_all:
            ax.set_title(f'{name}\n(no data)'); continue
        m = merged_all[name]
        ax.scatter(m['log2FC_t'], m['log2FC_v'], s=3, alpha=0.15, c='gray',
                    label=f'all (n={len(m):,})')
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
                     f'Directional = '
                     f'{mtr["Directional_agree_strong_pct"]:.1f}%',
                     transform=ax.transAxes, va='top', fontsize=9,
                     bbox=dict(boxstyle='round', fc='white', ec='gray',
                                alpha=0.9))
        ax.set_xlabel('log$_2$FC  GSE126848 (training)')
        ax.set_ylabel('log$_2$FC  GSE89632 (validation 2)')
        ax.set_title(name)
        ax.legend(loc='lower right', fontsize=8)
    plt.suptitle('Second-cohort validation: GSE89632 (Illumina) vs '
                  'GSE126848 (RNA-seq)', fontsize=12, fontweight='bold')
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
