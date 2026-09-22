#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
step3_DEG.py
================================================================
Step 3 of GSE48452 independent validation pipeline.

For each of three pairwise comparisons (Healthy obese, Steatosis, Nash vs
Control), compute Welch's t-test per probe with Benjamini-Hochberg FDR.
Collapse probe-level results to gene-level by selecting the probe with the
smallest padj per gene.

Samples with bariatric_surgery == 'after surgery' are excluded by default
(post-treatment epigenetic remodeling, not a clean disease-state sample).

Note on method:
  Welch's t-test lacks the empirical-Bayes shrinkage of limma. For small
  groups (Steatosis n=9 in this dataset) it produces conservative padj. For
  ODE-model validation this is acceptable because the validation signal is
  log2FC directional concordance, not stand-alone significance. For final
  manuscript submission you may upgrade to limma (R via rpy2, or the pure-
  Python `inmoose` package).

Inputs (from steps 1 and 2):
  - GSE48452_sample_metadata.csv
  - GSE48452_expression_probe.csv
  - GPL11532_probe_to_gene.csv

Outputs:
  - GSE48452_DEG_probe_level.xlsx     (3 sheets: per-probe statistics)
  - GSE48452_DEG_gene_level.xlsx      (3 sheets: collapsed to gene symbol)
================================================================
"""
import os
import numpy as np
import pandas as pd
from scipy import stats

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
META_CSV  = os.path.join(SCRIPT_DIR, 'GSE48452_sample_metadata.csv')
EXPR_CSV  = os.path.join(SCRIPT_DIR, 'GSE48452_expression_probe.csv')
ANNOT_CSV = os.path.join(SCRIPT_DIR, 'GPL11532_probe_to_gene.csv')
OUT_PROBE = os.path.join(SCRIPT_DIR, 'GSE48452_DEG_probe_level.xlsx')
OUT_GENE  = os.path.join(SCRIPT_DIR, 'GSE48452_DEG_gene_level.xlsx')

COMPARISONS = [
    ('Obese_vs_Normal', 'Healthy obese', 'Control'),
    ('NAFL_vs_Normal',  'Steatosis',     'Control'),
    ('NASH_vs_Normal',  'Nash',          'Control'),
]


def bh_fdr(p):
    """Benjamini-Hochberg FDR adjustment. Returns adjusted p-values."""
    p = np.asarray(p, float)
    n = len(p)
    order = np.argsort(p)
    p_sorted = p[order]
    adj = np.minimum.accumulate(
        (p_sorted * n / (np.arange(n) + 1))[::-1])[::-1]
    adj = np.clip(adj, 0, 1)
    out = np.empty(n)
    out[order] = adj
    return out


def deg_one(expr_use, samples_t, samples_r):
    """Welch's t-test per probe. Returns DataFrame."""
    X_t = expr_use[samples_t].values
    X_r = expr_use[samples_r].values
    mean_t = X_t.mean(axis=1)
    mean_r = X_r.mean(axis=1)
    log2FC = mean_t - mean_r          # already log2 (RMA-normalised input)
    t_stat, p_val = stats.ttest_ind(X_t, X_r, axis=1, equal_var=False)
    padj = bh_fdr(np.nan_to_num(p_val, nan=1.0))
    return pd.DataFrame({
        'probe_id': expr_use.index,
        'mean_treatment': mean_t,
        'mean_reference': mean_r,
        'log2FC': log2FC,
        't_stat': t_stat,
        'pvalue': p_val,
        'padj':   padj,
    })


def collapse_to_gene(probe_df, annot):
    """For each gene, keep the probe with smallest padj."""
    m = probe_df.merge(annot, on='probe_id', how='inner')
    m = m.sort_values('padj').drop_duplicates('gene_symbol', keep='first')
    return m[['gene_symbol', 'probe_id', 'log2FC', 'pvalue',
              'padj', 'mean_treatment', 'mean_reference']]


def main():
    if not (os.path.exists(META_CSV) and os.path.exists(EXPR_CSV)):
        raise FileNotFoundError('Run step1_parse_series_matrix.py first.')

    meta = pd.read_csv(META_CSV)
    expr = pd.read_csv(EXPR_CSV, index_col=0)
    expr.index.name = 'probe_id'

    meta_use = meta[meta['usable']].copy()
    print(f'Usable samples: {len(meta_use)} / {len(meta)} '
          '(excluded "after surgery")')
    print(meta_use['group'].value_counts().to_string())

    expr_use = expr[meta_use['GSM'].values]
    print(f'\nExpression matrix (usable): {expr_use.shape}')

    # --- DEG per comparison ---
    results = {}
    for name, grp_t, grp_r in COMPARISONS:
        s_t = meta_use.loc[meta_use['group'] == grp_t, 'GSM'].tolist()
        s_r = meta_use.loc[meta_use['group'] == grp_r, 'GSM'].tolist()
        print(f'\n[{name}]  n_treat={len(s_t)} ({grp_t})  '
              f'n_ref={len(s_r)} ({grp_r})')
        df = deg_one(expr_use, s_t, s_r)
        n1 = ((df['padj'] < 0.05) & (df['log2FC'].abs() > 1.0)).sum()
        n2 = ((df['padj'] < 0.05) & (df['log2FC'].abs() > 0.5)).sum()
        n3 = (df['pvalue'] < 0.01).sum()
        print(f'  padj<0.05, |log2FC|>1.0 : {n1}')
        print(f'  padj<0.05, |log2FC|>0.5 : {n2}')
        print(f'  pvalue<0.01             : {n3}')
        results[name] = df

    # --- Save probe-level ---
    with pd.ExcelWriter(OUT_PROBE, engine='openpyxl') as xw:
        for n, df in results.items():
            df.to_excel(xw, sheet_name=n, index=False)
    print(f'\nProbe-level DEG saved: {OUT_PROBE}')

    # --- Collapse to gene level (if annotation exists) ---
    if not os.path.exists(ANNOT_CSV):
        print(f'\n[WARN] Annotation file not found: {ANNOT_CSV}')
        print('Run step2_annotate_probes.py first to enable gene-level output.')
        return

    annot = pd.read_csv(ANNOT_CSV)
    annot['probe_id'] = annot['probe_id'].astype(int)
    print(f'\nAnnotation: {len(annot):,} probes -> '
          f'{annot["gene_symbol"].nunique():,} unique genes')

    gene_results = {n: collapse_to_gene(df, annot)
                     for n, df in results.items()}
    with pd.ExcelWriter(OUT_GENE, engine='openpyxl') as xw:
        for n, df in gene_results.items():
            df.to_excel(xw, sheet_name=n, index=False)
            ns = ((df['padj'] < 0.05) & (df['log2FC'].abs() > 1.0)).sum()
            print(f'  [{n}]  N_genes={len(df):,}  '
                  f'sig (padj<0.05, |log2FC|>1.0)={ns}')
    print(f'\nGene-level DEG saved: {OUT_GENE}')

    print('\nDone.')


if __name__ == '__main__':
    main()
