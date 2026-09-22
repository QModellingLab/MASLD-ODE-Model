#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
step3_DEG.py  (GSE89632 version)
================================================================
Step 3 of GSE89632 independent validation pipeline.

Two pairwise comparisons (no obese-without-steatosis group exists):
  - NAFL_vs_Normal : SS (Simple Steatosis) vs HC (Healthy Control)
  - NASH_vs_Normal : NASH vs HC

Welch's t-test per probe with Benjamini-Hochberg FDR; gene-level
collapse by smallest padj.

Inputs:
  - GSE89632_sample_metadata.csv
  - GSE89632_expression_probe.csv
  - GPL14951_probe_to_gene.csv

Outputs:
  - GSE89632_DEG_probe_level.xlsx
  - GSE89632_DEG_gene_level.xlsx
================================================================
"""
import os
import numpy as np
import pandas as pd
from scipy import stats

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
META_CSV  = os.path.join(SCRIPT_DIR, 'GSE89632_sample_metadata.csv')
EXPR_CSV  = os.path.join(SCRIPT_DIR, 'GSE89632_expression_probe.csv')
ANNOT_CSV = os.path.join(SCRIPT_DIR, 'GPL14951_probe_to_gene.csv')
OUT_PROBE = os.path.join(SCRIPT_DIR, 'GSE89632_DEG_probe_level.xlsx')
OUT_GENE  = os.path.join(SCRIPT_DIR, 'GSE89632_DEG_gene_level.xlsx')

# Comparison name (matched to GSE126848 column suffix) -> (treatment, reference)
# Standard group labels (auto-normalised in step1): 'HC', 'SS', 'NASH'
COMPARISONS = [
    ('NAFL_vs_Normal', 'SS',   'HC'),
    ('NASH_vs_Normal', 'NASH', 'HC'),
]


def bh_fdr(p):
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
    X_t = expr_use[samples_t].values
    X_r = expr_use[samples_r].values
    mean_t = X_t.mean(axis=1)
    mean_r = X_r.mean(axis=1)
    log2FC = mean_t - mean_r          # input expected to be log2-scale
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

    if 'group' not in meta.columns:
        raise SystemExit(
            'No "group" column in metadata.  Edit step1 group detection.')

    # Check that the expression matrix uses log2 scale; if raw intensity, log2-transform
    # (Illumina series matrix files are typically already log2-normalised by submitters)
    if expr.values.max() > 50:
        print(f'[INFO] Detected non-log2 scale (max={expr.values.max():.1f}); '
              'applying log2(x + 1) transformation.')
        expr = np.log2(expr + 1)
    else:
        print(f'[INFO] Expression appears log2-scaled (max={expr.values.max():.2f})')

    print(f'\nGroup counts:')
    print(meta['group'].value_counts().to_string())

    # Align expression columns to metadata samples
    available = [s for s in meta['GSM'].values if s in expr.columns]
    if len(available) < len(meta):
        print(f'[WARN] {len(meta) - len(available)} samples not in expression matrix')
    expr_use = expr[available]
    meta_use = meta.set_index('GSM').loc[available].reset_index()
    print(f'\nUsable samples: {expr_use.shape[1]}')

    # --- DEG per comparison ---
    results = {}
    for name, grp_t, grp_r in COMPARISONS:
        s_t = meta_use.loc[meta_use['group'] == grp_t, 'GSM'].tolist()
        s_r = meta_use.loc[meta_use['group'] == grp_r, 'GSM'].tolist()
        if len(s_t) < 3 or len(s_r) < 3:
            print(f'\n[{name}] SKIP - insufficient samples '
                  f'(treatment={len(s_t)}, reference={len(s_r)})')
            continue
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

    if not results:
        raise SystemExit('No comparisons could be computed.')

    # --- Save probe-level ---
    with pd.ExcelWriter(OUT_PROBE, engine='openpyxl') as xw:
        for n, df in results.items():
            df.to_excel(xw, sheet_name=n, index=False)
    print(f'\nProbe-level DEG saved: {OUT_PROBE}')

    # --- Collapse to gene level ---
    if not os.path.exists(ANNOT_CSV):
        print(f'\n[WARN] Annotation file not found: {ANNOT_CSV}')
        print('Run step2_annotate_probes.py first for gene-level output.')
        return

    annot = pd.read_csv(ANNOT_CSV)
    annot['probe_id'] = annot['probe_id'].astype(str)
    # Ensure probe IDs in DEG tables are string for join
    for n in results:
        results[n]['probe_id'] = results[n]['probe_id'].astype(str)
    print(f'\nAnnotation: {len(annot):,} probes -> '
          f'{annot["gene_symbol"].nunique():,} unique genes')

    gene_results = {n: collapse_to_gene(df, annot) for n, df in results.items()}
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
