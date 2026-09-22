#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
node_ratio_core.py
================================================================
三個新資料集（GSE130970 / GSE162694 / GSE213621）共用的核心函式。
與 compute_node_ratios_from_expression.py 邏輯一致，但：
  1. 支援用 entrez_id 或 gene_symbol 對應 node_name_table
  2. 支援 RNA-seq normalization（log2 CPM / log2 TPM / log2 FPKM）
  3. 統一輸出格式，與 fig_ValidationCrossCohort_..._3datasets.py
     讀取的 xlsx 格式（sheet = f'{stage}_ratios', 欄位 node/ratio）完全相容

放在跟三支 dataset script 同一個資料夾。
================================================================
"""
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


def load_node_gene_map(node_table_path, id_type='symbol'):
    """
    讀取 node_name_table_hsa04932.xlsx -> node : [gene id list]
    id_type: 'symbol' 用 gene_symbols 欄；'entrez' 用 entrez_ids 欄；
             'ensembl' 用 ensembl_ids 欄
    """
    df = pd.read_excel(node_table_path, 'NodeNameTable')
    col_map = {'symbol': 'gene_symbols', 'entrez': 'entrez_ids', 'ensembl': 'ensembl_ids'}
    col = col_map[id_type]
    node_map = {}
    for _, r in df.iterrows():
        ids = []
        for g in str(r[col]).split(','):
            g = g.strip().rstrip('*').strip()
            if '(' in g:
                g = g.split('(')[0].strip()
            if id_type in ('entrez', 'ensembl'):
                # entrez/ensembl id 統一成不含小數點/版本號的字串
                g = g.split('.')[0]
            if g and g.lower() != 'nan' and g not in ids:
                ids.append(g)
        node_map[r['ode_variable_name']] = ids
    return node_map


def compute_node_ratios(node_gene_map, gm, disease_grp, ref_grp, log_scale=True):
    """
    gm: DataFrame, index=gene id (symbol/entrez), columns=group name,
        值為「組平均」。
    log_scale=True  -> gm 內已是 log2 值，ratio = 2^(disease-ref)
    log_scale=False -> gm 內是線性值（如 mean TPM），ratio = disease/ref（zero-safe）
    """
    rows = []
    for node, genes in node_gene_map.items():
        ratios_found, genes_found = [], []
        for g in genes:
            if g not in gm.index:
                continue
            if disease_grp not in gm.columns or ref_grp not in gm.columns:
                continue
            d, r = gm.loc[g, disease_grp], gm.loc[g, ref_grp]
            if pd.isna(d) or pd.isna(r):
                continue
            if log_scale:
                ratio = float(2 ** (d - r))
            else:
                if r <= 0:
                    continue
                ratio = float(d / r)
            ratios_found.append(ratio)
            genes_found.append(g)
        ratio = float(np.mean(ratios_found)) if ratios_found else 1.0
        rows.append({
            'node': node,
            'n_genes_total': len(genes),
            'n_genes_found': len(ratios_found),
            'genes_found': ','.join(genes_found) if genes_found else 'NaN',
            'ratio': ratio,
            'log2_ratio': float(np.log2(ratio)) if ratio > 0 else np.nan,
            'mapped': len(ratios_found) > 0,
        })
    return pd.DataFrame(rows)


def direction(x, up_thr=1.05, dn_thr=0.95):
    if x > up_thr:
        return 'up'
    if x < dn_thr:
        return 'dn'
    return 'nc'


def print_and_save(dataset_name, node_gene_map, gm, reference_group,
                    comparisons, log_scale, out_xlsx, out_fig):
    """
    comparisons: {'NAFL':'<group name in gm.columns>', 'NASH': ...}
    輸出格式與 GSE48452/GSE89632 一致，可直接被
    fig_ValidationCrossCohort_..._3datasets.py 的 load_ratios() 讀取。
    """
    all_results = {}
    for label, grp in comparisons.items():
        if grp not in gm.columns:
            print(f'  SKIP {label} ({grp}): group not in matrix')
            continue
        df = compute_node_ratios(node_gene_map, gm, grp, reference_group, log_scale)
        all_results[label] = df
        n_mapped = df['mapped'].sum()
        print(f'  {label} ({grp} vs {reference_group}): {n_mapped}/{len(df)} nodes mapped')

    if not all_results:
        print('  [!] 沒有任何比較成功，請檢查 group 名稱拼字')
        return all_results

    labels = list(all_results.keys())
    summary = all_results[labels[0]][
        ['node', 'n_genes_total', 'n_genes_found', 'genes_found', 'mapped']].copy()
    for label, df in all_results.items():
        summary[f'ratio_{label}'] = df['ratio'].values
        summary[f'log2_{label}'] = df['log2_ratio'].values
        summary[f'dir_{label}'] = df['ratio'].apply(direction)

    print()
    print('=' * 80)
    print(f'{dataset_name}: Node Initial Ratios (condition / {reference_group})')
    print('=' * 80)
    for _, r in summary.iterrows():
        vals = '  '.join(f'{lb}={r[f"ratio_{lb}"]:.3f}({r[f"dir_{lb}"]})' for lb in labels)
        print(f'  {r.node:<14} {vals}   mapped={int(r.n_genes_found)}/{int(r.n_genes_total)}')

    unmapped = summary[~summary['mapped']]
    if len(unmapped):
        print(f'\n  [!] {len(unmapped)} unmapped nodes (ratio 強制設為 1.0):')
        for _, r in unmapped.iterrows():
            print(f'      {r.node}')

    for label in labels:
        n_up = (summary[f'dir_{label}'] == 'up').sum()
        n_dn = (summary[f'dir_{label}'] == 'dn').sum()
        n_nc = (summary[f'dir_{label}'] == 'nc').sum()
        print(f'  {label}: up={n_up} dn={n_dn} nc={n_nc} mapped={summary["mapped"].sum()}/{len(summary)}')

    out_cols = ['node', 'genes_found', 'n_genes_total', 'n_genes_found', 'mapped']
    for lb in labels:
        out_cols += [f'ratio_{lb}', f'log2_{lb}', f'dir_{lb}']

    with pd.ExcelWriter(out_xlsx, engine='openpyxl') as writer:
        summary[out_cols].to_excel(writer, sheet_name='AllNodes', index=False)
        for label, df in all_results.items():
            df.to_excel(writer, sheet_name=f'{label}_ratios', index=False)
    print(f'\n  Excel saved: {out_xlsx}')

    # heatmap
    mapped_summary = summary[summary['mapped']].reset_index(drop=True)
    n_nodes, n_comps = len(mapped_summary), len(labels)
    if n_nodes > 0:
        mat = np.array([[mapped_summary[f'log2_{lb}'].values[i] for lb in labels]
                        for i in range(n_nodes)])
        vmax = min(np.nanmax(np.abs(mat)), 2.5) if np.isfinite(mat).any() else 1.0
        fig, ax = plt.subplots(figsize=(2 + n_comps * 0.8, max(6, n_nodes * 0.22)),
                                constrained_layout=True)
        im = ax.imshow(mat, cmap='RdBu_r', vmin=-vmax, vmax=vmax,
                        aspect='auto', interpolation='nearest')
        ax.set_xticks(range(n_comps)); ax.set_xticklabels(labels, fontsize=9)
        ax.set_yticks(range(n_nodes)); ax.set_yticklabels(mapped_summary['node'].tolist(), fontsize=7)
        ax.set_title(f'{dataset_name}: Node Initial Ratios (log2 vs {reference_group})\n'
                     'Red=up  Blue=dn  White=nc', fontsize=9)
        plt.colorbar(im, ax=ax, label='log2(ratio)', fraction=0.05, pad=0.02)
        plt.savefig(out_fig, dpi=200, bbox_inches='tight')
        plt.close()
        print(f'  Heatmap saved: {out_fig}')

    return all_results
