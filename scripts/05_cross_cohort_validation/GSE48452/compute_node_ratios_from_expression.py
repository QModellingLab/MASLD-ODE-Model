#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
compute_node_ratios_from_expression.py
================================================================
從基因表現 matrix 計算 ODE 模型的 node initial ratio
（condition / Normal 的 fold change）

適用於：
  - GSE48452（Affymetrix HuGene 1.1 ST，RMA log2 normalized）
  - GSE89632（Illumina HT-12 DASL，log2 normalized）
  - 任何 log2 scale 的基因表現 matrix

輸入：
  1. expression_matrix.csv   — rows = probe/gene, cols = sample GSM IDs
                               值為 log2 normalized intensity
  2. sample_metadata.csv     — 必要欄位: GSM, group
                               選用欄位: usable（True/False）
  3. probe_to_gene.csv       — 欄位: probe_id, gene_symbol
                               （如果 expression matrix 已是 gene level 則可略）
  4. node_name_table_hsa04932.xlsx — ODE node → gene list 對應表

輸出：
  node_ratios_<DATASET>.xlsx   — 所有 nodes 的 ratio（各比較 vs Reference）
  node_ratios_<DATASET>.png    — Heatmap

設定（修改下方 CONFIG 區塊）：
  DATASET           : 資料集名稱（用於輸出檔名）
  EXPR_CSV          : expression matrix 路徑
  META_CSV          : sample metadata 路徑
  ANNOT_CSV         : probe→gene annotation 路徑（gene level 可設 None）
  NODE_TABLE        : node_name_table_hsa04932.xlsx 路徑
  REFERENCE_GROUP   : 對照組的 group 名稱
  COMPARISONS       : 要計算的比較（輸出標籤 → metadata group 名稱）
  USE_USABLE_FLAG   : 是否使用 metadata 的 usable 欄位過濾樣本
  IS_GENE_LEVEL     : expression matrix 是否已是 gene level（不需 probe 對應）

================================================================
"""
import os
import sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))

# ================================================================
# ★ CONFIG — 修改這裡 ★
# ================================================================
DATASET = 'GSE48452'   # 輸出檔名用，可改為 GSE89632 等

# 輸入檔案路徑（相對或絕對路徑均可）
EXPR_CSV   = os.path.join(SCRIPT_DIR, 'GSE48452_expression_probe.csv')
META_CSV   = os.path.join(SCRIPT_DIR, 'GSE48452_sample_metadata.csv')
ANNOT_CSV  = os.path.join(SCRIPT_DIR, 'GPL11532_probe_to_gene.csv')   # probe→gene
NODE_TABLE = os.path.join(SCRIPT_DIR, 'node_name_table_hsa04932.xlsx')

# 對照組名稱（metadata group 欄位中的值）
REFERENCE_GROUP = 'Control'

# 要計算的比較：輸出標籤 → metadata group 名稱
# 格式：{'輸出名稱': 'metadata中的group名稱'}
COMPARISONS = {
    'Obese':     'Healthy obese',
    'Steatosis': 'Steatosis',
    'NASH':      'Nash',
}

# 是否使用 usable 欄位（True = 排除 usable==False 的樣本）
USE_USABLE_FLAG = True

# expression matrix 是否已是 gene level（True = 不需要 ANNOT_CSV）
IS_GENE_LEVEL = False

# 輸出檔案
OUT_XLSX = os.path.join(SCRIPT_DIR, f'node_ratios_{DATASET}.xlsx')
OUT_FIG  = os.path.join(SCRIPT_DIR, f'node_ratios_{DATASET}_heatmap.png')
# ================================================================


def load_node_gene_map(node_table_path):
    """讀取 node_name_table，建立 node → gene list 對應。"""
    df = pd.read_excel(node_table_path, 'NodeNameTable')
    node_map = {}
    for _, r in df.iterrows():
        genes = []
        for g in str(r['gene_symbols']).split(','):
            g = g.strip().rstrip('*').strip()
            if '(' in g:
                g = g.split('(')[0].strip()
            if g and g.lower() != 'nan' and g not in genes:
                genes.append(g)
        node_map[r['ode_variable_name']] = genes
    return node_map


def build_gene_group_means(expr_csv, meta_csv, annot_csv,
                            reference_group, comparisons,
                            use_usable_flag, is_gene_level):
    """
    從 expression matrix 計算各 group 的 per-gene 平均表現值。
    回傳 DataFrame，index=gene_symbol，columns=group names。
    """
    print('[1/4] 讀取 expression matrix...')
    expr = pd.read_csv(expr_csv, index_col=0)
    expr.index = expr.index.astype(str)
    print(f'  Expression matrix: {expr.shape[0]:,} probes × {expr.shape[1]} samples')

    print('[2/4] 讀取 sample metadata...')
    meta = pd.read_csv(meta_csv)
    if use_usable_flag and 'usable' in meta.columns:
        n_before = len(meta)
        meta = meta[meta['usable'] == True]
        print(f'  Filtered by usable: {len(meta)}/{n_before} samples kept')

    all_groups = [reference_group] + list(comparisons.values())
    print(f'  Groups needed: {all_groups}')
    print(f'  Groups in metadata: {sorted(meta["group"].unique())}')

    # 確認所有 group 都存在
    missing_groups = [g for g in all_groups if g not in meta['group'].values]
    if missing_groups:
        print(f'  [WARNING] Missing groups: {missing_groups}')
        print('  請確認 COMPARISONS 和 REFERENCE_GROUP 的拼字與 metadata 一致')

    if not is_gene_level:
        print('[3/4] 讀取 probe→gene annotation 並 collapse to gene level...')
        annot = pd.read_csv(annot_csv)
        annot['probe_id'] = annot['probe_id'].astype(str)
        # 確認 gene_symbol 欄位
        if 'gene_symbol' not in annot.columns:
            raise ValueError(f'annotation 需要 gene_symbol 欄位，現有: {annot.columns.tolist()}')
        # merge
        e2 = expr.merge(annot, left_index=True, right_on='probe_id', how='inner')
        print(f'  Probes matched: {len(e2):,} / {len(expr):,}')
        # collapse: per gene, 取平均（也可改 max）
        sample_cols = [c for c in expr.columns if c in e2.columns]
        ge = e2.groupby('gene_symbol')[sample_cols].mean()
        print(f'  Gene-level matrix: {ge.shape[0]:,} genes × {ge.shape[1]} samples')
    else:
        print('[3/4] Expression matrix 已是 gene level，跳過 probe collapse...')
        ge = expr.copy()
        ge.index.name = 'gene_symbol'

    print('[4/4] 計算各 group 平均...')
    gm = {}
    for grp in all_groups:
        gsms = meta.loc[meta['group'] == grp, 'GSM'].tolist()
        gsms_present = [g for g in gsms if g in ge.columns]
        if not gsms_present:
            print(f'  [WARNING] Group "{grp}": 0 samples found in expression matrix')
            continue
        gm[grp] = ge[gsms_present].mean(axis=1)
        print(f'  {grp}: {len(gsms_present)} samples')

    return pd.DataFrame(gm)


def compute_node_ratios(node_gene_map, gm, disease_grp, ref_grp):
    """
    計算每個 ODE node 的 ratio = 2^(disease_mean - ref_mean)，
    多基因 node 取算術平均。
    """
    rows = []
    for node, genes in node_gene_map.items():
        ratios_found = []
        genes_found  = []
        for g in genes:
            if g in gm.index:
                if disease_grp in gm.columns and ref_grp in gm.columns:
                    r = float(2 ** (gm.loc[g, disease_grp] - gm.loc[g, ref_grp]))
                    ratios_found.append(r)
                    genes_found.append(g)
        ratio  = float(np.mean(ratios_found)) if ratios_found else 1.0
        mapped = len(ratios_found) > 0
        rows.append({
            'node':          node,
            'n_genes_total': len(genes),
            'n_genes_found': len(ratios_found),
            'genes_found':   ','.join(genes_found) if genes_found else 'NaN',
            'ratio':         ratio,
            'log2_ratio':    float(np.log2(ratio)) if ratio > 0 else np.nan,
            'mapped':        mapped,
        })
    return pd.DataFrame(rows)


def direction(x, up_thr=1.05, dn_thr=0.95):
    if x > up_thr:   return 'up'
    if x < dn_thr:   return 'dn'
    return 'nc'


def main():
    print('=' * 65)
    print(f'compute_node_ratios_from_expression.py  [{DATASET}]')
    print('=' * 65)

    # --- File check ---
    required = [EXPR_CSV, META_CSV, NODE_TABLE]
    if not IS_GENE_LEVEL:
        required.append(ANNOT_CSV)
    for f in required:
        exists = os.path.exists(f)
        print(f'  {"OK  " if exists else "MISS"} {os.path.basename(f)}: {f}')
        if not exists:
            raise FileNotFoundError(f'Required file not found: {f}')

    # --- Load node gene map ---
    print('\n[Node gene map]')
    node_gene_map = load_node_gene_map(NODE_TABLE)
    n_multi = sum(1 for v in node_gene_map.values() if len(v) > 1)
    print(f'  {len(node_gene_map)} nodes, {n_multi} multi-gene nodes')

    # --- Build gene group means ---
    print()
    gm = build_gene_group_means(
        EXPR_CSV, META_CSV, ANNOT_CSV if not IS_GENE_LEVEL else None,
        REFERENCE_GROUP, COMPARISONS, USE_USABLE_FLAG, IS_GENE_LEVEL)

    print(f'\n  Gene×group matrix: {gm.shape}')
    print(f'  Groups: {list(gm.columns)}')

    # --- Compute node ratios for each comparison ---
    print('\n[Computing node ratios]')
    all_results = {}
    for label, grp in COMPARISONS.items():
        if grp not in gm.columns:
            print(f'  SKIP {label} ({grp}): group not in expression matrix')
            continue
        df = compute_node_ratios(node_gene_map, gm, grp, REFERENCE_GROUP)
        all_results[label] = df
        n_mapped = df['mapped'].sum()
        print(f'  {label} ({grp} vs {REFERENCE_GROUP}): '
              f'{n_mapped}/{len(df)} nodes mapped')

    # --- Build summary table ---
    summary = all_results[list(all_results.keys())[0]][
        ['node','n_genes_total','n_genes_found','genes_found','mapped']].copy()

    for label, df in all_results.items():
        summary[f'ratio_{label}']    = df['ratio'].values
        summary[f'log2_{label}']     = df['log2_ratio'].values
        summary[f'dir_{label}']      = df['ratio'].apply(direction)

    # --- Console output ---
    print()
    print('=' * 90)
    print('Node Initial Ratios (condition / Normal baseline)')
    print('  up = ratio > 1.05  |  dn = ratio < 0.95  |  nc = 0.95–1.05')
    print('=' * 90)

    labels = list(all_results.keys())
    header = f'{"Node":<12} '
    for lb in labels:
        header += f'{lb[:8]:>9} '
    header += '  '
    for lb in labels:
        header += f'{lb[:4]:>5} '
    header += ' mapped'
    print(header)
    print('-' * len(header))

    for _, r in summary.iterrows():
        line = f'{r.node:<12} '
        for lb in labels:
            line += f'{r[f"ratio_{lb}"]:>9.3f} '
        line += '  '
        for lb in labels:
            line += f'{r[f"dir_{lb}"]:>5} '
        mapped_str = f'{int(r.n_genes_found)}/{int(r.n_genes_total)}'
        line += f' {mapped_str}'
        print(line)

    # Unmapped
    unmapped = summary[~summary['mapped']]
    if len(unmapped):
        print(f'\n[!] {len(unmapped)} unmapped nodes (ratio imputed to 1.0):')
        for _, r in unmapped.iterrows():
            print(f'    {r.node} — genes: {r.genes_found}')

    # Summary stats
    print()
    print('=' * 50)
    for label in labels:
        n_up = (summary[f'dir_{label}'] == 'up').sum()
        n_dn = (summary[f'dir_{label}'] == 'dn').sum()
        n_nc = (summary[f'dir_{label}'] == 'nc').sum()
        print(f'{label}: up={n_up}  dn={n_dn}  nc={n_nc}  '
              f'mapped={summary["mapped"].sum()}/{len(summary)}')

    # --- Save Excel ---
    out_cols = ['node', 'genes_found', 'n_genes_total', 'n_genes_found', 'mapped']
    for lb in labels:
        out_cols += [f'ratio_{lb}', f'log2_{lb}', f'dir_{lb}']

    with pd.ExcelWriter(OUT_XLSX, engine='openpyxl') as writer:
        summary[out_cols].to_excel(writer, sheet_name='AllNodes', index=False)
        for label, df in all_results.items():
            df.to_excel(writer, sheet_name=f'{label}_ratios', index=False)
    print(f'\nExcel saved: {OUT_XLSX}')

    # --- Heatmap ---
    mapped_summary = summary[summary['mapped']].reset_index(drop=True)
    n_nodes = len(mapped_summary)
    n_comps = len(labels)

    mat = np.array([[mapped_summary[f'log2_{lb}'].values[i]
                     for lb in labels]
                    for i in range(n_nodes)])

    vmax = min(np.nanmax(np.abs(mat)), 2.5)

    fig, ax = plt.subplots(figsize=(2 + n_comps * 0.8,
                                    max(6, n_nodes * 0.22)),
                           constrained_layout=True)
    im = ax.imshow(mat, cmap='RdBu_r', vmin=-vmax, vmax=vmax,
                   aspect='auto', interpolation='nearest')

    ax.set_xticks(range(n_comps))
    ax.set_xticklabels(labels, fontsize=9)
    ax.set_yticks(range(n_nodes))
    ax.set_yticklabels(mapped_summary['node'].tolist(), fontsize=7)
    ax.set_title(f'{DATASET}: Node Initial Ratios (log2 vs {REFERENCE_GROUP})\n'
                 'Red=up  Blue=dn  White=nc', fontsize=9)
    plt.colorbar(im, ax=ax, label='log2(ratio)', fraction=0.05, pad=0.02)

    plt.savefig(OUT_FIG, dpi=200, bbox_inches='tight')
    plt.close()
    print(f'Heatmap saved: {OUT_FIG}')
    print('\nDone.')


if __name__ == '__main__':
    main()
