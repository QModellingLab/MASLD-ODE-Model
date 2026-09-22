#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
compute_node_ratios_GSE213621.py
================================================================
GSE213621（Chen et al. 2023）— FPKM matrix
+ GSE213621_series_matrix.txt 內的 !Sample_characteristics_ch1
  （ballooning / inflammation / fibrosis 分組資訊通常放在這裡，
   實際欄位內容請先用 --dry-run 確認）

流程：
  1. 解析 series_matrix.txt 取得每個 GSM 的所有 characteristics 欄位
  2. --dry-run 印出每個 characteristics 欄位的唯一值分布，供核對
  3. 依 GROUP_RULES 規則（在 characteristics 文字中比對關鍵字）
     把樣本分類成 Normal / NAFL / NASH
  4. FPKM -> log2(FPKM+1) -> 各組平均 -> node ratio
================================================================
"""
import os
import re
import sys
import numpy as np
import pandas as pd
from node_ratio_core import load_node_gene_map, print_and_save

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))

# ================================================================
# ★ CONFIG ★
# ================================================================
DATASET = 'GSE213621'

FPKM_TXT      = os.path.join(SCRIPT_DIR, 'GSE213621_FPKMs_allsamples.txt')
SERIES_MATRIX = os.path.join(SCRIPT_DIR, 'GSE213621_series_matrix.txt')
NODE_TABLE    = os.path.join(SCRIPT_DIR, 'node_name_table_hsa04932.xlsx')

ID_TYPE = 'symbol'   # FPKM 矩陣 row index 型態：'symbol' 或 'entrez'

REFERENCE_GROUP = 'Normal'
COMPARISONS = {'NAFL': 'NAFL', 'NASH': 'NASH'}

# 分類規則：依 char_2 (fibrotic stage) 實際內容核對後鎖定如下
#   Control     -> Normal
#   F0F1        -> NAFL  （無/輕度 fibrosis）
#   F2 / F3F4   -> NASH  （顯著 fibrosis，與 NAFL 區隔）
GROUP_RULES = [
    (re.compile(r'control', re.I), 'Normal'),
    (re.compile(r'F0F1', re.I), 'NAFL'),
    (re.compile(r'F2|F3F4', re.I), 'NASH'),
]

OUT_XLSX = os.path.join(SCRIPT_DIR, f'{DATASET}_node_initial_ratios.xlsx')
OUT_FIG  = os.path.join(SCRIPT_DIR, f'{DATASET}_node_initial_ratios_heatmap.png')
# ================================================================


def parse_series_matrix(path):
    """取出 GSM、title、所有 characteristics_ch1 欄位。"""
    gsm_line, title_line = None, None
    char_lines = []
    with open(path, encoding='utf-8', errors='ignore') as fh:
        for line in fh:
            if line.startswith('!Sample_geo_accession'):
                gsm_line = line
            elif line.startswith('!Sample_title'):
                title_line = line
            elif line.startswith('!Sample_characteristics_ch1'):
                char_lines.append(line)

    def split_row(line):
        return [v.strip('"') for v in line.strip().split('\t')[1:]]

    gsms = split_row(gsm_line)
    titles = split_row(title_line) if title_line else [''] * len(gsms)
    df = pd.DataFrame({'GSM': gsms, 'title': titles})
    for i, line in enumerate(char_lines, start=1):
        df[f'char_{i}'] = split_row(line)
    return df


def classify_row(row, char_cols):
    text = ' | '.join(str(row[c]) for c in char_cols) + ' | ' + str(row['title'])
    for pat, grp in GROUP_RULES:
        if pat.search(text):
            return grp
    return None


def main():
    dry_run = '--dry-run' in sys.argv

    print('=' * 65)
    print(f'compute_node_ratios_GSE213621.py  [{DATASET}]'
          + ('  (DRY RUN)' if dry_run else ''))
    print('=' * 65)

    for f, lb in [(FPKM_TXT, 'FPKM_TXT'), (SERIES_MATRIX, 'SERIES_MATRIX'),
                  (NODE_TABLE, 'NODE_TABLE')]:
        ok = os.path.exists(f)
        print(f'  {"OK  " if ok else "MISS"} {lb}: {f}')
        if not ok:
            raise FileNotFoundError(f)

    print('\n[Parsing series_matrix.txt]')
    meta = parse_series_matrix(SERIES_MATRIX)
    char_cols = [c for c in meta.columns if c.startswith('char_')]
    print(f'  {len(meta)} samples, {len(char_cols)} characteristics 欄位')

    if dry_run:
        for c in char_cols:
            print(f'\n  -- {c} 唯一值 --')
            print(meta[c].value_counts(dropna=False).to_string())
        meta.to_csv(os.path.join(SCRIPT_DIR, f'{DATASET}_characteristics_preview.csv'),
                    index=False)
        print(f'\n[Dry run] 完整 characteristics 已存成 '
              f'{DATASET}_characteristics_preview.csv，請核對後調整 GROUP_RULES，'
              '再正式執行（不加 --dry-run）。')
        return

    meta['group'] = meta.apply(lambda r: classify_row(r, char_cols), axis=1)
    print(meta['group'].value_counts(dropna=False).to_string())
    unmatched = meta[meta['group'].isna()]
    if len(unmatched):
        print(f'  [!] {len(unmatched)} 樣本未分類（將被排除）')

    print('\n[Node gene map]')
    node_gene_map = load_node_gene_map(NODE_TABLE, id_type=ID_TYPE)
    n_multi = sum(1 for v in node_gene_map.values() if len(v) > 1)
    print(f'  {len(node_gene_map)} nodes, {n_multi} multi-gene nodes')

    print('\n[Loading FPKM matrix]')
    # FPKM 檔通常是 tab-delimited .txt；自動偵測分隔符號
    fpkm = pd.read_csv(FPKM_TXT, sep=None, engine='python', index_col=0)
    fpkm.index = fpkm.index.astype(str)
    if ID_TYPE == 'entrez':
        fpkm.index = fpkm.index.str.split('.').str[0]
    print(f'  {fpkm.shape[0]:,} genes x {fpkm.shape[1]} samples')

    # ID 比對規則：已用 diagnose_sample_id_matching.py 窮舉驗證，
    # 'Sample_' + title 最後一段 為唯一能 100% 對上矩陣欄名的規則
    # (matched 368/368)，title 格式如 'human liver sample 200707904'
    # -> 矩陣欄名 'Sample_200707904'。
    meta['sample_id'] = 'Sample_' + meta['title'].str.split().str[-1]
    matched_meta = meta[meta['group'].notna() & meta['sample_id'].isin(fpkm.columns)]
    id_col = 'sample_id'
    print(f'  {len(matched_meta)}/{len(meta)} samples matched (Sample_+title最後一段)')
    if len(matched_meta) == 0:
        print('\n  [DEBUG] FPKM 矩陣前 10 個欄名:')
        for c in fpkm.columns[:10]:
            print(f'    {repr(c)}')
        print('\n  [DEBUG] series_matrix 前 10 個 sample_id:')
        for s in meta['sample_id'].tolist()[:10]:
            print(f'    {repr(s)}')
        raise ValueError('Sample_+title最後一段 對不上 FPKM 矩陣欄名，命名規則可能已變動，請重跑 diagnose_sample_id_matching.py')

    log2_fpkm = np.log2(fpkm[matched_meta[id_col].tolist()] + 1.0)
    log2_fpkm.columns = matched_meta[id_col].tolist()

    all_groups = [REFERENCE_GROUP] + list(COMPARISONS.values())
    gm = {}
    for grp in all_groups:
        ids = matched_meta.loc[matched_meta['group'] == grp, id_col].tolist()
        if not ids:
            print(f'  [WARNING] Group "{grp}": 0 samples')
            continue
        gm[grp] = log2_fpkm[ids].mean(axis=1)
        print(f'  {grp}: {len(ids)} samples')
    gm = pd.DataFrame(gm)

    print('\n[Computing node ratios]')
    print_and_save(DATASET, node_gene_map, gm, REFERENCE_GROUP, COMPARISONS,
                   log_scale=True, out_xlsx=OUT_XLSX, out_fig=OUT_FIG)
    print('\nDone.')


if __name__ == '__main__':
    main()
