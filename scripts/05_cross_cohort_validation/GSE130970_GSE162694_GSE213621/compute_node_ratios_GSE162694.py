#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
compute_node_ratios_GSE162694.py
================================================================
GSE162694 — raw_counts.csv + GSE162694_series_matrix.txt

與 GSE130970 / GSE213621 腳本邏輯一致：分組資訊一律從 GEO
series_matrix.txt 的 !Sample_characteristics_ch1（family file 的
metadata 內容）讀取，不靠解析/猜測 sample title 字串。

流程：
  1. 解析 series_matrix.txt 取得每個 GSM 的 title + 所有
     !Sample_characteristics_ch1 欄位
  2. --dry-run 印出每個 characteristics 欄位的唯一值分布，供核對
  3. 依 GROUP_RULES（在 --dry-run 核對實際內容後鎖定）分類成
     Normal / NAFL / NASH
  4. raw counts -> CPM -> log2(CPM+1) -> 各組平均 -> node ratio
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
DATASET = 'GSE162694'

COUNTS_CSV    = os.path.join(SCRIPT_DIR, 'GSE162694_raw_counts.csv')
SERIES_MATRIX = os.path.join(SCRIPT_DIR, 'GSE162694_series_matrix.txt')
NODE_TABLE    = os.path.join(SCRIPT_DIR, 'node_name_table_hsa04932.xlsx')

ID_TYPE = 'ensembl'   # 已用 diagnose_gse162694_gene_id.py 驗證確認：
                       # raw_counts.csv row index 為 Ensembl ID（ENSG...），
                       # 與 node table ensembl_ids 重疊 145/152，非 symbol/entrez

REFERENCE_GROUP = 'Normal'
COMPARISONS = {'NAFL': 'NAFL', 'NASH': 'NASH'}

# 分類規則：依查證文獻鎖定（非猜測）。
# 來源：Identification of disease-related genes and construction of a
#       gene co-expression database in NAFLD (Frontiers in Genetics,
#       PMC10083285)，該文分析同一個 GSE162694 dataset，明確記載
#       對此 dataset 的比較方式是 "control vs. higher-grade fibrosis"
#       （即用 fibrosis stage 分組，而非 NAS-based NAFL/NASH，因為本
#       dataset 沒有獨立的 steatosis/ballooning grade，只有 fibrosis
#       stage + NAS total score）。
# 與本資料夾 --dry-run 印出的 char_4 (fibrosis stage) 欄位完全吻合：
#   normal liver histology (31) / 0 (35) / 1 (30) / 2 (27) / 3 (8) / 4 (12)
# 為與 GSE213621（同樣只有 fibrosis stage）方法學對齊，採用：
#   Normal : fibrosis stage = "normal liver histology"
#   NAFL   : fibrosis stage = 0 或 1（無/輕度 fibrosis）
#   NASH   : fibrosis stage = 2、3 或 4（顯著 fibrosis）
FIBROSIS_COL = 'char_4'

FIBROSIS_GROUP_MAP = {
    'normal liver histology': 'Normal',
    '0': 'NAFL', '1': 'NAFL',
    '2': 'NASH', '3': 'NASH', '4': 'NASH',
}

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
    val = str(row.get(FIBROSIS_COL, ''))
    m = re.search(r':\s*(.+)$', val)
    key = m.group(1).strip().lower() if m else None
    return FIBROSIS_GROUP_MAP.get(key)


def main():
    dry_run = '--dry-run' in sys.argv

    print('=' * 65)
    print(f'compute_node_ratios_GSE162694.py  [{DATASET}]'
          + ('  (DRY RUN)' if dry_run else ''))
    print('=' * 65)

    for f, lb in [(COUNTS_CSV, 'COUNTS_CSV'), (SERIES_MATRIX, 'SERIES_MATRIX'),
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
        print('\n  -- 前 30 筆 title 原始內容（供參考比對）--')
        for tt in meta['title'].tolist()[:30]:
            print(f'    {tt}')
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

    print('\n[Loading raw counts]')
    counts = pd.read_csv(COUNTS_CSV, index_col=0)
    counts.index = counts.index.astype(str)
    if ID_TYPE in ('entrez', 'ensembl'):
        counts.index = counts.index.str.split('.').str[0]
    print(f'  {counts.shape[0]:,} genes x {counts.shape[1]} samples')

    # ID 比對規則：已用 diagnose_sample_id_matching.py 窮舉驗證，
    # 'title 最後一段(空白後)' 為唯一能 100% 對上矩陣欄名的規則
    # (matched 143/143)，title 格式如 'nash1_F0 548nash1'。
    meta['sample_id'] = meta['title'].str.split().str[-1]
    matched_meta = meta[meta['group'].notna() & meta['sample_id'].isin(counts.columns)]
    id_col = 'sample_id'
    print(f'  {len(matched_meta)}/{len(meta)} samples matched (title 最後一段)')
    if len(matched_meta) == 0:
        print('\n  [DEBUG] counts 矩陣前 10 個欄名:')
        for c in counts.columns[:10]:
            print(f'    {repr(c)}')
        print('\n  [DEBUG] series_matrix 前 10 個 sample_id:')
        for s in meta['sample_id'].tolist()[:10]:
            print(f'    {repr(s)}')
        raise ValueError('title 最後一段對不上 counts 矩陣欄名，命名規則可能已變動，請重跑 diagnose_sample_id_matching.py')

    # raw counts -> CPM -> log2(CPM+1)
    counts_used = counts[matched_meta[id_col].tolist()]
    libsize = counts_used.sum(axis=0)
    cpm = counts_used.div(libsize, axis=1) * 1e6
    log2_cpm = np.log2(cpm + 1.0)
    log2_cpm.columns = matched_meta[id_col].tolist()

    all_groups = [REFERENCE_GROUP] + list(COMPARISONS.values())
    gm = {}
    for grp in all_groups:
        ids = matched_meta.loc[matched_meta['group'] == grp, id_col].tolist()
        if not ids:
            print(f'  [WARNING] Group "{grp}": 0 samples')
            continue
        gm[grp] = log2_cpm[ids].mean(axis=1)
        print(f'  {grp}: {len(ids)} samples')
    gm = pd.DataFrame(gm)

    print('\n[Computing node ratios]')
    print_and_save(DATASET, node_gene_map, gm, REFERENCE_GROUP, COMPARISONS,
                   log_scale=True, out_xlsx=OUT_XLSX, out_fig=OUT_FIG)
    print('\nDone.')


if __name__ == '__main__':
    main()
