#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
compute_node_ratios_GSE130970.py
================================================================
GSE130970（Hoang et al. 2019, Sci Rep）— salmon/tximport TPM,
entrez gene ID 為 row index + GSE130970_series_matrix.txt

流程（與 GSE213621 腳本邏輯一致）：
  1. 解析 series_matrix.txt 取得每個 GSM 的 title + 所有
     !Sample_characteristics_ch1 欄位
  2. --dry-run 印出每個 characteristics 欄位的唯一值分布，供核對
     （Hoang 2019 通常會把 NAS score / fibrosis stage / diagnosis
      放在 characteristics 裡）
  3. 依 GROUP_RULES 規則（在 characteristics 文字中比對關鍵字）
     把樣本分類成 Normal / NAFL / NASH
  4. TPM -> log2(TPM+1) -> 各組平均 -> node ratio（用 entrez_ids 對應）
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
DATASET = 'GSE130970'

TPM_CSV       = os.path.join(SCRIPT_DIR, 'GSE130970_all_sample_salmon_tximport_TPM_entrez_gene_ID.csv')
SERIES_MATRIX = os.path.join(SCRIPT_DIR, 'GSE130970_series_matrix.txt')
NODE_TABLE    = os.path.join(SCRIPT_DIR, 'node_name_table_hsa04932.xlsx')

REFERENCE_GROUP = 'Normal'
COMPARISONS = {'NAFL': 'NAFL', 'NASH': 'NASH'}

# 分類規則：依文獻查證鎖定（非猜測）。
# 來源：Dynamic co-expression modular network analysis in NAFLD
#       (PMC8380347)，該文分析同一個 GSE130970 dataset，明確記載：
#   "According to the NIDDK NASH CRN criteria, 26 samples with NAS >= 5
#    were regarded as NASH, 10 samples with steatosis and NAS < 3 were
#    considered as NAFL, 4 samples with NAS = 0 was chosen as normal
#    tissues. ... the rest 36 samples could not be defined and were not
#    enrolled."
# 與本資料夾 --dry-run 印出的 NAS 分布完全吻合（NAS=0: 4例；
# NAS=5+6: 18+8=26例），故採用此標準：
#   Normal : NAS == 0
#   NAFL   : steatosis >= 1 且 NAS < 3
#   NASH   : NAS >= 5
#   其餘（NAS 3-4，或 NAS 1-2 但 steatosis=0）：排除（不分類）
BALLOON_COL = 'char_5'     # cytological ballooning grade（本版未使用，保留供參考）
STEATOSIS_COL = 'char_6'   # steatosis grade
NAS_COL = 'char_7'         # nafld activity score


def parse_score(text):
    """從 'xxx grade: N' 或 'xxx score: N' 取出整數 N。"""
    m = re.search(r':\s*(\d+)', str(text))
    return int(m.group(1)) if m else None

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
    steatosis = parse_score(row.get(STEATOSIS_COL))
    nas = parse_score(row.get(NAS_COL))
    if steatosis is None or nas is None:
        return None
    if nas == 0:
        return 'Normal'
    if nas >= 5:
        return 'NASH'
    if steatosis >= 1 and nas < 3:
        return 'NAFL'
    return None   # NAS 3-4 borderline，或 NAS 1-2 但無 steatosis -> 排除


def main():
    dry_run = '--dry-run' in sys.argv

    print('=' * 65)
    print(f'compute_node_ratios_GSE130970.py  [{DATASET}]'
          + ('  (DRY RUN)' if dry_run else ''))
    print('=' * 65)

    for f, lb in [(TPM_CSV, 'TPM_CSV'), (SERIES_MATRIX, 'SERIES_MATRIX'),
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

    print('\n[Node gene map] (entrez_id)')
    node_gene_map = load_node_gene_map(NODE_TABLE, id_type='entrez')
    n_multi = sum(1 for v in node_gene_map.values() if len(v) > 1)
    print(f'  {len(node_gene_map)} nodes, {n_multi} multi-gene nodes')

    print('\n[Loading TPM matrix]')
    tpm = pd.read_csv(TPM_CSV, index_col=0)
    tpm.index = tpm.index.astype(str).str.split('.').str[0]   # entrez id 去除小數
    print(f'  {tpm.shape[0]:,} genes x {tpm.shape[1]} samples')

    matched_meta = meta[meta['group'].notna() & meta['GSM'].isin(tpm.columns)]
    print(f'  {len(matched_meta)}/{len(meta)} samples matched to TPM 欄名')
    if len(matched_meta) == 0:
        matched_meta = meta[meta['group'].notna() & meta['title'].isin(tpm.columns)]
        id_col = 'title'
        print(f'  改用 title 比對: {len(matched_meta)}/{len(meta)} samples matched')
        if len(matched_meta) == 0:
            raise ValueError('GSM 與 title 都對不上 TPM 矩陣欄名，請手動檢查欄名格式')
    else:
        id_col = 'GSM'

    log2_tpm = np.log2(tpm[matched_meta[id_col].tolist()] + 1.0)
    log2_tpm.columns = matched_meta[id_col].tolist()

    all_groups = [REFERENCE_GROUP] + list(COMPARISONS.values())
    gm = {}
    for grp in all_groups:
        ids = matched_meta.loc[matched_meta['group'] == grp, id_col].tolist()
        if not ids:
            print(f'  [WARNING] Group "{grp}": 0 samples')
            continue
        gm[grp] = log2_tpm[ids].mean(axis=1)
        print(f'  {grp}: {len(ids)} samples')
    gm = pd.DataFrame(gm)

    print('\n[Computing node ratios]')
    print_and_save(DATASET, node_gene_map, gm, REFERENCE_GROUP, COMPARISONS,
                   log_scale=True, out_xlsx=OUT_XLSX, out_fig=OUT_FIG)
    print('\nDone.')


if __name__ == '__main__':
    main()
