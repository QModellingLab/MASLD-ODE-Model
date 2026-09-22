#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
diagnose_sample_id_matching.py
================================================================
用途：對 GSE162694 / GSE213621，窮舉常見的 GSM<->矩陣欄名 轉換規則，
      自動算出每種規則的 matched 樣本數，印出結果排序，
      不用人工肉眼猜 pattern。

跑法：
  %run diagnose_sample_id_matching.py GSE162694
  %run diagnose_sample_id_matching.py GSE213621
================================================================
"""
import os
import re
import sys
import pandas as pd

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))


def parse_series_matrix(path):
    gsm_line, title_line = None, None
    with open(path, encoding='utf-8', errors='ignore') as fh:
        for line in fh:
            if line.startswith('!Sample_geo_accession'):
                gsm_line = line
            elif line.startswith('!Sample_title'):
                title_line = line
            if gsm_line and title_line:
                break

    def split_row(line):
        return [v.strip('"') for v in line.strip().split('\t')[1:]]

    gsms = split_row(gsm_line)
    titles = split_row(title_line) if title_line else [''] * len(gsms)
    return pd.DataFrame({'GSM': gsms, 'title': titles})


def load_matrix_columns(dataset):
    if dataset == 'GSE162694':
        path = os.path.join(SCRIPT_DIR, 'GSE162694_raw_counts.csv')
        cols = pd.read_csv(path, index_col=0, nrows=0).columns.tolist()
    elif dataset == 'GSE213621':
        path = os.path.join(SCRIPT_DIR, 'GSE213621_FPKMs_allsamples.txt')
        cols = pd.read_csv(path, sep=None, engine='python', index_col=0, nrows=0).columns.tolist()
    else:
        raise ValueError(dataset)
    return cols


def candidate_transforms(title_series):
    """回傳 {規則名稱: 轉換後的 Series}"""
    cands = {}
    cands['title (原樣)'] = title_series
    cands['title 第一段(空白前)'] = title_series.str.split().str[0]
    cands['title 最後一段(空白後)'] = title_series.str.split().str[-1]
    cands['title 去空白'] = title_series.str.replace(r'\s+', '', regex=True)
    cands['Sample_+title最後一段'] = 'Sample_' + title_series.str.split().str[-1]
    cands['title最後一段數字'] = title_series.str.extract(r'(\d+)\s*$')[0]
    cands['Sample_+title最後數字'] = 'Sample_' + title_series.str.extract(r'(\d+)\s*$')[0].astype(str)
    cands['title第一段去底線後綴'] = title_series.str.split().str[0].str.split('_').str[0]
    return cands


def main():
    if len(sys.argv) < 2 or sys.argv[1] not in ('GSE162694', 'GSE213621'):
        print('用法: %run diagnose_sample_id_matching.py GSE162694 (或 GSE213621)')
        return
    dataset = sys.argv[1]

    sm_path = os.path.join(SCRIPT_DIR, f'{dataset}_series_matrix.txt')
    meta = parse_series_matrix(sm_path)
    cols = load_matrix_columns(dataset)
    cols_set = set(cols)

    print(f'[{dataset}] series_matrix: {len(meta)} samples | 矩陣欄位: {len(cols)} columns')
    print(f'  矩陣欄名前 5 個: {cols[:5]}')
    print(f'  GSM 前 5 個: {meta["GSM"].tolist()[:5]}')
    print(f'  title 前 5 個: {meta["title"].tolist()[:5]}')

    print('\n[窮舉比對規則，依 matched 數排序]')
    results = []
    results.append(('GSM (原樣)', meta['GSM'].isin(cols_set).sum()))
    for name, series in candidate_transforms(meta['title']).items():
        n = series.isin(cols_set).sum()
        results.append((name, n))

    results.sort(key=lambda x: -x[1])
    for name, n in results:
        print(f'  {name:<28} matched = {n}/{len(meta)}')

    best_name, best_n = results[0]
    print(f'\n[結論] 最佳規則: "{best_name}"，matched {best_n}/{len(meta)}')
    if best_n < len(meta):
        print(f'  仍有 {len(meta)-best_n} 筆未對上，列出前 10 筆對不上的 title 供人工核對：')
        cands = candidate_transforms(meta['title'])
        cands['GSM (原樣)'] = meta['GSM']
        best_series = cands[best_name] if best_name in cands else meta['GSM']
        unmatched = meta[~best_series.isin(cols_set)]
        for t in unmatched['title'].tolist()[:10]:
            print(f'    {t}')


if __name__ == '__main__':
    main()
