#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
diagnose_gse48452_89632_sample_n.py
================================================================
查證 GSE48452 / GSE89632 的實際樣本分組數，直接從 GEO series_matrix.txt
的 !Sample_characteristics_ch1 解析，不依賴文獻引用的總樣本數。

請把 GSE48452_series_matrix.txt 和/或 GSE89632_series_matrix.txt
放進本腳本所在資料夾（跟你先前驗證 GSE130970 等三個資料集時做法相同）。
若檔案不存在，腳本會列出 MISS 並跳過該資料集，不會中斷。

跑法：
  %run diagnose_gse48452_89632_sample_n.py
================================================================
"""
import os
import pandas as pd

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))


def parse_series_matrix(path):
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


def diagnose(dataset):
    path = os.path.join(SCRIPT_DIR, f'{dataset}_series_matrix.txt')
    print('=' * 65)
    print(f'[{dataset}]')
    if not os.path.exists(path):
        print(f'  MISS {path}')
        print('  -> 請把這個檔案放進本資料夾後重跑')
        return
    print(f'  OK   {path}')

    meta = parse_series_matrix(path)
    char_cols = [c for c in meta.columns if c.startswith('char_')]
    print(f'  {len(meta)} samples total, {len(char_cols)} characteristics 欄位')

    for c in char_cols:
        print(f'\n  -- {c} 唯一值分布 --')
        print(meta[c].value_counts(dropna=False).to_string())

    print('\n  -- 前 10 筆 title（供核對命名規則）--')
    for tt in meta['title'].tolist()[:10]:
        print(f'    {tt}')

    out_csv = os.path.join(SCRIPT_DIR, f'{dataset}_characteristics_preview.csv')
    meta.to_csv(out_csv, index=False)
    print(f'\n  完整 characteristics 已存成: {out_csv}')


def main():
    for ds in ['GSE48452', 'GSE89632']:
        diagnose(ds)
        print()


if __name__ == '__main__':
    main()
