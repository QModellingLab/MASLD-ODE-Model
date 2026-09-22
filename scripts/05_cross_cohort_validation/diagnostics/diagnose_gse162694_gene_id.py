#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
diagnose_gse162694_gene_id.py
================================================================
GSE162694 node mapping = 0/58，問題出在 gene id 格式假設錯誤。
本腳本印出 raw_counts.csv 實際 row index 前 20 筆，並對 node table
的 symbol / entrez_ids / ensembl_ids 三種 ID 分別計算重疊率，
找出真正對應的 ID 型態，不用猜。
================================================================
"""
import os
import pandas as pd

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))

COUNTS_CSV = os.path.join(SCRIPT_DIR, 'GSE162694_raw_counts.csv')
NODE_TABLE = os.path.join(SCRIPT_DIR, 'node_name_table_hsa04932.xlsx')


def main():
    counts = pd.read_csv(COUNTS_CSV, index_col=0, nrows=2000)  # 只讀前 2000 列加速
    idx = counts.index.astype(str)
    print(f'[raw_counts.csv] row index dtype 範例（前 20 筆）:')
    for v in idx[:20]:
        print(f'    {repr(v)}')
    print(f'\n  總列數（檔案全長，重新計算）: ', end='')
    full_idx = pd.read_csv(COUNTS_CSV, usecols=[0]).iloc[:, 0].astype(str)
    print(len(full_idx))

    df = pd.read_excel(NODE_TABLE, 'NodeNameTable')

    def collect_ids(col):
        ids = set()
        for v in df[col].dropna():
            for g in str(v).split(','):
                g = g.strip().rstrip('*').strip()
                if '(' in g:
                    g = g.split('(')[0].strip()
                if g and g.lower() != 'nan':
                    ids.add(g)
        return ids

    symbol_ids = collect_ids('gene_symbols')
    entrez_ids = collect_ids('entrez_ids')
    ensembl_ids = collect_ids('ensembl_ids')

    full_idx_set = set(full_idx)
    full_idx_set_nover = set(full_idx.str.split('.').str[0])

    print(f'\n[比對 node table 各 ID 型態 與 raw_counts.csv index 的重疊數]')
    print(f'  gene_symbols  : {len(symbol_ids)} 個唯一值, 重疊 = {len(symbol_ids & full_idx_set)}')
    print(f'  entrez_ids    : {len(entrez_ids)} 個唯一值, 重疊 = {len(entrez_ids & full_idx_set)}'
          f' (去版本號重疊 = {len(entrez_ids & full_idx_set_nover)})')
    print(f'  ensembl_ids   : {len(ensembl_ids)} 個唯一值, 重疊 = {len(ensembl_ids & full_idx_set)}'
          f' (去版本號重疊 = {len(ensembl_ids & full_idx_set_nover)})')

    print(f'\n  node table gene_symbols 範例: {list(symbol_ids)[:10]}')
    print(f'  node table entrez_ids 範例:   {list(entrez_ids)[:10]}')
    print(f'  node table ensembl_ids 範例:  {list(ensembl_ids)[:10]}')


if __name__ == '__main__':
    main()
