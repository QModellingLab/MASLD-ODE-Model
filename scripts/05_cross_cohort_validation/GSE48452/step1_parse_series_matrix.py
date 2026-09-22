#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
step1_parse_series_matrix.py
================================================================
Step 1 of GSE48452 independent validation pipeline.

Parse a GEO series matrix file to extract:
  (a) sample metadata (group, bariatric surgery status, NAS, fibrosis, etc.)
  (b) probe-level expression matrix (RMA log2-normalised by submitter)

Outputs CSV next to this script:
  - GSE48452_sample_metadata.csv
  - GSE48452_expression_probe.csv

How to run (Spyder, deseq2_env)
  1. Download GSE48452_series_matrix.txt from
       https://ftp.ncbi.nlm.nih.gov/geo/series/GSE48nnn/GSE48452/matrix/
     and place in the same folder as this script.
  2. Open in Spyder, F5.
================================================================
"""
import os
import re
import pandas as pd
import numpy as np

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
INPUT  = os.path.join(SCRIPT_DIR, 'GSE48452_series_matrix.txt')
OUT_META = os.path.join(SCRIPT_DIR, 'GSE48452_sample_metadata.csv')
OUT_EXPR = os.path.join(SCRIPT_DIR, 'GSE48452_expression_probe.csv')


def split_quoted_tabs(line):
    parts = line.rstrip('\n').split('\t')
    return [p.strip('"') for p in parts[1:]]


def main():
    if not os.path.exists(INPUT):
        raise FileNotFoundError(
            f'Series matrix not found: {INPUT}\n'
            'Download from https://ftp.ncbi.nlm.nih.gov/geo/series/'
            'GSE48nnn/GSE48452/matrix/')

    # --- Locate data block start ---
    data_start = None
    with open(INPUT, encoding='utf-8') as f:
        for i, line in enumerate(f):
            if line.startswith('!series_matrix_table_begin'):
                data_start = i + 1
                break
    print(f'Data block starts at line {data_start + 1}')

    # --- Parse sample metadata (header section) ---
    meta_keys = {
        '!Sample_geo_accession':   'GSM',
        '!Sample_source_name_ch1': 'source_name',
    }
    char_keys = ['group:', 'bariatric surgery:', 'Sex:', 'age:', 'bmi:',
                  'nas:', 'fibrosis:', 'fat:', 'inflammation:', 'other_id:']

    meta = {v: None for v in meta_keys.values()}
    chars = {k.rstrip(':'): None for k in char_keys}

    with open(INPUT, encoding='utf-8') as f:
        for line in f:
            if line.startswith('!series_matrix_table_begin'):
                break
            for prefix, label in meta_keys.items():
                if line.startswith(prefix):
                    meta[label] = split_quoted_tabs(line)
            if line.startswith('!Sample_characteristics_ch1'):
                for ck in char_keys:
                    if f'"{ck}' in line:
                        vals = split_quoted_tabs(line)
                        chars[ck.rstrip(':')] = [
                            v.replace(ck + ' ', '') for v in vals]
                        break

    meta_df = pd.DataFrame({**meta, **chars})

    # Infer hybridisation batch (A1359 vs A1649) from other_id
    def parse_batch(x):
        if isinstance(x, str) and x not in ('NA', ''):
            m = re.match(r'(A\d+)-', x)
            return m.group(1) if m else None
        return None
    meta_df['batch'] = meta_df['other_id'].map(parse_batch)

    # Numeric coercion
    for c in ('age', 'bmi', 'nas', 'fibrosis', 'fat', 'inflammation'):
        meta_df[c] = pd.to_numeric(meta_df[c], errors='coerce')

    # 'usable' flag = exclude post-bariatric-surgery samples
    meta_df['usable'] = meta_df['bariatric surgery'] != 'after surgery'

    print('\nSample metadata parsed:')
    print(meta_df.groupby(['group', 'bariatric surgery']).size()
                   .unstack(fill_value=0))
    print(f'\nUsable (excl. after surgery): '
          f'{meta_df.usable.sum()}/{len(meta_df)}')
    meta_df.to_csv(OUT_META, index=False)
    print(f'Saved: {OUT_META}')

    # --- Parse expression matrix ---
    print('\nReading expression matrix...')
    expr = pd.read_csv(INPUT, sep='\t', skiprows=data_start,
                        comment='!', index_col=0)
    expr.index.name = 'probe_id'
    print(f'Expression matrix: {expr.shape[0]:,} probes × {expr.shape[1]} samples')
    print(f'Value range: {expr.values.min():.2f}-{expr.values.max():.2f} '
          f'(log2 RMA; expected ~1-14)')
    expr.to_csv(OUT_EXPR)
    print(f'Saved: {OUT_EXPR}')


if __name__ == '__main__':
    main()
