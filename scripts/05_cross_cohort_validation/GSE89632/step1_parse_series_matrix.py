#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
step1_parse_series_matrix.py  (GSE89632 version)
================================================================
Step 1 of GSE89632 independent validation pipeline (second cohort
for cross-cohort robustness of the MASLD ODE model).

Parse GSE89632 series matrix to extract:
  (a) sample metadata (group, sex, age, BMI, fibrosis, NAS, ...)
  (b) probe-level expression matrix
      (Illumina HumanHT-12 WG-DASL V4.0 R2, GPL14951)

GSE89632 sample design (after QC by Arendt et al. 2015):
  HC   (Healthy controls)         n=24
  SS   (Simple Steatosis / NAFL)  n=20
  NASH (Steatohepatitis)          n=19
  Total: 63 samples

Outputs (saved next to this script):
  - GSE89632_sample_metadata.csv
  - GSE89632_expression_probe.csv

How to run (Anaconda Spyder, Windows, deseq2_env)
  1. Download GSE89632_series_matrix.txt from
       https://ftp.ncbi.nlm.nih.gov/geo/series/GSE89nnn/GSE89632/matrix/
     and place in this folder.
  2. F5 in Spyder.
================================================================
"""
import os
import re
import pandas as pd
import numpy as np

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
INPUT    = os.path.join(SCRIPT_DIR, 'GSE89632_series_matrix.txt')
OUT_META = os.path.join(SCRIPT_DIR, 'GSE89632_sample_metadata.csv')
OUT_EXPR = os.path.join(SCRIPT_DIR, 'GSE89632_expression_probe.csv')


def split_quoted_tabs(line):
    parts = line.rstrip('\n').split('\t')
    return [p.strip('"') for p in parts[1:]]


def main():
    if not os.path.exists(INPUT):
        raise FileNotFoundError(
            f'Series matrix not found: {INPUT}\n'
            'Download from https://ftp.ncbi.nlm.nih.gov/geo/series/'
            'GSE89nnn/GSE89632/matrix/')

    # Locate data block
    data_start = None
    with open(INPUT, encoding='utf-8') as f:
        for i, line in enumerate(f):
            if line.startswith('!series_matrix_table_begin'):
                data_start = i + 1
                break
    print(f'Data block starts at line {data_start + 1}')

    # --- Parse sample metadata ---
    # Header fields we always want
    meta_keys = {
        '!Sample_geo_accession':   'GSM',
        '!Sample_title':           'sample_title',
        '!Sample_source_name_ch1': 'source_name',
    }
    meta = {v: None for v in meta_keys.values()}

    # All characteristics fields → store as a dict keyed by prefix
    # GSE89632 typically has fields like "diagnosis:", "Sex:", "age:",
    # "bmi:", "fibrosis:", "steatosis:", "lobular inflammation:", etc.
    # We parse all and let user filter downstream.
    char_data = {}

    with open(INPUT, encoding='utf-8') as f:
        for line in f:
            if line.startswith('!series_matrix_table_begin'):
                break
            for prefix, label in meta_keys.items():
                if line.startswith(prefix):
                    meta[label] = split_quoted_tabs(line)
            if line.startswith('!Sample_characteristics_ch1'):
                vals = split_quoted_tabs(line)
                # Determine field name from first non-empty value
                key = None
                for v in vals:
                    if ':' in v:
                        key = v.split(':', 1)[0].strip()
                        break
                if key is None:
                    continue
                stripped_vals = [
                    v.split(':', 1)[1].strip() if ':' in v else v
                    for v in vals]
                # Make column name safe (replace spaces)
                col = key.lower().replace(' ', '_').replace('-', '_')
                # If duplicated key, suffix _2, _3, ...
                base = col; n = 2
                while col in char_data:
                    col = f'{base}_{n}'; n += 1
                char_data[col] = stripped_vals

    meta_df = pd.DataFrame({**meta, **char_data})

    # Try to identify the "disease group" column heuristically
    # Common labels: 'diagnosis', 'group', 'phenotype', 'disease_state',
    # 'condition', 'sample_type'
    GROUP_CANDIDATES = ['diagnosis', 'group', 'phenotype', 'disease_state',
                         'condition', 'sample_type', 'disease',
                         'tissue_group', 'patient_group']
    group_col = None
    for cand in GROUP_CANDIDATES:
        if cand in meta_df.columns:
            # Check if values look like HC/SS/NASH or similar
            unique_vals = meta_df[cand].dropna().unique()
            if 2 <= len(unique_vals) <= 6:
                group_col = cand
                break
    # Fallback: if no group_col found, try source_name parsing
    if group_col is None and meta_df['source_name'].notna().all():
        # Heuristic: source_name may contain "HC", "SS", "NASH"
        meta_df['source_group'] = meta_df['source_name'].apply(
            lambda s: ('HC'   if re.search(r'\bHC\b|healthy', s, re.I) else
                       'NASH' if re.search(r'\bNASH\b', s, re.I)        else
                       'SS'   if re.search(r'\bSS\b|steatosis', s, re.I) else
                       None))
        group_col = 'source_group'

    if group_col is None:
        print('WARNING: Could not auto-detect disease group column.')
        print(f'Available characteristic columns: {list(char_data.keys())}')
        print('Edit step3_DEG.py manually to specify GROUP_COL.')
    else:
        print(f'Detected disease group column: "{group_col}"')
        meta_df['group'] = meta_df[group_col]
        # Standardise group labels
        STD_MAP = {
            'hc': 'HC', 'healthy': 'HC', 'healthy control': 'HC',
            'normal': 'HC', 'control': 'HC',
            'ss': 'SS', 'simple steatosis': 'SS', 'steatosis': 'SS',
            'nafl': 'SS',
            'nash': 'NASH', 'steatohepatitis': 'NASH',
            'nonalcoholic steatohepatitis': 'NASH',
        }
        meta_df['group'] = meta_df['group'].astype(str).str.strip().str.lower().map(
            lambda x: STD_MAP.get(x, x))

    # Numeric coercion for any numeric-like fields
    for c in meta_df.columns:
        if c in ('GSM', 'sample_title', 'source_name', 'group',
                 group_col or ''):
            continue
        # Try numeric
        try:
            numeric = pd.to_numeric(meta_df[c], errors='coerce')
            if numeric.notna().sum() > len(meta_df) * 0.3:
                meta_df[c] = numeric
        except Exception:
            pass

    print('\nSample metadata parsed:')
    if 'group' in meta_df.columns:
        print(meta_df['group'].value_counts().to_string())

    meta_df.to_csv(OUT_META, index=False)
    print(f'\nSaved: {OUT_META}')

    # --- Parse expression matrix ---
    print('\nReading expression matrix...')
    expr = pd.read_csv(INPUT, sep='\t', skiprows=data_start,
                        comment='!', index_col=0)
    expr.index.name = 'probe_id'
    print(f'Expression matrix: {expr.shape[0]:,} probes × '
          f'{expr.shape[1]} samples')
    print(f'Value range: {expr.values.min():.2f}-{expr.values.max():.2f}')
    print(f'Probe ID format examples: {list(expr.index[:5])}')

    expr.to_csv(OUT_EXPR)
    print(f'Saved: {OUT_EXPR}')


if __name__ == '__main__':
    main()
