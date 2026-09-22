#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
step2_annotate_probes.py
================================================================
Step 2 of GSE48452 independent validation pipeline.

Download GPL11532 (Affymetrix HuGene 1.1 ST) platform annotation from NCBI GEO
and build a probe_id -> gene_symbol mapping table.

Requires:
  pip install GEOparse

First run downloads ~5 MB SOFT file from
  https://ftp.ncbi.nlm.nih.gov/geo/platforms/GPL11nnn/GPL11532/
Subsequent runs use the cached file in `geo_cache/`.

Output:
  - GPL11532_probe_to_gene.csv

How to run
  1. Ensure internet access is available.
  2. Open in Spyder, F5.
================================================================
"""
import os
import pandas as pd

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
CACHE_DIR  = os.path.join(SCRIPT_DIR, 'geo_cache')
OUT_ANNOT  = os.path.join(SCRIPT_DIR, 'GPL11532_probe_to_gene.csv')


def extract_symbol(s):
    """Extract gene symbol from Affymetrix gene_assignment field.

    Format example:
      "NM_001011892 // SCN9A // sodium channel ... /// NM_002977 // SCN9A // ..."
    We take the first " /// "-separated block, then the 2nd " // "-separated
    sub-field as the gene symbol.
    """
    if not isinstance(s, str) or s in ('---', ''):
        return None
    primary = s.split(' /// ')[0]
    parts = primary.split(' // ')
    if len(parts) >= 2:
        sym = parts[1].strip()
        return sym if sym and sym != '---' else None
    return None


def main():
    try:
        import GEOparse
    except ImportError:
        raise SystemExit(
            'GEOparse not installed.  Run: pip install GEOparse')

    os.makedirs(CACHE_DIR, exist_ok=True)
    print('Downloading GPL11532 from GEO (first run only)...')
    gpl = GEOparse.get_GEO(geo='GPL11532', destdir=CACHE_DIR, silent=False)

    df = gpl.table.copy()
    print(f'\nAnnotation rows: {len(df):,}')
    print(f'Columns sample: {list(df.columns)[:8]}')

    if 'gene_assignment' not in df.columns:
        raise SystemExit(
            "GPL11532 table missing 'gene_assignment' column. "
            f"Available: {list(df.columns)}")

    df['gene_symbol'] = df['gene_assignment'].map(extract_symbol)

    out = df[['ID', 'gene_symbol']].rename(columns={'ID': 'probe_id'})
    out = out.dropna(subset=['gene_symbol'])
    out['probe_id'] = pd.to_numeric(out['probe_id'], errors='coerce').dropna().astype(int)
    out = out.drop_duplicates()

    print(f'\nProbes with gene symbol : {len(out):,} / {len(df):,}  '
          f'({len(out)/len(df)*100:.1f}%)')
    print(f'Unique gene symbols     : {out["gene_symbol"].nunique():,}')

    out.to_csv(OUT_ANNOT, index=False)
    print(f'\nSaved: {OUT_ANNOT}')


if __name__ == '__main__':
    main()
