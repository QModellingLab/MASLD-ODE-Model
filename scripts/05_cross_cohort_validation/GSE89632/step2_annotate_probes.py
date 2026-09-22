#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
step2_annotate_probes.py  (GSE89632 version)
================================================================
Step 2 of GSE89632 independent validation pipeline.

Download GPL14951 (Illumina HumanHT-12 WG-DASL V4.0 R2) platform
annotation from NCBI GEO and build a probe_id -> gene_symbol mapping.

GPL14951 annotation table includes columns:
  ID            : Illumina probe ID (e.g., ILMN_xxxxxxx)
  Symbol        : HGNC gene symbol  <- primary
  ILMN_Gene     : Illumina gene name
  Entrez_Gene_ID: NCBI Entrez ID
  RefSeq_ID     : RefSeq transcript
  Synonyms      : alternative symbols

Output:
  - GPL14951_probe_to_gene.csv
================================================================
"""
import os
import pandas as pd

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
CACHE_DIR  = os.path.join(SCRIPT_DIR, 'geo_cache')
OUT_ANNOT  = os.path.join(SCRIPT_DIR, 'GPL14951_probe_to_gene.csv')


def main():
    try:
        import GEOparse
    except ImportError:
        raise SystemExit('GEOparse not installed.  Run: pip install GEOparse')

    os.makedirs(CACHE_DIR, exist_ok=True)
    print('Downloading GPL14951 from GEO (first run only)...')
    gpl = GEOparse.get_GEO(geo='GPL14951', destdir=CACHE_DIR, silent=False)

    df = gpl.table.copy()
    print(f'\nAnnotation rows: {len(df):,}')
    print(f'Columns sample: {list(df.columns)[:12]}')

    # Illumina GPL14951 typically uses 'Symbol' as gene symbol column.
    # Some legacy tables may use 'ILMN_Gene' or 'Gene Symbol'.
    SYMBOL_COL = None
    for candidate in ['Symbol', 'Gene Symbol', 'Gene_Symbol',
                       'ILMN_Gene', 'SYMBOL']:
        if candidate in df.columns:
            SYMBOL_COL = candidate
            break
    if SYMBOL_COL is None:
        print(f'Available columns: {list(df.columns)}')
        raise SystemExit(
            'Could not find gene symbol column in GPL14951 table.')
    print(f'Using gene symbol column: "{SYMBOL_COL}"')

    out = df[['ID', SYMBOL_COL]].rename(
        columns={'ID': 'probe_id', SYMBOL_COL: 'gene_symbol'})
    # Clean: drop empty, drop '---' placeholders
    out = out.dropna(subset=['gene_symbol'])
    out = out[out['gene_symbol'].astype(str).str.strip() != '']
    out = out[out['gene_symbol'] != '---']
    out['probe_id'] = out['probe_id'].astype(str)
    out['gene_symbol'] = out['gene_symbol'].astype(str).str.strip()
    out = out.drop_duplicates()

    print(f'\nProbes with gene symbol : {len(out):,} / {len(df):,}  '
          f'({len(out)/len(df)*100:.1f}%)')
    print(f'Unique gene symbols     : {out["gene_symbol"].nunique():,}')

    out.to_csv(OUT_ANNOT, index=False)
    print(f'\nSaved: {OUT_ANNOT}')


if __name__ == '__main__':
    main()
