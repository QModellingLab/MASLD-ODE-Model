#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
compare_node_ratios_vs_Ahrens2013_TableS2_v2.py
================================================================
v2 新增：從 GSE48452_DEG_gene_level.xlsx 補上
         not-in-model 基因的 log2FC，讓全部 59 個 Table S2 基因
         都有對應的「我們算出的 DEG 數值」可以比較。

輸入：
  node_ratios_GSE48452.xlsx          (由 compute_node_ratios_from_expression.py 產出)
  GSE48452_DEG_gene_level.xlsx       (由 step3_DEG.py 產出)

輸出：
  Table_GSE48452_NodeRatios_vs_Ahrens2013_TableS2_v2.xlsx
================================================================
"""
import os
import numpy as np
import pandas as pd

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
RATIOS_XLSX = os.path.join(SCRIPT_DIR, 'node_ratios_GSE48452.xlsx')
DEG_XLSX    = os.path.join(SCRIPT_DIR, 'GSE48452_DEG_gene_level.xlsx')
OUT_XLSX    = os.path.join(SCRIPT_DIR,
    'Table_GSE48452_NodeRatios_vs_Ahrens2013_TableS2_v2.xlsx')

# ----------------------------------------------------------------
# Ahrens 2013 Table S2 EXPRESSION 欄位（從圖片讀出）
# D.MED.NC = NASH - Control median log2 expression difference
# ----------------------------------------------------------------
PAPER_DATA = [
    ('ZNF274',   -0.27, 8.0E-03, 4.4E-02),
    ('PPARGC1A', -0.51, 1.9E-02, 1.7E-01),
    ('SREBF2',    0.45, 1.6E-02, 2.8E-02),
    ('ESRRA',    -0.08, 5.3E-01, 8.1E-01),
    ('GRP20',    -0.12, 7.4E-02, 9.4E-03),
    ('HSF1',     -0.12, 1.4E-01, 3.3E-01),
    ('NFE2L2',   -0.18, 1.8E-01, 1.7E-02),
    ('CEBPZ',    -0.19, 6.9E-03, 6.0E-03),
    ('SREBF1',    0.45, 1.0E-01, 1.0E-01),
    ('ZEB1',     -0.08, 6.2E-01, 6.6E-03),
    ('BHLHE40',  -0.20, 5.6E-01, 8.1E-01),
    ('JUN',      -1.68, 9.2E-05, 8.4E-01),
    ('TCF12',     0.01, 7.7E-01, 2.2E-02),
    ('NR2C2',     0.09, 1.8E-01, 3.1E-01),
    ('HNF4A',     0.02, 7.1E-01, 6.7E-01),
    ('ZBTB33',    0.30, 3.2E-01, 3.1E-03),
    ('NRF1',      0.01, 8.8E-01, 8.7E-01),
    ('RXRA',      0.08, 5.1E-01, 1.1E-01),
    ('SRF',      -0.25, 1.3E-02, 3.7E-01),
    ('CEBPB',    -0.32, 4.8E-02, 2.7E-01),
    ('SP2',      -0.06, 3.0E-01, 2.7E-01),
    ('USF2',     -0.01, 5.0E-01, 6.6E-01),
    ('HNF4G',    -0.11, 9.8E-01, 5.5E-02),
    ('FOXA1',     0.11, 2.4E-01, 2.2E-01),
    ('FOXA2',     0.11, 3.0E-01, 5.3E-01),
    ('EP300',     0.02, 7.1E-01, 4.8E-01),
    ('SP1',       0.01, 8.4E-01, 1.0E+00),
    ('FOSL2',     0.00, 1.0E+00, 8.0E-01),
    ('BRCA1',     0.26, 1.3E-01, 1.0E-01),
    ('CREB1',     0.06, 2.5E-01, 1.5E-03),
    ('MYC',      -0.61, 1.8E-01, 1.2E-01),
    ('PRICKLE1',  0.45, 2.7E-01, 1.6E-03),
    ('JUND',     -0.13, 4.3E-04, 7.0E-01),
    ('RFX5',      0.30, 9.2E-03, 1.2E-02),
    ('CEBPD',    -0.48, 4.2E-03, 9.9E-02),
    ('TAF1',     -0.09, 8.2E-01, 5.5E-01),
    ('SAP18',    -0.11, 2.9E-02, 6.0E-02),
    ('NFIC',     -0.06, 6.3E-01, 7.7E-01),
    ('YY1',      -0.03, 1.3E-01, 9.6E-01),
    ('MXI1',     -0.14, 1.8E-01, 9.7E-02),
    ('RCOR1',     0.03, 6.3E-01, 3.5E-02),
    ('TBP',      -0.09, 9.5E-01, 3.1E-01),
    ('ELF1',     -0.13, 2.3E-02, 2.5E-01),
    ('HDAC2',    -0.03, 7.6E-04, 5.2E-01),
    ('MBD4',     -0.11, 2.1E-01, 8.2E-01),
    ('MAZ',       0.01, 9.8E-01, 5.7E-01),
    ('MYBL2',     0.04, 4.0E-01, 3.6E-01),
    ('RAD21',     0.05, 5.3E-01, 3.5E-01),
    ('TEAD4',     0.10, 5.6E-01, 1.2E-01),
    ('GABPA',    -0.12, 7.4E-02, 3.4E-03),
    ('ZBTB7A',   -0.05, 1.9E-01, 2.4E-01),
    ('CHD2',      0.03, 9.4E-01, 4.1E-02),
    ('ARID3A',    0.04, 3.0E-01, 4.5E-01),
    ('USF1',      0.07, 2.1E-01, 1.4E-01),
    ('SMC3',      0.22, 3.0E-01, 3.7E-01),
    ('TCF7L2',   -0.35, 5.8E-04, 8.9E-02),
    ('ATF3',     -0.84, 6.1E-05, 6.3E+00),
    ('IRF3',     -0.18, 9.3E-02, 2.0E+00),
    ('CTCF',      0.12, 6.5E-02, 8.2E-01),
]

def dir_sym(val, up=0.05, dn=-0.05):
    if val > up:  return 'up'
    if val < dn:  return 'dn'
    return 'nc'

def main():
    # ---- Load node ratios ----
    if not os.path.exists(RATIOS_XLSX):
        raise FileNotFoundError(f'Not found: {RATIOS_XLSX}')
    df_nodes = pd.read_excel(RATIOS_XLSX, sheet_name='AllNodes')

    # Build gene → node+ratio lookup
    gene_to_node = {}
    for _, row in df_nodes.iterrows():
        genes = [g.strip() for g in str(row['genes_found']).split(',')
                 if g.strip() not in ('', 'nan', 'NaN')]
        for g in genes:
            gene_to_node[g.upper()] = {
                'node':            row['node'],
                'ratio_NASH':      row.get('ratio_NASH', np.nan),
                'log2_NASH':       row.get('log2_NASH', np.nan),
            }

    # ---- Load DEG log2FC for all genes ----
    deg_lookup = {}
    if os.path.exists(DEG_XLSX):
        xls = pd.ExcelFile(DEG_XLSX)
        # try NASH_vs_Normal sheet
        sheet = 'NASH_vs_Normal' if 'NASH_vs_Normal' in xls.sheet_names else xls.sheet_names[0]
        df_deg = pd.read_excel(xls, sheet)
        # detect gene and log2FC columns
        gene_col = next((c for c in df_deg.columns
                         if 'gene' in c.lower()), df_deg.columns[0])
        lfc_col  = next((c for c in df_deg.columns
                         if 'log2' in c.lower() or 'fc' in c.lower()), df_deg.columns[1])
        padj_col = next((c for c in df_deg.columns
                         if 'padj' in c.lower() or 'adj' in c.lower()), None)
        for _, row in df_deg.iterrows():
            g = str(row[gene_col]).strip().upper()
            lfc  = row[lfc_col]
            padj = row[padj_col] if padj_col else np.nan
            deg_lookup[g] = {'log2FC': lfc, 'padj': padj}
        print(f'DEG loaded: {len(deg_lookup):,} genes from sheet "{sheet}"')
    else:
        print(f'[WARN] DEG file not found: {DEG_XLSX}')
        print('  Only node ratios will be reported for in-model genes.')

    # ---- Build comparison table ----
    rows = []
    for gene, dmed_nc, p_nc, pgkw in PAPER_DATA:
        paper_dir = dir_sym(dmed_nc)
        sig = '★' if p_nc < 0.05 else ''

        # ODE node ratio
        node_info = gene_to_node.get(gene.upper())
        if node_info:
            in_model  = 'Yes'
            node      = node_info['node']
            node_ratio = node_info['ratio_NASH']
            node_l2    = node_info['log2_NASH']
            node_dir   = dir_sym(node_l2) if np.isfinite(node_l2) else 'N/A'
        else:
            in_model   = 'No'
            node       = '—'
            node_ratio = np.nan
            node_l2    = np.nan
            node_dir   = 'N/A'

        # DEG log2FC (our calculation from raw expression)
        deg_info = deg_lookup.get(gene.upper())
        if deg_info:
            our_l2fc = deg_info['log2FC']
            our_padj = deg_info['padj']
            our_ratio = 2 ** our_l2fc if np.isfinite(our_l2fc) else np.nan
            our_dir  = dir_sym(our_l2fc)
            deg_match = '✅' if our_dir == paper_dir else '❌'
        else:
            our_l2fc  = np.nan
            our_padj  = np.nan
            our_ratio = np.nan
            our_dir   = 'N/A'
            deg_match = 'N/A'

        # Direction match using best available source
        if in_model == 'Yes':
            final_match = '✅' if node_dir == paper_dir else '❌'
            source = 'node_ratio'
        elif deg_info:
            final_match = deg_match
            source = 'DEG_log2FC'
        else:
            final_match = 'N/A'
            source = '—'

        rows.append({
            'Gene':             gene,
            'In_ODE_model':     in_model,
            'ODE_node':         node,
            'Paper_D.MED.NC':   dmed_nc,
            'Paper_P.NC':       p_nc,
            'Paper_sig':        sig,
            'Paper_dir':        paper_dir,
            # ODE node ratio
            'Node_ratio_NASH':  round(node_ratio, 3) if np.isfinite(node_ratio) else np.nan,
            'Node_log2_NASH':   round(node_l2, 3)    if np.isfinite(node_l2)    else np.nan,
            'Node_dir':         node_dir,
            # Our DEG
            'Our_log2FC_NASH':  round(our_l2fc, 3)   if np.isfinite(our_l2fc)  else np.nan,
            'Our_ratio_NASH':   round(our_ratio, 3)  if np.isfinite(our_ratio) else np.nan,
            'Our_padj':         our_padj,
            'Our_dir':          our_dir,
            # Final
            'Best_source':      source,
            'Direction_match':  final_match,
        })

    result = pd.DataFrame(rows)

    # ---- Console output ----
    print()
    print('=' * 95)
    print('GSE48452 vs Ahrens 2013 Table S2 — Full comparison (v2)')
    print('D.MED.NC = NASH-Control median log2 diff  |  ★ = P.NC < 0.05')
    print('Our_log2FC = from step3 DEG (Welch t-test)  |  Node_ratio = ODE model initial value')
    print('=' * 95)
    print(f'{"Gene":<12} {"ODE_node":<12} {"D.MED.NC":>10} {"P.NC":>10} {"★":>2} '
          f'{"P_dir":>6} {"Our_L2FC":>9} {"Our_dir":>8} {"Match":>7} {"Source"}')
    print('-' * 95)

    for _, r in result.iterrows():
        our_val = (f'{r["Our_log2FC_NASH"]:>+9.3f}'
                   if np.isfinite(r['Our_log2FC_NASH']) else '      N/A')
        print(f'{r.Gene:<12} {r.ODE_node:<12} '
              f'{r["Paper_D.MED.NC"]:>+10.3f} {r["Paper_P.NC"]:>10.2E} '
              f'{r.Paper_sig:>2} {r.Paper_dir:>6} '
              f'{our_val} {r.Our_dir:>8} '
              f'{r.Direction_match:>7}  {r.Best_source}')

    # Summary stats
    has_match = result[result['Direction_match'].isin(['✅','❌'])]
    n_match   = (has_match['Direction_match'] == '✅').sum()
    n_total   = len(has_match)
    sig_genes = result[result['Paper_sig'] == '★']
    sig_match = sig_genes[sig_genes['Direction_match'].isin(['✅','❌'])]

    print()
    print('=' * 60)
    print(f'Overall direction concordance:  {n_match}/{n_total} '
          f'({n_match/n_total*100:.0f}%)')
    sig_n_m = (sig_match['Direction_match'] == '✅').sum()
    print(f'Significant genes (P.NC<0.05):  {sig_n_m}/{len(sig_match)} '
          f'concordant ({sig_n_m/len(sig_match)*100:.0f}% if any)')
    print(f'Not found in our DEG:           '
          f'{result["Our_dir"].eq("N/A").sum()} genes')

    # ---- Save Excel ----
    with pd.ExcelWriter(OUT_XLSX, engine='openpyxl') as writer:
        result.to_excel(writer, sheet_name='Full_comparison', index=False)
        result[result['In_ODE_model']=='Yes'].to_excel(
            writer, sheet_name='InModel_genes', index=False)
        sig_genes.to_excel(
            writer, sheet_name='Significant_P05', index=False)
        # mismatch only
        mismatch = result[result['Direction_match'] == '❌']
        mismatch.to_excel(writer, sheet_name='Mismatch', index=False)

    print(f'\nExcel saved: {OUT_XLSX}')
    print('Done.')


if __name__ == '__main__':
    main()
