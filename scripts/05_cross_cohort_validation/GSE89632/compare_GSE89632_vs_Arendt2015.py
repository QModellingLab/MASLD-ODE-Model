#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
compare_GSE89632_vs_Arendt2015.py
================================================================
驗證 GSE89632 分析結果與 Arendt 2015 (Hepatology 61:1565) 之一致性

分析目標：
  1. Node initial ratios 是否與 paper 報告方向一致
  2. 基因層級 DEG 是否與 paper Table 2/3 一致

輸入：
  - GSE89632_node_initial_ratios.xlsx    (step5 產出，NASH/HC ratio)
  - GSE89632_DEG_gene_level.xlsx         (step3 產出，Welch t-test log2FC)

輸出：
  - Table_GSE89632_vs_Arendt2015_GeneLevel.xlsx
  - Table_GSE89632_vs_Arendt2015_NodeRatios.xlsx
  - (console 詳細報告)

設定（修改 CONFIG）：
================================================================
"""
import os
import sys
import numpy as np
import pandas as pd

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))

# ================================================================
# ★ CONFIG — 修改這裡 ★
# ================================================================
NODE_RATIOS_XLSX  = os.path.join(SCRIPT_DIR, 'GSE89632_node_initial_ratios.xlsx')
DEG_XLSX          = os.path.join(SCRIPT_DIR, 'GSE89632_DEG_gene_level.xlsx')
OUT_GENE_XLSX     = os.path.join(SCRIPT_DIR, 'Table_GSE89632_vs_Arendt2015_GeneLevel.xlsx')
OUT_NODE_XLSX     = os.path.join(SCRIPT_DIR, 'Table_GSE89632_vs_Arendt2015_NodeRatios.xlsx')
# ================================================================

# ----------------------------------------------------------------
# Arendt 2015 paper 數值（Table 3，NASH vs HC）
# fold_change > 0 = up, < 0 = down（已統一為 linear fold change）
# 僅收錄 Table 3 中有明確 fold change 的基因
# ----------------------------------------------------------------
ARENDT_TABLE3 = {
    # Fibrosis set
    'IL6':      {'table': 'Table3-Fibrosis',    'NASH_vs_HC': -8.96},
    'MYC':      {'table': 'Table3-Fibrosis',    'NASH_vs_HC': -7.43},
    'THBS1':    {'table': 'Table3-Fibrosis',    'NASH_vs_HC': -3.50},
    'TGFB3':    {'table': 'Table3-Fibrosis',    'NASH_vs_HC': -3.42},
    'CCL2':     {'table': 'Table3-Fibrosis',    'NASH_vs_HC': -3.23},
    'IL1B':     {'table': 'Table3-Fibrosis',    'NASH_vs_HC': -3.17},
    'SERPINE1': {'table': 'Table3-Fibrosis',    'NASH_vs_HC': -2.86},
    'IL10':     {'table': 'Table3-Fibrosis',    'NASH_vs_HC': -2.29},
    'JUN':      {'table': 'Table3-Fibrosis',    'NASH_vs_HC': -2.28},
    'CCL3':     {'table': 'Table3-Fibrosis',    'NASH_vs_HC': -2.26},
    'SERPINA1': {'table': 'Table3-Fibrosis',    'NASH_vs_HC': -2.05},
    'CXCR4':    {'table': 'Table3-Fibrosis',    'NASH_vs_HC': -2.00},
    # Chronic inflammatory response
    'THBS1_b':  {'table': 'Table3-ChronicInflam', 'NASH_vs_HC': -3.50, 'gene': 'THBS1'},
    'IL1RN':    {'table': 'Table3-ChronicInflam', 'NASH_vs_HC': -3.58},
    'S100A8':   {'table': 'Table3-ChronicInflam', 'NASH_vs_HC': -2.98},
    # Oxidative stress
    'FOS':      {'table': 'Table3-OxidStress',  'NASH_vs_HC': -7.39},
    'KLF4':     {'table': 'Table3-OxidStress',  'NASH_vs_HC': -3.66},
    'NR4A2':    {'table': 'Table3-OxidStress',  'NASH_vs_HC': -3.35},
    'FOXO1':    {'table': 'Table3-OxidStress',  'NASH_vs_HC': -2.18},
    # Lipogenesis
    'SIK1':     {'table': 'Table3-Lipogenesis', 'NASH_vs_HC': -5.29},
    'PTGS2':    {'table': 'Table3-Lipogenesis', 'NASH_vs_HC': -4.50},
    'CYR61':    {'table': 'Table3-Lipogenesis', 'NASH_vs_HC': -3.54},
    'AGPAT9':   {'table': 'Table3-Lipogenesis', 'NASH_vs_HC': -3.17},
    'GCK':      {'table': 'Table3-Lipogenesis', 'NASH_vs_HC':  4.35},
    'CYP7A1':   {'table': 'Table3-Lipogenesis', 'NASH_vs_HC':  8.34},
    'FADS1':    {'table': 'Table3-Lipogenesis', 'NASH_vs_HC':  2.49},
    'FADS2':    {'table': 'Table3-Lipogenesis', 'NASH_vs_HC':  2.81},
    'PNPLA3':   {'table': 'Table3-Lipogenesis', 'NASH_vs_HC':  2.38},
    'HNF4A':    {'table': 'Table3-Lipogenesis', 'NASH_vs_HC':  2.20},
    # Long-chain & unsaturated FA metabolism
    'TNFRSF1A': {'table': 'Table3-LongChainFA', 'NASH_vs_HC': -2.16},
    # Table 2 (NASH vs SS and HC)
    'AKR1B10':  {'table': 'Table2-NASHvsHC',   'NASH_vs_HC':  9.95},
    'EEF1A2':   {'table': 'Table2-NASHvsHC',   'NASH_vs_HC':  4.88},
    'STMN2':    {'table': 'Table2-NASHvsHC',   'NASH_vs_HC':  2.05},
    'PEG10':    {'table': 'Table2-NASHvsHC',   'NASH_vs_HC':  4.43},
    'FMO1':     {'table': 'Table2-NASHvsHC',   'NASH_vs_HC':  5.07},
    # Additional key genes from Supp Table S4a
    'FOSB':     {'table': 'SuppS4a',            'NASH_vs_HC': -25.8},
    'JUNB':     {'table': 'SuppS4a',            'NASH_vs_HC': -7.2},
    'SOCS3':    {'table': 'SuppS4a',            'NASH_vs_HC': -9.3},
    'GADD45G':  {'table': 'SuppS4a',            'NASH_vs_HC': -6.7},
    'KLF5':     {'table': 'SuppS4a',            'NASH_vs_HC': -5.3},  # KLF5 = FOSL1?
    'NR4A1':    {'table': 'SuppS4a',            'NASH_vs_HC': -5.8},
    'IL1RL1':   {'table': 'SuppS4a',            'NASH_vs_HC': -8.7},
    'MT1A':     {'table': 'SuppS4a',            'NASH_vs_HC': -5.1},
    'MT1B':     {'table': 'SuppS4a',            'NASH_vs_HC': -3.3},
    # IL8 = CXCL8 in paper, stored as 'IL8' in our DEG file
    'IL8':      {'table': 'SuppS4a',            'NASH_vs_HC': -7.4},
}


def dir_label(fc):
    """Linear fold change → direction string."""
    if fc > 1.05:
        return 'up'
    if fc < 0.95:
        return 'dn'
    return 'nc'


def fc_to_log2(fc):
    """Paper linear FC → log2 (paper uses ±x as x-fold up/down)."""
    if fc > 0:
        return np.log2(fc)
    else:  # negative means down-regulation
        return np.log2(abs(fc)) * -1


def direction_from_log2(l2):
    """log2FC → direction string."""
    if l2 > 0.07:
        return 'up'
    if l2 < -0.07:
        return 'dn'
    return 'nc'


def main():
    print('=' * 70)
    print('GSE89632 vs Arendt 2015 — 驗證比對')
    print('=' * 70)

    # --- Load files ---
    for f in [NODE_RATIOS_XLSX, DEG_XLSX]:
        exists = os.path.exists(f)
        print(f'  {"OK  " if exists else "MISS"} {os.path.basename(f)}: {f}')
        if not exists:
            raise FileNotFoundError(f)

    node_df = pd.read_excel(NODE_RATIOS_XLSX)
    deg_df  = pd.read_excel(DEG_XLSX)

    # Normalize gene symbol column
    if 'gene_symbol' not in deg_df.columns:
        print('[ERROR] DEG file needs gene_symbol column')
        sys.exit(1)

    deg_df['gene_symbol'] = deg_df['gene_symbol'].astype(str).str.strip().str.upper()
    deg_lookup = {}
    for _, row in deg_df.iterrows():
        g = row['gene_symbol']
        if g not in deg_lookup:  # keep first occurrence (best sorted by padj ideally)
            deg_lookup[g] = row

    # ============================================================
    # PART 1: Gene-level comparison vs Arendt 2015
    # ============================================================
    print('\n' + '=' * 70)
    print('PART 1: 基因層級比對 — GSE89632 DEG vs Arendt 2015 Table 2/3')
    print('  方向判斷：|log2FC| > 0.07 為 up/dn；否則 nc')
    print('=' * 70)

    rows = []
    for gene_key, info in ARENDT_TABLE3.items():
        # resolve actual gene name (handle duplicates like THBS1_b)
        gene = info.get('gene', gene_key)
        gene_upper = gene.upper()

        paper_fc    = info['NASH_vs_HC']
        paper_log2  = fc_to_log2(paper_fc)
        paper_dir   = direction_from_log2(paper_log2)

        if gene_upper in deg_lookup:
            row = deg_lookup[gene_upper]
            our_log2  = float(row['log2FC'])
            our_padj  = float(row['padj'])
            our_dir   = direction_from_log2(our_log2)
            our_fc    = 2**our_log2 if our_log2 >= 0 else -(2**abs(our_log2))
            found     = True
        else:
            our_log2  = np.nan
            our_padj  = np.nan
            our_dir   = 'NF'
            our_fc    = np.nan
            found     = False

        match = (paper_dir == our_dir) if (found and paper_dir != 'nc') else None

        rows.append({
            'Gene':          gene,
            'Table_source':  info['table'],
            'Paper_FC':      paper_fc,
            'Paper_log2FC':  round(paper_log2, 3),
            'Paper_dir':     paper_dir,
            'Our_log2FC':    round(our_log2, 3) if found else np.nan,
            'Our_padj':      our_padj,
            'Our_dir':       our_dir,
            'Concordant':    match,
            'Found_in_DEG':  found,
        })

    gene_df = pd.DataFrame(rows)

    # Deduplicate (THBS1_b → THBS1)
    gene_df = gene_df.drop_duplicates(subset='Gene', keep='first')

    # Print table
    print(f'\n{"Gene":<12} {"Paper_FC":>10} {"Paper_dir":>10} {"Our_log2":>10} {"Our_padj":>12} {"Our_dir":>8} {"Match":>7}')
    print('-' * 75)
    n_match = 0
    n_mismatch = 0
    n_nf = 0
    for _, r in gene_df.iterrows():
        match_str = '✓' if r['Concordant'] is True else ('✗' if r['Concordant'] is False else 'NF')
        if r['Concordant'] is True:  n_match += 1
        elif r['Concordant'] is False: n_mismatch += 1
        else: n_nf += 1

        padj_str = f'{r.Our_padj:.2e}' if not pd.isna(r.Our_padj) else 'N/A'
        log2_str = f'{r.Our_log2FC:.3f}' if not pd.isna(r.Our_log2FC) else 'N/A'
        print(f'{r.Gene:<12} {r.Paper_FC:>10.2f} {r.Paper_dir:>10} {log2_str:>10} {padj_str:>12} {r.Our_dir:>8} {match_str:>7}')

    total_testable = n_match + n_mismatch
    print('-' * 75)
    print(f'\n  方向一致（Concordant）：{n_match}/{total_testable} = {100*n_match/total_testable:.1f}%  (未找到基因 = {n_nf})')
    print(f'  方向不一致：{n_mismatch}')

    if n_mismatch > 0:
        print('\n  [!] 方向不一致基因：')
        for _, r in gene_df[gene_df['Concordant'] == False].iterrows():
            print(f'      {r.Gene}: Paper={r.Paper_FC:+.2f}({r.Paper_dir}), Ours log2={r.Our_log2FC:.3f}({r.Our_dir})')

    # ============================================================
    # PART 2: Node initial ratios — key nodes vs Arendt 2015 direction
    # ============================================================
    print('\n' + '=' * 70)
    print('PART 2: Node Initial Ratios — 關鍵 ODE nodes 與 paper 方向比對')
    print('  依據 Arendt 2015 Table 3 報告的 NASH vs HC 方向')
    print('=' * 70)

    # Key node → expected direction from paper
    NODE_EXPECTATIONS = {
        'AP_1':    {'genes': 'FOS, JUN', 'paper_dir': 'dn',  'paper_note': 'FOS -7.4x, JUN -2.3x'},
        'NF_kB':   {'genes': 'NFKB1',   'paper_dir': 'dn',  'paper_note': 'NFKB2 -2.1x (Supp S4a)'},
        'PPAR_a':  {'genes': 'PPARA',    'paper_dir': 'up',  'paper_note': 'FADS1/2 up; PPARA not in Table3 directly'},
        'TGF_b1':  {'genes': 'TGFB1',   'paper_dir': 'dn',  'paper_note': 'TGFB3 -3.4x; TGFB1 mRNA-protein discordance'},
        'TNFR1':   {'genes': 'TNFRSF1A','paper_dir': 'dn',  'paper_note': 'TNFRSF1A -2.16x in Table3'},
        'IL_6':    {'genes': 'IL6',      'paper_dir': 'dn',  'paper_note': 'IL6 -8.96x in Table3'},
        'IL_1':    {'genes': 'IL1A,IL1B','paper_dir': 'dn',  'paper_note': 'IL1B -3.17x in Table3'},
        'CASP7':   {'genes': 'CASP7',    'paper_dir': 'nc',  'paper_note': 'Not in Table3; microarray sensitivity limitation'},
        'CASP3':   {'genes': 'CASP3',    'paper_dir': 'up',  'paper_note': 'CASP3 up trend (Supp S4c)'},
        'IL_8':    {'genes': 'CXCL8',    'paper_dir': 'dn',  'paper_note': 'IL8 -7.4x (Supp S4a); direction reversal vs main model'},
        'TNFa':    {'genes': 'TNF',      'paper_dir': 'nc',  'paper_note': 'TNF not significant in Arendt 2015'},
        'SOCS3':   {'genes': 'SOCS3',    'paper_dir': 'dn',  'paper_note': 'SOCS3 -9.3x (Supp S4a)'},
        'FOXO1':   {'genes': 'FOXO1',    'paper_dir': 'dn',  'paper_note': 'FOXO1 -2.18x in Table3'},
        'ACDC':    {'genes': 'ADIPOQ',   'paper_dir': 'dn',  'paper_note': 'ADIPOQ expected down in NASH'},
        'SREBP_1c':{'genes': 'SREBF1',   'paper_dir': 'up',  'paper_note': 'PNPLA3 +2.38x; lipogenesis up'},
        'PPAR_g':  {'genes': 'PPARG',    'paper_dir': 'nc',  'paper_note': 'NR5A2 down (Supp S4c) but PPARG not directly in tables'},
    }

    node_lookup = {}
    for _, r in node_df.iterrows():
        node_lookup[r['node']] = r

    node_rows = []
    print(f'\n{"Node":<12} {"Ratio":>8} {"Our_dir":>8} {"Paper_dir":>10} {"Match":>7}  Paper note')
    print('-' * 80)
    n_match_n = 0
    n_mismatch_n = 0

    for node, exp in NODE_EXPECTATIONS.items():
        if node in node_lookup:
            r = node_lookup[node]
            ratio = float(r['ratio'])
            our_dir = dir_label(ratio)
        else:
            ratio = np.nan
            our_dir = 'NF'

        paper_dir = exp['paper_dir']
        # For nc paper expectation, we call it informational only
        if paper_dir != 'nc' and our_dir not in ('NF', 'nc'):
            match = (our_dir == paper_dir)
            match_str = '✓' if match else '✗'
            if match: n_match_n += 1
            else: n_mismatch_n += 1
        elif paper_dir == 'nc':
            match_str = '—'  # not tested
            match = None
        else:
            match_str = '?'
            match = None

        ratio_str = f'{ratio:.3f}' if not np.isnan(ratio) else 'N/A'
        print(f'{node:<12} {ratio_str:>8} {our_dir:>8} {paper_dir:>10} {match_str:>7}  {exp["paper_note"]}')

        node_rows.append({
            'Node':         node,
            'Genes':        exp['genes'],
            'Our_ratio':    round(ratio, 4) if not np.isnan(ratio) else np.nan,
            'Our_dir':      our_dir,
            'Paper_dir':    paper_dir,
            'Concordant':   match,
            'Paper_note':   exp['paper_note'],
        })

    total_testable_n = n_match_n + n_mismatch_n
    print('-' * 80)
    if total_testable_n > 0:
        print(f'\n  節點方向一致：{n_match_n}/{total_testable_n} = {100*n_match_n/total_testable_n:.1f}%')

    # ============================================================
    # PART 3: ODE validation-critical nodes
    # ============================================================
    print('\n' + '=' * 70)
    print('PART 3: ODE 驗證關鍵節點快查')
    print('=' * 70)

    critical_genes = {
        'CASP7':    'P_Hepatocyte_injury 主要 upstream（跨 platform sensitivity 問題）',
        'TNF':      'P_Cell_death upstream',
        'FASLG':    'P_Cell_death upstream (FasL)',
        'IL8':      'P_Inflammation upstream（已知方向反轉）',
        'TGFB1':    'P_Inflammation upstream（mRNA-protein discordance）',
        'CXCR4':    'P_Inflammation 相關（Arendt Table3 dn 佐證）',
    }

    print(f'\n{"Gene":<10} {"Our_log2FC":>12} {"Our_ratio":>12} {"Our_dir":>8}  ODE relevance')
    print('-' * 80)
    for g, note in critical_genes.items():
        gu = g.upper()
        if gu in deg_lookup:
            row = deg_lookup[gu]
            l2 = float(row['log2FC'])
            ratio = 2**l2
            d = direction_from_log2(l2)
        else:
            l2 = np.nan
            ratio = np.nan
            d = 'NF'
        l2_str = f'{l2:.3f}' if not np.isnan(l2) else 'N/A'
        r_str = f'{ratio:.4f}' if not np.isnan(ratio) else 'N/A'
        print(f'{g:<10} {l2_str:>12} {r_str:>12} {d:>8}  {note}')

    # ============================================================
    # Save Excel
    # ============================================================
    with pd.ExcelWriter(OUT_GENE_XLSX, engine='openpyxl') as writer:
        gene_df.to_excel(writer, sheet_name='GeneLevel_Comparison', index=False)
        summary = pd.DataFrame([
            ['Script', os.path.basename(__file__)],
            ['DEG file', os.path.basename(DEG_XLSX)],
            ['Comparison', 'NASH vs HC (treatment=NASH)'],
            ['Concordance', f'{n_match}/{total_testable} ({100*n_match/total_testable:.1f}%)'],
            ['Not found', str(n_nf)],
            ['Method', 'Direction agreement: |log2FC|>0.07 = up/dn'],
        ], columns=['Parameter', 'Value'])
        summary.to_excel(writer, sheet_name='Summary', index=False)
    print(f'\nGene-level Excel saved: {OUT_GENE_XLSX}')

    node_out = pd.DataFrame(node_rows)
    with pd.ExcelWriter(OUT_NODE_XLSX, engine='openpyxl') as writer:
        node_out.to_excel(writer, sheet_name='NodeRatio_Comparison', index=False)
    print(f'Node ratio Excel saved: {OUT_NODE_XLSX}')

    print('\n[DONE]')


if __name__ == '__main__':
    main()
