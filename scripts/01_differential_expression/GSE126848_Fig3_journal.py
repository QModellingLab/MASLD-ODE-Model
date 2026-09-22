#!/usr/bin/env python3
"""
GSE126848_Fig3_enrichment_dotplot.py
=====================================
Fig 3: Functional Enrichment Analysis (pyDESeq2 DEGs, no BC)

一支腳本完成：
  1. 讀取 DEG_pydeseq2_NASH_vs_Normal.csv
  2. Enrichr API 跑 KEGG + GO BP（需要網路）
  3. 產出 1×2 dot plot（左 GO BP、右 KEGG）
  4. 儲存 Enrichment Excel（含所有數據）

Input:  DEG_pydeseq2_NASH_vs_Normal.csv（與本腳本同資料夾）
Output: Enrichment_NASH_vs_Normal_pyDESeq2.xlsx
        Fig3_dotplot_pyDESeq2.png / .pdf

在 Spyder (deseq2_env) 中執行。
Requirements: pip install gseapy pandas matplotlib openpyxl

Author: Yu-Yao Tseng  |  Date: 2026-04-15
"""

import os, sys
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import warnings; warnings.filterwarnings('ignore')

plt.rcParams.update({
    'font.family': 'sans-serif',
    'font.sans-serif': ['Calibri', 'Arial', 'DejaVu Sans'],
    'font.size': 8,
    'savefig.dpi': 300,
})

try:
    import gseapy as gp
    print("✓ gseapy loaded")
except ImportError:
    raise ImportError("請先安裝: pip install gseapy")

# ============================================================
# 設定區
# ============================================================

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
DEG_FILE   = os.path.join(SCRIPT_DIR, "DEG_pydeseq2_NASH_vs_Normal.csv")
OUTPUT_DIR = SCRIPT_DIR

P_THR  = 0.05
FC_THR = 1.0

KEGG_DB  = 'KEGG_2021_Human'
GO_BP_DB = 'GO_Biological_Process_2023'

TOP_N       = 10
TERM_MAXLEN = 42


# ============================================================
# Part 1: Enrichment
# ============================================================

def run_enrichr(gene_list, gene_sets, desc=""):
    if len(gene_list) < 5:
        print(f"    Skip {desc}: only {len(gene_list)} genes")
        return pd.DataFrame()
    try:
        enr = gp.enrichr(
            gene_list=gene_list, gene_sets=gene_sets,
            organism='human', outdir=None, no_plot=True, cutoff=0.1,
        )
        res = enr.results.copy().sort_values('Adjusted P-value')
        sig = res[res['Adjusted P-value'] < 0.05]
        if len(sig) == 0:
            sig = res[res['Adjusted P-value'] < 0.10]
        return sig
    except Exception as e:
        print(f"    Error in {desc}: {e}")
        return pd.DataFrame()


def do_enrichment():
    print(f"\n{'='*60}")
    print(f"  Part 1: Enrichment (Enrichr API)")
    print(f"{'='*60}")

    res = pd.read_csv(DEG_FILE).dropna(subset=['padj'])
    sig = res[(res['padj'] < P_THR) & (res['log2FoldChange'].abs() > FC_THR)]
    up_genes = sig[sig['log2FoldChange'] > 0]['Gene'].dropna().tolist()
    dn_genes = sig[sig['log2FoldChange'] < 0]['Gene'].dropna().tolist()
    up_genes = [g for g in up_genes if g and str(g) != 'nan' and len(str(g)) > 1]
    dn_genes = [g for g in dn_genes if g and str(g) != 'nan' and len(str(g)) > 1]
    print(f"  DEGs: Up={len(up_genes):,}, Down={len(dn_genes):,}")

    results = {}
    for name, genes, db in [
        ('KEGG_Up',   up_genes, KEGG_DB),
        ('KEGG_Down', dn_genes, KEGG_DB),
        ('GO_BP_Up',  up_genes, GO_BP_DB),
        ('GO_BP_Down',dn_genes, GO_BP_DB),
    ]:
        print(f"  Running {name} ({len(genes)} genes)...")
        results[name] = run_enrichr(genes, db, name)
        print(f"    → {len(results[name])} significant terms")

    return results, up_genes, dn_genes


# ============================================================
# Part 2: Save Excel
# ============================================================

def save_excel(results, up_genes, dn_genes):
    print(f"\n{'='*60}")
    print(f"  Part 2: Save Excel")
    print(f"{'='*60}")

    xlsx = os.path.join(OUTPUT_DIR, "Enrichment_NASH_vs_Normal_pyDESeq2.xlsx")
    with pd.ExcelWriter(xlsx, engine='openpyxl') as w:
        for name, df in results.items():
            if len(df) > 0:
                df.to_excel(w, sheet_name=name, index=False)

        summary_rows = [
            ('DEG source',      os.path.basename(DEG_FILE)),
            ('Method',          'pyDESeq2 (no batch correction)'),
            ('Thresholds',      f'|log2FC| > {FC_THR}, padj < {P_THR}'),
            ('Up genes',        len(up_genes)),
            ('Down genes',      len(dn_genes)),
            ('KEGG DB',         KEGG_DB),
            ('GO BP DB',        GO_BP_DB),
            ('KEGG Up terms',   len(results.get('KEGG_Up', []))),
            ('KEGG Down terms', len(results.get('KEGG_Down', []))),
            ('GO BP Up terms',  len(results.get('GO_BP_Up', []))),
            ('GO BP Down terms',len(results.get('GO_BP_Down', []))),
        ]

        nafld_found = []
        for name, df in results.items():
            if len(df) == 0: continue
            nafld = df[df['Term'].str.contains(
                'fatty liver|NAFLD|NASH|Non-alcoholic|lipid|insulin|fibrosis',
                case=False, na=False)]
            for _, r in nafld.iterrows():
                nafld_found.append(f"{name}: {r['Term']} (padj={r['Adjusted P-value']:.2e})")
        if nafld_found:
            summary_rows.append(('', ''))
            summary_rows.append(('NAFLD-related terms', ''))
            for f in nafld_found:
                summary_rows.append(('', f))

        pd.DataFrame(summary_rows, columns=['Item', 'Value']).to_excel(
            w, sheet_name='Summary', index=False)
    print(f"  ✓ {xlsx}")


# ============================================================
# Part 3: Dot Plot — 字體放大 2×，dot size 正比 gene count
# ============================================================

def truncate(s, maxlen=TERM_MAXLEN):
    s = str(s)
    if '(' in s and s.endswith(')'):
        s = s[:s.rfind('(')].strip()
    return s[:maxlen] + '...' if len(s) > maxlen else s


def prep_panel(df, n=TOP_N):
    if len(df) == 0:
        return pd.DataFrame()
    df = df.sort_values('Adjusted P-value').head(n)
    out = []
    for _, r in df.iterrows():
        genes_str = str(r.get('Genes', ''))
        n_genes = len(genes_str.split(';')) if genes_str and genes_str != 'nan' else 0
        out.append({
            'term': truncate(r['Term']),
            'nlog10': -np.log10(max(r['Adjusted P-value'], 1e-300)),
            'n_genes': n_genes,
        })
    return pd.DataFrame(out)


def make_dotplot(results):
    print(f"\n{'='*60}")
    print(f"  Part 3: Dot Plot")
    print(f"{'='*60}")

    go_up   = prep_panel(results.get('GO_BP_Up', pd.DataFrame()))
    go_dn   = prep_panel(results.get('GO_BP_Down', pd.DataFrame()))
    kegg_up = prep_panel(results.get('KEGG_Up', pd.DataFrame()))
    kegg_dn = prep_panel(results.get('KEGG_Down', pd.DataFrame()))

    panels = [
        ('GO Biological Process', go_up, go_dn),
        ('KEGG Pathways',         kegg_up, kegg_dn),
    ]

    # Global max gene count for proportional sizing
    all_ng = []
    for _, up, dn in panels:
        if len(up) > 0: all_ng.extend(up['n_genes'].tolist())
        if len(dn) > 0: all_ng.extend(dn['n_genes'].tolist())
    max_ng = max(all_ng) if all_ng else 1

    # --- Dot size: area proportional to gene count ---
    # s = (n / max_n) * MAX_DOT_AREA, with minimum MIN_DOT_AREA
    MIN_DOT_AREA = 20
    MAX_DOT_AREA = 350

    def gene_count_to_size(n):
        return max(MIN_DOT_AREA, (n / max_ng) * MAX_DOT_AREA)

    # Figure — npj SBA 雙欄寬度 183mm = 7.205 inch
    n_rows_L = len(go_up) + len(go_dn)
    n_rows_R = len(kegg_up) + len(kegg_dn)
    max_rows = max(n_rows_L, n_rows_R, 1)
    fig_h = max_rows * 0.21 + 1.45  # 依寬度縮放比例(7.205/11.0≈0.655)同步調整

    fig, axes = plt.subplots(1, 2, figsize=(7.205, fig_h))
    fig.subplots_adjust(wspace=2.00, top=0.90, bottom=0.10, left=0.24, right=0.96)

    for ax, (title, up_data, dn_data) in zip(axes, panels):
        up_r = up_data.iloc[::-1].copy() if len(up_data) > 0 else pd.DataFrame()
        dn_r = dn_data.iloc[::-1].copy() if len(dn_data) > 0 else pd.DataFrame()
        n_dn = len(dn_r); n_up = len(up_r)

        all_terms=[]; all_nlog=[]; all_ng_list=[]; all_dir=[]
        for _, r in dn_r.iterrows():
            all_terms.append(r['term']); all_nlog.append(r['nlog10'])
            all_ng_list.append(r['n_genes']); all_dir.append('down')
        for _, r in up_r.iterrows():
            all_terms.append(r['term']); all_nlog.append(r['nlog10'])
            all_ng_list.append(r['n_genes']); all_dir.append('up')

        n_total = len(all_terms)
        if n_total == 0:
            ax.text(0.5, 0.5, 'No significant terms', transform=ax.transAxes,
                    ha='center', va='center', fontsize=9, color='#888')
            ax.set_title(title, fontsize=9, fontweight='bold', loc='left')
            continue

        y = list(range(n_total))
        nlog = np.array(all_nlog)
        sizes = np.array([gene_count_to_size(n) for n in all_ng_list])
        colors = ['#D94F4F' if d == 'up' else '#4F7FD9' for d in all_dir]

        ax.scatter(nlog, y, s=sizes, c=colors, edgecolors='#333',
                   linewidths=0.4, zorder=3, alpha=0.85)
        ax.set_yticks(y)
        ax.set_yticklabels(all_terms, fontsize=7)
        ax.set_title(title, fontsize=9, fontweight='bold', loc='left', pad=6)
        ax.axvline(-np.log10(0.05), color='#999', ls='--', lw=0.5, zorder=0)
        ax.tick_params(axis='y', length=0, pad=3)
        ax.tick_params(axis='x', labelsize=7)
        ax.grid(axis='x', alpha=0.10, zorder=0)

        xmax = nlog.max() if len(nlog) > 0 else 1
        ax.set_xlim(0, xmax * 1.28)

        for i in range(n_total):
            ax.text(nlog[i] + xmax * 0.03, i, str(all_ng_list[i]),
                    fontsize=6, va='center', ha='left', color='#444', fontweight='bold')

        if n_dn > 0 and n_up > 0:
            ax.axhline(n_dn - 0.5, color='#AAA', ls='-', lw=0.7, zorder=1)

    axes[0].set_xlabel('-log₁₀(Adjusted P-value)', fontsize=8)
    axes[1].set_xlabel('-log₁₀(Adjusted P-value)', fontsize=8)

    fig.suptitle('Functional Enrichment: NASH vs. Normal (pyDESeq2 DEGs)',
                 fontsize=10, fontweight='bold', x=0.5, y=0.97)

    # Legend: use same sizing function for consistency
    legend_el = [
        Line2D([0],[0], marker='o', color='w', markerfacecolor='#D94F4F',
               markeredgecolor='#333', markersize=6, label='Upregulated'),
        Line2D([0],[0], marker='o', color='w', markerfacecolor='#4F7FD9',
               markeredgecolor='#333', markersize=6, label='Downregulated'),
    ]
    # Size legend: use actual gene_count_to_size → convert to markersize
    for ns in [5, 20, 50]:
        s_area = gene_count_to_size(ns)
        ms = np.sqrt(s_area) * 0.655  # 依雙欄縮放比例同步縮小圖例 marker
        legend_el.append(
            Line2D([0],[0], marker='o', color='w', markerfacecolor='#aaa',
                   markeredgecolor='#333', markersize=ms,
                   label=f'n = {ns}'))

    fig.legend(handles=legend_el, loc='lower center', ncol=5, fontsize=7,
               frameon=False, bbox_to_anchor=(0.5, -0.03))

    for ext in ['png', 'pdf']:
        out = os.path.join(OUTPUT_DIR, f'Fig3_dotplot_journal.{ext}')
        plt.savefig(out, dpi=300 if ext == 'png' else None, bbox_inches='tight')
        print(f"  Saved: {out}")
    plt.close()


# ============================================================
# Main
# ============================================================

def main():
    print("=" * 60)
    print("  Fig 3: Enrichment + Dot Plot (all-in-one)")
    print("=" * 60)

    results, up_genes, dn_genes = do_enrichment()
    save_excel(results, up_genes, dn_genes)
    make_dotplot(results)

    print(f"\n{'='*60}")
    print(f"  All done!")
    print(f"  Output:")
    print(f"    Enrichment_NASH_vs_Normal_pyDESeq2.xlsx")
    print(f"    Fig3_dotplot_pyDESeq2.png / .pdf")
    print(f"{'='*60}")


if __name__ == '__main__':
    main()
