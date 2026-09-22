#!/usr/bin/env python3
"""
GSE126848_Fig2_volcano_pyDESeq2.py
================================
Fig 2 Volcano + three-group pyDESeq2 comparison (Python-only, no batch correction)

Pipeline:
  1. pyDESeq2: NASH vs Normal, NAFL vs Normal, Obese vs Normal
  2. Merge results -> Excel (includes ratio column, used for ODE initial values + Fig3 enrichment)
  3. Volcano plot (NASH vs Normal, S1 genes highlighted)
  4. Concordance check vs Suppli Fig.1C

Input:  GSE126848_Gene_counts_mapped.xlsx
        Fig1C_pixel_values.xlsx (optional)
        suppli_S1_genes.csv (optional)
Output: DEG_pydeseq2_3comparisons.xlsx (three-group comparison + Merged + Summary)
        DEG_pydeseq2_NASH_vs_Normal.csv (NASH vs Normal, kept for backward compatibility with Fig3)
        Fig2_volcano_pyDESeq2.png / .pdf

Run in Spyder (deseq2_env).
Requirements: pip install pydeseq2 adjustText openpyxl

Author: Yu-Yao Tseng  |  Date: 2026-04-15
"""

import sys, os
import pandas as pd
import numpy as np
from collections import Counter
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import warnings; warnings.filterwarnings('ignore')

plt.rcParams.update({
    'font.family': 'sans-serif',
    'font.sans-serif': ['Calibri', 'Arial', 'DejaVu Sans'],
})
try:
    from adjustText import adjust_text
    HAS_AT = True
except ImportError:
    HAS_AT = False
    print("⚠ adjustText not found → pip install adjustText")

# ============================================================
# Configuration
# ============================================================

SCRIPT_DIR  = os.path.dirname(os.path.abspath(__file__))
COUNTS_FILE = os.path.join(SCRIPT_DIR, "GSE126848_Gene_counts_mapped.xlsx")
PIXEL_XLSX  = os.path.join(SCRIPT_DIR, "Fig1C_pixel_values.xlsx")
S1_CSV      = os.path.join(SCRIPT_DIR, "suppli_S1_genes.csv")
OUTPUT_DIR  = SCRIPT_DIR

P_THR  = 0.05
FC_THR = 1.0

COMPARISONS = [
    ('NASH',  'Normal'),
    ('NAFL',  'Normal'),
    ('Obese', 'Normal'),
]

GS_COLORS = {
    'Monocyte recruitment':      '#8E24AA',
    'Inflammation signaling':    '#C62828',
    'Hepatocellular cell death': '#E65100',
    'Stellate cell activation':  '#2E7D32',
    'Lipid metabolism':          '#1565C0',
    'Insulin signaling':         '#00838F',
    'FXR signaling':             '#6D4C41',
}

PIX2HGNC = {
    'A-SMA':'ACTA2','APOCII':'APOC2','APOCIII':'APOC3','BACS':'ABCB4',
    'CPT1':'CPT1A','F4/80_EMR1':'ADGRE1','FXR':'NR1H4','G6PASE':'G6PC',
    'IKK':'CHUK','JNK':'MAPK8','LGALS3/MAC-2':'LGALS3','MCP-1':'CCL2',
    'MDR3':'ABCB4','MEK2':'MAP2K2','NLRP1B':'RIPK1','OAT2':'SLC22A7',
    'OSTB':'SLC51B','PEPCK':'PCK1','PYG':'PYGL','RIP1':'RIPK1','SCD1':'SCD',
    'TGFB':'TGFB1','TGFB_stellate':'TGFB1','TNFR':'TNFRSF1A',
}


# ============================================================
# Part 1: pyDESeq2
# ============================================================

def load_counts():
    df = pd.read_excel(COUNTS_FILE)
    sample_cols = [c for c in df.columns
                   if any(c.startswith(p) for p in ['Normal','Obese','NAFL','NASH'])]
    gene_sym = df['Gene_Symbol'].values
    cnt = Counter(gene_sym); dups = {g for g,c in cnt.items() if c>1}; seen = Counter()
    gs_uniq = []
    for g in gene_sym:
        if g in dups: seen[g]+=1; gs_uniq.append(f'{g}_{seen[g]}')
        else: gs_uniq.append(g)
    counts_df = pd.DataFrame(df[sample_cols].values.astype(int), index=gs_uniq, columns=sample_cols)
    return counts_df, gs_uniq, sample_cols


def get_group_cols(sample_cols, group_name):
    prefix_map = {'Normal':'Normal', 'Obese':'Obese', 'NAFL':'NAFL_', 'NASH':'NASH_'}
    prefix = prefix_map.get(group_name, group_name)
    return [c for c in sample_cols if c.startswith(prefix)]


def run_one_comparison(counts_df, sample_cols, disease, control):
    from pydeseq2.dds import DeseqDataSet
    from pydeseq2.ds import DeseqStats

    dis_cols = get_group_cols(sample_cols, disease)
    ctl_cols = get_group_cols(sample_cols, control)
    use_cols = ctl_cols + dis_cols
    sub_counts = counts_df[use_cols]
    condition = [control]*len(ctl_cols) + [disease]*len(dis_cols)
    metadata = pd.DataFrame({'condition': condition}, index=use_cols)
    print(f"    {disease} (n={len(dis_cols)}) vs {control} (n={len(ctl_cols)})")

    dds = DeseqDataSet(counts=sub_counts.T, metadata=metadata, design="~condition", refit_cooks=True)
    dds.fit_size_factors()
    dds.fit_genewise_dispersions()
    dds.fit_dispersion_trend()
    dds.fit_dispersion_prior()
    dds.fit_MAP_dispersions()
    dds.fit_LFC()
    dds.calculate_cooks()
    dds.refit()

    stat = DeseqStats(dds, contrast=["condition", disease, control])
    stat.summary()

    res = stat.results_df.copy()
    res['Gene'] = res.index
    res = res.reset_index(drop=True)
    res['ratio'] = np.power(2, res['log2FoldChange'])

    valid = res.dropna(subset=['padj'])
    sig = valid[(valid['padj']<P_THR) & (valid['log2FoldChange'].abs()>FC_THR)]
    n_up = (sig['log2FoldChange']>0).sum(); n_dn = (sig['log2FoldChange']<0).sum()
    print(f"    → {len(valid)} tested, {len(sig)} DEGs (↑{n_up} ↓{n_dn})")
    return res


def run_all_comparisons():
    print(f"\n{'='*60}")
    print(f"  Part 1: pyDESeq2 — three-group comparison")
    print(f"{'='*60}")
    counts_df, gs_uniq, sample_cols = load_counts()
    print(f"  Loaded: {counts_df.shape[0]} genes × {counts_df.shape[1]} samples")
    results = {}
    for disease, control in COMPARISONS:
        label = f"{disease}_vs_{control}"
        print(f"\n  [{label}]")
        results[label] = run_one_comparison(counts_df, sample_cols, disease, control)
    return results


# ============================================================
# Part 2: Save Excel + CSV
# ============================================================

def save_results(results):
    print(f"\n{'='*60}")
    print(f"  Part 2: Save Excel + CSV")
    print(f"{'='*60}")

    xlsx = os.path.join(OUTPUT_DIR, "DEG_pydeseq2_3comparisons.xlsx")
    with pd.ExcelWriter(xlsx, engine='openpyxl') as w:
        for label, df in results.items():
            df.to_excel(w, sheet_name=label, index=False)

        # Merged sheet
        base = results['NASH_vs_Normal'][['Gene','baseMean']].copy()
        for label, df in results.items():
            short = label.split('_vs_')[0]
            base[f'log2FC_{short}'] = df['log2FoldChange'].values
            base[f'padj_{short}'] = df['padj'].values
            base[f'ratio_{short}'] = df['ratio'].values
        base.to_excel(w, sheet_name='Merged_3grp', index=False)

        # Summary
        rows = [
            ('Source', os.path.basename(COUNTS_FILE)),
            ('Method', 'pyDESeq2 (no batch correction)'),
            ('Thresholds', f'|log2FC| > {FC_THR}, padj < {P_THR}'),
            ('', ''),
        ]
        for label, df in results.items():
            valid = df.dropna(subset=['padj'])
            sig = valid[(valid['padj']<P_THR) & (valid['log2FoldChange'].abs()>FC_THR)]
            n_up = (sig['log2FoldChange']>0).sum(); n_dn = (sig['log2FoldChange']<0).sum()
            rows.append((label, f'{len(sig)} DEGs (↑{n_up} ↓{n_dn})'))
        pd.DataFrame(rows, columns=['Item','Value']).to_excel(w, sheet_name='Summary', index=False)

    print(f"  ✓ {xlsx}")

    # CSV for backward compatibility
    csv_out = os.path.join(OUTPUT_DIR, "DEG_pydeseq2_NASH_vs_Normal.csv")
    nash = results['NASH_vs_Normal'].copy()
    nash[['baseMean','log2FoldChange','lfcSE','stat','pvalue','padj','Gene']].to_csv(csv_out, index=False)
    print(f"  ✓ {csv_out}")


# ============================================================
# Part 3: Concordance
# ============================================================

def concordance_check(res):
    if not os.path.exists(PIXEL_XLSX):
        return 0, 0
    pixel = pd.read_excel(PIXEL_XLSX, sheet_name=0)
    cm, ct = 0, 0
    for _, pr in pixel.iterrows():
        pd_ = pr['Direction']
        if pd_ not in ('↑','↓'): continue
        hgnc = PIX2HGNC.get(pr['Gene'], pr['Gene'])
        m = res[res['Gene']==hgnc]
        if len(m)==0: continue
        ct += 1
        our = '↑' if m.iloc[0]['log2FoldChange']>0 else '↓'
        if our==pd_: cm += 1
    return cm, ct


# ============================================================
# Part 4: Volcano
# ============================================================

def make_volcano(res, n_up, n_dn, conc_m, conc_t):
    print(f"\n{'='*60}")
    print(f"  Part 4: Volcano Plot")
    print(f"{'='*60}")

    res = res.copy()
    res['nlog10'] = -np.log10(res['padj'].clip(lower=1e-300))
    res['nlog10'] = res['nlog10'].replace([np.inf], 50)
    res['cat'] = 'NS'
    res.loc[(res['padj']<P_THR) & (res['log2FoldChange']>FC_THR), 'cat'] = 'Up'
    res.loc[(res['padj']<P_THR) & (res['log2FoldChange']<-FC_THR), 'cat'] = 'Down'

    Y_MAX = 17.5
    res_plot = res[res['nlog10'] <= Y_MAX]    # filter points beyond y-axis

    fig, ax = plt.subplots(figsize=(7.5, 3.8))
    for cat,c,s,a in [('NS','#D0D0D0',3,0.15),('Down','#4F7FD9',4,0.35),('Up','#D94F4F',4,0.35)]:
        sub = res_plot[res_plot['cat']==cat]
        ax.scatter(sub['log2FoldChange'], sub['nlog10'], c=c, s=s, alpha=a, edgecolors='none', rasterized=True)

    ax.axhline(-np.log10(P_THR), color='#888', ls='--', lw=0.8, zorder=1)
    ax.axvline(-FC_THR, color='#888', ls=':', lw=0.6, zorder=1)
    ax.axvline(FC_THR, color='#888', ls=':', lw=0.6, zorder=1)

    texts_auto=[]; auto_x=[]; auto_y=[]
    if os.path.exists(S1_CSV):
        s1 = pd.read_csv(S1_CSV)
        s1m = s1.merge(res[['Gene','log2FoldChange','padj','nlog10','cat']],
                       left_on='Symbol', right_on='Gene', how='inner')
        for _, row in s1m.iterrows():
            x, y = row['log2FoldChange'], row['nlog10']
            # Skip genes outside cropped y-axis (CD163, etc.)
            if y > Y_MAX:
                continue
            gc = GS_COLORS.get(row['GeneSet'], '#555')
            is_sig = row['cat'] != 'NS'
            ax.scatter(x, y, c=gc, s=30 if is_sig else 10, marker='D',
                       edgecolors='black' if is_sig else '#aaa',
                       linewidths=0.6 if is_sig else 0.25,
                       alpha=1.0 if is_sig else 0.45,
                       zorder=6 if is_sig else 4)
            if is_sig:
                txt = ax.text(x, y, row['Symbol'], fontsize=6.5, fontweight='bold', color=gc,
                              bbox=dict(boxstyle='round,pad=0.08', fc='white', ec='none', alpha=0.85), zorder=7)
                texts_auto.append(txt); auto_x.append(x); auto_y.append(y)

        if HAS_AT and texts_auto:
            adjust_text(texts_auto, x=auto_x, y=auto_y, ax=ax,
                arrowprops=dict(arrowstyle='->', color='#888', lw=0.35, shrinkA=4, shrinkB=3),
                force_text=(0.4,0.5), force_points=(0.3,0.4),
                expand_text=(1.2,1.3), expand_points=(1.3,1.3), only_move={'text':'xy'})

    ax.set_xlabel('log₂ Fold Change (NASH / Normal)', fontsize=10)
    ax.set_ylabel('-log₁₀(adjusted P-value)', fontsize=10)
    ax.set_xlim(-6, 6)
    ax.set_ylim(-0.5, Y_MAX)

    legend_el = [Line2D([0],[0], marker='D', color='w', markerfacecolor=gc,
                        markeredgecolor='k', markeredgewidth=0.3, markersize=5.5,
                        label=gs) for gs, gc in GS_COLORS.items()]
    ax.legend(handles=legend_el, loc='upper right', framealpha=0.92, fontsize=5.5,
              handletextpad=0.3, borderpad=0.5, labelspacing=0.35)

    ax.text(0.02, 0.98, f'DEGs: {n_up+n_dn} (↑{n_up} ↓{n_dn})\n|log₂FC|>{FC_THR}, padj<{P_THR}',
            transform=ax.transAxes, fontsize=6.5, va='top',
            bbox=dict(boxstyle='round', fc='white', alpha=0.85))

    for sp in ax.spines.values(): sp.set_linewidth(0.4)
    ax.grid(True, alpha=0.06)

    for ext in ['png','pdf']:
        out = os.path.join(OUTPUT_DIR, f'Fig2_volcano_pyDESeq2.{ext}')
        fig.savefig(out, dpi=300 if ext=='png' else None, bbox_inches='tight', facecolor='white')
        print(f"  Saved: {out}")
    plt.close(fig)


# ============================================================
# Main
# ============================================================

def main():
    print("="*60)
    print("  Fig 2: pyDESeq2 three-group comparison + Volcano")
    print("="*60)

    results = run_all_comparisons()
    save_results(results)

    nash_res = results['NASH_vs_Normal'].dropna(subset=['padj'])
    conc_m, conc_t = concordance_check(nash_res)
    if conc_t > 0:
        print(f"\n  Concordance vs Suppli Fig.1C: {conc_m}/{conc_t} ({100*conc_m/conc_t:.0f}%)")

    sig = nash_res[(nash_res['padj']<P_THR) & (nash_res['log2FoldChange'].abs()>FC_THR)]
    n_up = (sig['log2FoldChange']>0).sum(); n_dn = (sig['log2FoldChange']<0).sum()
    make_volcano(nash_res, n_up, n_dn, conc_m, conc_t)

    print(f"\n{'='*60}")
    print(f"  All done!")
    print(f"  Output:")
    print(f"    DEG_pydeseq2_3comparisons.xlsx")
    print(f"    DEG_pydeseq2_NASH_vs_Normal.csv")
    print(f"    Fig2_volcano_pyDESeq2.png / .pdf")
    print(f"{'='*60}")

if __name__ == '__main__':
    main()
