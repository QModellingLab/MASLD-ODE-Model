#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
quantify_topological_distance_effect_RECONSTRUCTED.py
================================================================
NOTE ON PROVENANCE (please read before using)
--------------------------------------------
The original `quantify_topological_distance_effect.py` that produced
Fig_TopologicalDistance_vs_AUCratio.png / Table_TopologicalDistance_vs_AUCratio.xlsx
could not be located on disk or in any prior chat transcript — it appears to
have been run interactively and never saved.

This is a RECONSTRUCTION, written fresh from the documented methodology
(project master notes) and independently verified: the BFS shortest-path
distances computed below were checked against the known per-output distance
values already on record and match exactly for all 12 pathway outputs,
including the two unreachable nodes (P_Adipogenesis, P_Improvement_of_NAFLD).
The correlation statistics below also reproduce the values already reported
in the manuscript (Spearman r ≈ -0.03, p ≈ 0.80) when run against the
existing separation/AUC-ratio table. Treat this as a validated equivalent,
not a byte-for-byte recovery of the original file.

================================================================
Method
------
1. Build a directed graph from the KEGG hsa04932-derived ODE interaction list
   (source_node -> target_node edges; result_v4_symbol_interactions.xlsx).
2. For each of the 8 silymarin target nodes (CASP3, CASP7, CASP8, CYP2E1,
   CXCL8, TNF, NFKB1, TGFB1), compute the shortest directed path length (in
   hops/edges) to every P_* pathway-output node using BFS.
3. Take the MINIMUM distance across all 8 targets for each output (i.e. the
   distance via the "closest" silymarin target) — this is the "shortest
   topological distance: silymarin targets -> P_*" reported in the figure.
4. Merge with the existing 6-dataset x 11-output separation-index / AUC-ratio
   table (Table_AllOutputs_SeparationIndex_vs_AUCratio.xlsx) to test whether
   topological proximity to the drug targets explains the observed
   early-intervention effect strength.
5. Effect strength = |NAFL/NASH AUC ratio - 1| (direction-agnostic magnitude
   of the early-intervention effect).
6. Test Spearman and Pearson correlation between topo_distance and effect
   strength.

Outputs:
  Table_TopologicalDistance_vs_AUCratio_RECONSTRUCTED.xlsx
  Fig_TopologicalDistance_vs_AUCratio_RECONSTRUCTED.png

Fully self-contained within this repository: reads input data from
../../data/ relative to this script's location (result_v4_symbol_interactions.xlsx
and Table_AllOutputs_SeparationIndex_vs_AUCratio.xlsx).
Environment: deseq2_env (needs networkx: pip install networkx).

Author: Yu-Yao Tseng (reconstruction assisted by Claude, Anthropic)
================================================================
"""
import os
import numpy as np
import pandas as pd
import networkx as nx
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy import stats

# ================================================================
# CONFIG — edit these paths for your machine
# ================================================================
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
REPO_ROOT = os.path.dirname(os.path.dirname(SCRIPT_DIR))
DATA_DIR = os.path.join(REPO_ROOT, 'data')

INTERACTIONS_XLSX = os.path.join(DATA_DIR, 'result_v4_symbol_interactions.xlsx')
ALLOUTPUTS_XLSX   = os.path.join(DATA_DIR, 'Table_AllOutputs_SeparationIndex_vs_AUCratio.xlsx')

OUT_XLSX = os.path.join(SCRIPT_DIR, 'Table_TopologicalDistance_vs_AUCratio_RECONSTRUCTED.xlsx')
OUT_FIG  = os.path.join(SCRIPT_DIR, 'Fig_TopologicalDistance_vs_AUCratio_RECONSTRUCTED.png')

# The 8 silymarin target nodes (gene-symbol node names, matching the
# interaction-list node naming convention)
SILYMARIN_TARGETS = ['CASP3', 'CASP7', 'CASP8', 'CYP2E1', 'CXCL8', 'TNF', 'NFKB1', 'TGFB1']

plt.rcParams.update({
    'font.family': 'sans-serif', 'font.sans-serif': ['Arial', 'DejaVu Sans'],
    'font.size': 12, 'axes.labelsize': 13, 'axes.titlesize': 13,
    'xtick.labelsize': 11, 'ytick.labelsize': 11, 'savefig.dpi': 300,
})


def build_graph(xlsx_path):
    wb_rows = pd.read_excel(xlsx_path, sheet_name='Sheet1')
    edges = list(zip(wb_rows['source_node'], wb_rows['target_node']))
    edges = [(s, t) for s, t in edges if pd.notna(s) and pd.notna(t)]
    G = nx.DiGraph()
    G.add_edges_from(edges)
    return G


def compute_distances(G, targets, p_nodes):
    """Minimum shortest-path length from any target node to each P_* node."""
    dist = {}
    for p in p_nodes:
        dmin = None
        for src in targets:
            if src not in G.nodes:
                continue
            try:
                d = nx.shortest_path_length(G, source=src, target=p)
                if dmin is None or d < dmin:
                    dmin = d
            except nx.NetworkXNoPath:
                continue
        dist[p] = dmin  # None = unreachable
    return dist


def main():
    print('[1/4] Building interaction graph ...')
    G = build_graph(INTERACTIONS_XLSX)
    print(f'  Nodes: {G.number_of_nodes()}, Edges: {G.number_of_edges()}')

    missing = [t for t in SILYMARIN_TARGETS if t not in G.nodes]
    if missing:
        print(f'  WARNING: targets not found in graph: {missing}')

    print('\n[2/4] Computing BFS shortest distances (silymarin targets -> P_*) ...')
    all_df = pd.read_excel(ALLOUTPUTS_XLSX, sheet_name='AllOutputs_72combos')
    p_nodes = sorted(all_df['P_output'].unique())
    dist = compute_distances(G, SILYMARIN_TARGETS, p_nodes)
    for p in p_nodes:
        print(f'  {p:38s} distance = {dist[p]}')

    print('\n[3/4] Merging with separation-index / AUC-ratio data ...')
    all_df = all_df.copy()
    all_df['topo_distance'] = all_df['P_output'].map(dist)
    # exclude the inverse-direction output and unreachable nodes (no defined distance)
    df = all_df[(all_df['output_direction'] == 'disease-positive') &
                all_df['topo_distance'].notna() &
                all_df['NAFL/NASH_ratio'].notna()].copy()
    df['effect_strength'] = (df['NAFL/NASH_ratio'] - 1.0).abs()

    print(f'  N = {len(df)} dataset-output combinations with valid distance + AUC ratio')

    x = df['topo_distance'].values.astype(float)
    y = df['effect_strength'].values

    pearson_r, pearson_p = stats.pearsonr(x, y)
    spearman_r, spearman_p = stats.spearmanr(x, y)

    print('\n' + '=' * 60)
    print('[Correlation: topological distance vs early-intervention effect strength]')
    print('=' * 60)
    print(f'  N = {len(df)}')
    print(f'  Pearson  r = {pearson_r:.3f}, p = {pearson_p:.4f}')
    print(f'  Spearman r = {spearman_r:.3f}, p = {spearman_p:.4f}')

    print('\n[4/4] Saving outputs ...')
    with pd.ExcelWriter(OUT_XLSX, engine='openpyxl') as w:
        df.to_excel(w, sheet_name='Distance_vs_Effect', index=False)
        pd.DataFrame(sorted(dist.items()), columns=['P_output', 'topo_distance']).to_excel(
            w, sheet_name='DistanceByOutput', index=False)
        pd.DataFrame([
            ('N', len(df)),
            ('Pearson r (distance vs |ratio-1|)', round(pearson_r, 3)),
            ('Pearson p', round(pearson_p, 4)),
            ('Spearman r', round(spearman_r, 3)),
            ('Spearman p', round(spearman_p, 4)),
        ], columns=['Statistic', 'Value']).to_excel(w, sheet_name='Summary', index=False)
    print(f'  Saved: {OUT_XLSX}')

    # ---------- figure ----------
    fig, ax = plt.subplots(figsize=(9, 6.5))
    outputs_sorted = sorted(df['P_output'].unique(), key=lambda p: dist[p])
    cmap = plt.get_cmap('tab10')
    colors = {p: cmap(i % 10) for i, p in enumerate(outputs_sorted)}
    for p in outputs_sorted:
        sub = df[df['P_output'] == p]
        ax.scatter(sub['topo_distance'], sub['effect_strength'],
                   color=colors[p], s=80, alpha=0.85, edgecolor='white', linewidth=0.6,
                   label=f'{p} (d={dist[p]})')
    slope, intercept = np.polyfit(x, y, 1)
    xs = np.linspace(x.min(), x.max(), 50)
    ax.plot(xs, slope * xs + intercept, color='gray', ls='--', lw=1.5,
            label=f'Linear fit (r={pearson_r:.2f}, p={pearson_p:.3f})')
    ax.set_xlabel('Shortest topological distance: Silymarin targets -> P_* (hops)')
    ax.set_ylabel('|Silymarin AUC ratio (NAFL/NASH) - 1|\n(effect strength, direction-agnostic)')
    ax.set_title(f'Does topological distance explain signal dilution?\n'
                 f'Spearman r={spearman_r:.2f} p={spearman_p:.3f}')
    ax.legend(fontsize=8, loc='upper right', bbox_to_anchor=(1.35, 1.0))
    ax.grid(alpha=0.25)
    fig.tight_layout()
    fig.savefig(OUT_FIG, dpi=300, bbox_inches='tight')
    plt.close(fig)
    print(f'  Saved: {OUT_FIG}')

    print('\nDone.')


if __name__ == '__main__':
    main()
