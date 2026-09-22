#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
quantify_topological_distance_effect.py
================================================================
計算 Silymarin 8 個標靶基因（CASP3, CASP7, CASP8, CYP2E1, IL_8,
TNFa, NF_kB, TGF_b1）到每個 P_* output 的最短拓樸距離（BFS hop
數，依方向性邊：ACTS + INHS，與 build_fig4_v11_reference.py 的
邊清單完全一致，非另行假設）。

驗證：拓樸距離是否能解釋「擴大到 11 個 output 後，分離度 vs
早期介入效果的訊號被稀釋」這個現象——距離越近的 output，訊號是否
越強（NAFL/NASH ratio 越大、或跟 SI 的相關性越強）？

需先跑過 quantify_all_outputs_separation_vs_intervention.py，
本腳本讀取其輸出的 Table_AllOutputs_SeparationIndex_vs_AUCratio.xlsx。

輸出：
  Table_TopologicalDistance_vs_AUCratio.xlsx
  Fig_TopologicalDistance_vs_AUCratio.png
================================================================
"""
import os
from collections import deque, defaultdict
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy import stats

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
IN_XLSX = os.path.join(SCRIPT_DIR, 'Table_AllOutputs_SeparationIndex_vs_AUCratio.xlsx')

OUT_XLSX = os.path.join(SCRIPT_DIR, 'Table_TopologicalDistance_vs_AUCratio.xlsx')
OUT_FIG  = os.path.join(SCRIPT_DIR, 'Fig_TopologicalDistance_vs_AUCratio.png')

# 與 build_fig4_v11_reference.py 完全一致的邊清單（已查證，非另行假設）
ACTS = [("ACDC","adipoR"),("AMPK","P_Improvement_of_NAFLD"),("AMPK","p38"),("AP_1","IL_1"),("AP_1","IL_6"),("AP_1","TNFa"),("ASK1","JNK1_2"),("ATF4","CHOP"),("Bax","Cytc"),("Bid","Bax"),("Bim","Bax"),("CASP3","P_Hepatocyte_injury"),("CASP7","P_Hepatocyte_injury"),("CASP8","Bid"),("CHOP","Bim"),("CYP2E1","FasL"),("CYP2E1","IKKb"),("CYP2E1","IL_8"),("CYP2E1","JNK1_2"),("CYP2E1","TGF_b1"),("CYP2E1","TNFa"),("C_EBPa","P_Adipogenesis"),("Cdc42","MLK3"),("ChREBP","L_PK"),("ChREBP","P_De_novo_fatty_acid_synthesis"),("CxI","FasL"),("CxI","IKKb"),("CxI","IL_8"),("CxI","JNK1_2"),("CxI","TGF_b1"),("CxI","TNFa"),("CxII","FasL"),("CxII","IKKb"),("CxII","IL_8"),("CxII","JNK1_2"),("CxII","TGF_b1"),("CxII","TNFa"),("CxIII","FasL"),("CxIII","IKKb"),("CxIII","IL_8"),("CxIII","JNK1_2"),("CxIII","TGF_b1"),("CxIII","TNFa"),("CxIV","FasL"),("CxIV","IKKb"),("CxIV","IL_8"),("CxIV","JNK1_2"),("CxIV","TGF_b1"),("CxIV","TNFa"),("Cytc","CASP3"),("Cytc","CASP7"),("Fas","CASP8"),("FasL","Fas"),("FasL","P_Cell_death"),("GSK_3","INS"),("GSK_3","P_Hyperinsulinemia"),("IKKb","NF_kB"),("IL_1","P_Development_of_steatohepatitis"),("IL_6","IL_6R"),("IL_6","P_Development_of_steatohepatitis"),("IL_6R","SOCS3"),("IL_8","P_Inflammation"),("INS","INSR"),("INSR","IRS_1_2"),("INSR","LXR_a"),("IRE1a","TRAF2"),("IRE1a","XBP1"),("IRS_1_2","PI3K"),("ITCH","CASP8"),("JNK1_2","AP_1"),("JNK1_2","ITCH"),("JNK1_2","P_Apoptosis"),("JNK1_2","P_HCC_proliferation"),("LEP","ObR"),("LXR_a","RXR"),("MLK3","JNK1_2"),("NF_kB","IL_1"),("NF_kB","IL_6"),("NF_kB","TNFa"),("ObR","AMPK"),("PERK","eIF2a"),("PI3K","Akt"),("PPAR_a","P_Improvement_of_NAFLD"),("PPAR_g","P_Development_of_NAFLD"),("RXR","PPAR_g"),("RXR","SREBP_1c"),("Rac1","MLK3"),("SOCS3","SREBP_1c"),("SREBP_1c","P_De_novo_fatty_acid_synthesis"),("TGF_b1","P_Fibrosis"),("TGF_b1","P_Inflammation"),("TNFR1","NF_kB"),("TNFa","JNK1_2"),("TNFa","P_Cell_death"),("TNFa","P_Development_of_steatohepatitis"),("TNFa","TNFR1"),("TRAF2","ASK1"),("TRAF2","IKKb"),("XBP1","C_EBPa"),("XBP1","P_De_novo_fatty_acid_synthesis"),("adipoR","AMPK"),("eIF2a","ATF4"),("p38","PPAR_a")]
INHS = [("Akt","GSK_3"),("JNK1_2","IRS_1_2"),("SOCS3","IRS_1_2")]

SILYMARIN_TARGETS = ['CASP3', 'CASP7', 'CASP8', 'CYP2E1', 'IL_8', 'TNFa', 'NF_kB', 'TGF_b1']

ALL_P_OUTPUTS = [
    'P_Hyperinsulinemia', 'P_De_novo_fatty_acid_synthesis', 'P_Improvement_of_NAFLD',
    'P_Development_of_NAFLD', 'P_Adipogenesis', 'P_HCC_proliferation',
    'P_Apoptosis', 'P_Development_of_steatohepatitis', 'P_Cell_death',
    'P_Inflammation', 'P_Fibrosis', 'P_Hepatocyte_injury',
]

plt.rcParams.update({
    'font.family': 'sans-serif', 'font.sans-serif': ['Arial', 'DejaVu Sans'],
    'font.size': 12, 'axes.labelsize': 13, 'axes.titlesize': 14,
    'xtick.labelsize': 11, 'ytick.labelsize': 11, 'savefig.dpi': 300,
})


def build_graph():
    """建立有向圖（ACTS + INHS 都視為有向邊，因為 Silymarin 是透過抑制
    這些基因的活化速率來產生下游效應，不論 edge 是 activation 或
    inhibition，傳遞方向都是 source -> target）。"""
    graph = defaultdict(set)
    for s, t in ACTS + INHS:
        graph[s].add(t)
    return graph


def bfs_shortest_distance(graph, sources, target):
    """從多個 source node 之一出發，BFS 找到 target 的最短距離（hop 數）。"""
    if target in sources:
        return 0
    visited = set(sources)
    queue = deque((s, 0) for s in sources)
    while queue:
        node, dist = queue.popleft()
        for nxt in graph.get(node, []):
            if nxt == target:
                return dist + 1
            if nxt not in visited:
                visited.add(nxt)
                queue.append((nxt, dist + 1))
    return np.nan  # 無法到達


def main():
    print('[File check]')
    ok = os.path.exists(IN_XLSX)
    print(f'  {"OK  " if ok else "MISS"} {IN_XLSX}')
    if not ok:
        raise FileNotFoundError('請先跑 quantify_all_outputs_separation_vs_intervention.py')

    graph = build_graph()

    print('\n[Computing shortest topological distance: Silymarin targets -> each P_*]')
    distances = {}
    for pw in ALL_P_OUTPUTS:
        d = bfs_shortest_distance(graph, SILYMARIN_TARGETS, pw)
        distances[pw] = d
        print(f'  {pw:<32} distance = {d}')

    df = pd.read_excel(IN_XLSX, sheet_name='AllOutputs_72combos')
    df['topo_distance'] = df['P_output'].map(distances)

    # 排除已知反向節點與已知例外（沿用前一支腳本的判斷）
    df_main = df[df['P_output'] != 'P_Improvement_of_NAFLD'].copy()
    df_main = df_main.dropna(subset=['NAFL/NASH_ratio', 'topo_distance'])

    print(f'\n[Distance distribution among 11 main outputs]')
    dist_summary = df_main.groupby('P_output')['topo_distance'].first().sort_values()
    print(dist_summary.to_string())

    # ---------- 距離 vs |效果強度| 相關分析 ----------
    # 用 |ratio - 1| 當作「效果強度」（不管方向，純看效果有多明顯偏離1）
    df_main['effect_strength'] = (df_main['NAFL/NASH_ratio'] - 1).abs()

    valid = df_main.dropna(subset=['effect_strength', 'topo_distance'])
    x = valid['topo_distance'].values
    y = valid['effect_strength'].values
    sr, sp = stats.spearmanr(x, y)
    pr, pp = stats.pearsonr(x, y)

    print('\n' + '=' * 60)
    print('[Topological distance vs |NAFL/NASH ratio - 1| (effect strength)]')
    print('=' * 60)
    print(f'  N = {len(valid)}')
    print(f'  Pearson  r = {pr:.3f}, p = {pp:.4f}')
    print(f'  Spearman r = {sr:.3f}, p = {sp:.4f}')
    print('  (負相關 = 距離越近效果越強，支持拓樸距離稀釋假說)')

    # 按距離分組，看每組的 median |effect|
    print('\n  各距離分組的中位數效果強度:')
    for d in sorted(valid['topo_distance'].unique()):
        sub = valid[valid['topo_distance'] == d]
        print(f'    distance={int(d)}: N={len(sub)}, median |ratio-1| = {sub["effect_strength"].median():.3f}, '
              f'outputs = {sorted(sub["P_output"].unique())}')

    # ---------- 存 Excel ----------
    with pd.ExcelWriter(OUT_XLSX, engine='openpyxl') as w:
        df_main.to_excel(w, sheet_name='Distance_vs_Effect', index=False)
        dist_summary.reset_index().rename(columns={0: 'topo_distance'}).to_excel(
            w, sheet_name='DistanceByOutput', index=False)
        pd.DataFrame([
            ('N', len(valid)),
            ('Pearson r (distance vs |ratio-1|)', round(pr, 3)),
            ('Pearson p', round(pp, 4)),
            ('Spearman r', round(sr, 3)),
            ('Spearman p', round(sp, 4)),
        ], columns=['Statistic', 'Value']).to_excel(w, sheet_name='Summary', index=False)
    print(f'\n  Saved: {OUT_XLSX}')

    # ---------- 圖 ----------
    fig, ax = plt.subplots(figsize=(8, 6.5))
    cmap = plt.get_cmap('tab20')
    outputs_sorted = dist_summary.index.tolist()
    colors = {pw: cmap(i / len(outputs_sorted)) for i, pw in enumerate(outputs_sorted)}
    for pw in outputs_sorted:
        sub = valid[valid['P_output'] == pw]
        jitter = np.random.uniform(-0.08, 0.08, len(sub))
        ax.scatter(sub['topo_distance'] + jitter, sub['effect_strength'],
                   color=colors[pw], s=70, label=f'{pw.replace("P_","")} (d={int(dist_summary[pw])})',
                   edgecolor='white', linewidth=0.6)

    xs = np.linspace(x.min(), x.max(), 50)
    slope, intercept = np.polyfit(x, y, 1)
    ax.plot(xs, slope * xs + intercept, color='gray', ls='--', lw=1.5,
            label=f'Linear fit (r={pr:.2f}, p={pp:.3f})')

    ax.set_xlabel('Shortest topological distance: Silymarin targets -> P_* (hops)')
    ax.set_ylabel('|Silymarin AUC ratio (NAFL/NASH) - 1|\n(effect strength, direction-agnostic)')
    ax.set_title('Does topological distance explain signal dilution?\n'
                  f'Spearman r={sr:.2f} p={sp:.3f}')
    ax.legend(fontsize=8, loc='upper left', bbox_to_anchor=(1.01, 1.0))
    ax.grid(alpha=0.25)
    fig.tight_layout()
    fig.savefig(OUT_FIG, dpi=300, bbox_inches='tight')
    plt.close(fig)
    print(f'  Saved: {OUT_FIG}')

    print('\nDone.')


if __name__ == '__main__':
    main()
