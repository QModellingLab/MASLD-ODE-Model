#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Compare the Morris global-sensitivity ranking (from run_morris_gsa.py) against
simple topology-based rankings, addressing the second half of npj SBA
Reviewer 1, Major Comment 6:
  "Authors should ... compare the sensitivity results with simple
   topology-based rankings."

For each of the 58 molecular nodes, two purely topological metrics are
computed from the REAL hsa04932 interaction network (Supplementary Table
S2, 106 edges; see edges_table_s2.py):

  1. out_degree   -- number of direct downstream targets (a simple "hub-ness"
                      metric; a purely topological analogue of "importance").
  2. distance_to_<output> -- shortest directed path length (in hops) from
                      the node to each of the 3 core pathway-output nodes.
                      Closer nodes (fewer hops) might naively be expected
                      to have larger sensitivity if the ranking were purely
                      a function of network proximity.

Spearman rank correlation (with p-value) is computed between each core
output's Morris mu* ranking and (a) out_degree, (b) inverse distance to
that output, to test whether the global-sensitivity ranking is well
predicted by these simple topological properties alone, or whether it
carries additional information beyond topology (dynamics/kinetics/initial
conditions).

Usage:
    python step2_topology_comparison.py
(run AFTER run_morris_gsa.py has produced outputs/global_sensitivity_summary.xlsx)
"""
import os
import sys
import numpy as np
import pandas as pd
import networkx as nx
from scipy.stats import spearmanr

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from gsa_common import OUT_DIR, CORE_OUTPUTS
from edges_table_s2 import EDGES

GSA_XLSX = os.path.join(OUT_DIR, 'global_sensitivity_summary.xlsx')
OUT_XLSX = os.path.join(OUT_DIR, 'topology_comparison_summary.xlsx')


def build_graph():
    G = nx.DiGraph()
    G.add_edges_from(EDGES)
    return G


def topology_metrics(G, molecular_nodes):
    out_degree = {n: G.out_degree(n) if n in G else 0 for n in molecular_nodes}

    dist = {pw: {} for pw in CORE_OUTPUTS}
    for pw in CORE_OUTPUTS:
        # shortest path length TO the output node, from every node that can reach it
        try:
            rev_lengths = nx.single_target_shortest_path_length(G, pw)
        except nx.NodeNotFound:
            rev_lengths = {}
        for n in molecular_nodes:
            dist[pw][n] = rev_lengths.get(n, np.inf)  # inf = cannot reach this output at all

    return out_degree, dist


def main():
    if not os.path.exists(GSA_XLSX):
        print(f'[ERROR] {GSA_XLSX} not found -- run run_morris_gsa.py first.')
        return

    G = build_graph()
    print(f'Graph: {G.number_of_nodes()} nodes, {G.number_of_edges()} edges '
          f'(should be 106 edges from Table S2)')

    all_gsa = pd.read_excel(GSA_XLSX, sheet_name='all_combined')
    molecular_nodes = sorted(all_gsa['node'].unique())

    out_degree, dist = topology_metrics(G, molecular_nodes)

    summary_rows = []
    for pw in CORE_OUTPUTS:
        df = all_gsa[all_gsa['output'] == pw].copy()
        df['out_degree'] = df['node'].map(out_degree)
        df['dist_to_output'] = df['node'].map(dist[pw])
        # inverse distance (closer = higher "topological importance"); unreachable -> 0
        df['inv_dist_to_output'] = df['dist_to_output'].apply(
            lambda d: 0.0 if not np.isfinite(d) or d == 0 else 1.0 / d)

        finite = df[np.isfinite(df['dist_to_output'])]
        rho_deg, p_deg = spearmanr(df['mu_star'], df['out_degree'])
        rho_dist, p_dist = spearmanr(finite['mu_star'], finite['inv_dist_to_output'])

        print(f'\n=== {pw} ===')
        print(f'  Spearman(mu*, out_degree)        rho={rho_deg:.3f}  p={p_deg:.4f}  (n={len(df)})')
        print(f'  Spearman(mu*, 1/dist_to_output)  rho={rho_dist:.3f}  p={p_dist:.4f}  '
              f'(n={len(finite)}, {len(df)-len(finite)} unreachable excluded)')

        print(f'  Top 10 by mu* (global sensitivity) with topology annotations:')
        top10 = df.sort_values('mu_star', ascending=False).head(10)
        print(top10[['node', 'rank', 'mu_star', 'out_degree', 'dist_to_output']].to_string(index=False))

        summary_rows.append({
            'output': pw, 'spearman_rho_outdegree': rho_deg, 'spearman_p_outdegree': p_deg,
            'spearman_rho_invdist': rho_dist, 'spearman_p_invdist': p_dist,
            'n_total': len(df), 'n_reachable': len(finite),
        })

    stats_df = pd.DataFrame(summary_rows)
    full_df = all_gsa.copy()
    full_df['out_degree'] = full_df.apply(
        lambda r: out_degree.get(r['node'], np.nan), axis=1)
    full_df['dist_to_output'] = full_df.apply(
        lambda r: dist[r['output']].get(r['node'], np.nan), axis=1)

    with pd.ExcelWriter(OUT_XLSX, engine='openpyxl') as writer:
        stats_df.to_excel(writer, sheet_name='spearman_summary', index=False)
        full_df.to_excel(writer, sheet_name='full_annotated', index=False)
    print(f'\nExcel summary saved: {OUT_XLSX}')

    print('\nInterpretation guide:')
    print('  A LOW |rho| and/or high p-value means the Morris global-sensitivity')
    print('  ranking is NOT well explained by simple topology (out-degree or')
    print('  proximity to the output) alone -- i.e. dynamics/kinetics/initial')
    print('  conditions carry real information beyond "this gene is a hub" or')
    print('  "this gene is close to the output in the diagram". A HIGH |rho| with')
    print('  low p-value would support the reviewer\'s concern that the sensitivity')
    print('  ranking is largely predictable from network topology alone.')


if __name__ == '__main__':
    main()
