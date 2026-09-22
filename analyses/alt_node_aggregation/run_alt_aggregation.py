#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Alternative multi-gene node aggregation rules (R1-3).

Reviewer 1, Major Comment 3:
  "Multi-gene nodes, including large respiratory complexes and AMPK, are
   represented by unweighted arithmetic means, although their constituent
   genes have different functions and stoichiometry ... The authors should
   clarify this mapping and test whether alternative aggregation rules
   affect the results."

This script:
  1. Runs a SINGLE (non-resampled) PyDESeq2 fit on the full GSE126848
     cohort (57 samples) -- recovering the manuscript's actual per-gene
     log2FC point estimates for NASH_vs_Normal and NAFL_vs_Normal.
  2. Re-aggregates the 16 multi-gene nodes under 4 methods:
       mean            -- manuscript default (arithmetic mean of 2^log2FC)
       geomean         -- geometric mean (2^(mean(log2FC)))
       median          -- median of 2^log2FC
       max_abs_log2fc  -- value of whichever constituent gene has the
                          largest |log2FC| ("most differentially expressed
                          representative" -- a principled alternative to
                          KEGG's arbitrary gene-listing order)
  3. Re-runs the manuscript's Table-4-style silymarin AUC analysis (NASH/
     NAFL %reduction, advantage ratio) under each aggregation method and
     compares.

Single-gene nodes are UNCHANGED across all 4 methods (nothing to
aggregate), so only the 16 multi-gene nodes' initial ratios differ between
runs.

Requires (place in ./data/):
    GSE126848_Gene_counts_raw.txt
    GSE126848_Gene_counts_mapped.xlsx
(the SAME two files already used for patient_bootstrap_resampling)

Usage: python run_alt_aggregation.py (no parameters; just press F5)
Runtime: ~1-2 minutes (one DESeq2 fit + 4 methods x 4 ODE solves)
"""
import os
import numpy as np
import pandas as pd

from altagg_common import (
    RAW_COUNTS, MAPPED_XLSX, MODEL_NASH_PATH, MODEL_NAFL_PATH, OUT_DIR,
    CORE_OUTPUTS, KI, SILYMARIN_TARGETS, MULTI_GENE_NODES, AGGREGATION_METHODS,
    run_deseq2_single, aggregate_node, build_node_ratios, load_model,
    make_y0_from_node_ratios, build_ki_source, make_mod_from_source,
    run_sim, pct_reduction, output_active_index,
)
from node_gene_map import NODE_GENES


def main():
    for req, label in [(RAW_COUNTS, 'raw counts'), (MAPPED_XLSX, 'mapped xlsx')]:
        if not os.path.exists(req):
            print(f'[ERROR] Missing {label} file: {req}')
            print('        Place your GSE126848 raw-count and mapped-xlsx files in data/.')
            return

    print('[1/3] Running single (non-resampled) PyDESeq2 fit on the full cohort...')
    log2fc = run_deseq2_single()
    print(f'  Done. {len(log2fc["NASH"])} genes with NASH_vs_Normal log2FC, '
          f'{len(log2fc["NAFL"])} genes with NAFL_vs_Normal log2FC.')

    m_nash = load_model(MODEL_NASH_PATH)
    m_nafl = load_model(MODEL_NAFL_PATH)

    print('\n[2/3] Comparing multi-gene node ratios across aggregation methods '
          '(first 5 multi-gene nodes shown):')
    node_rows = []
    for node in MULTI_GENE_NODES:
        row = {'node': node, 'n_genes': len(NODE_GENES[node]),
               'manuscript_ratio_NASH': m_nash.GENE_METADATA[node]['initial_ratio'],
               'manuscript_ratio_NAFL': m_nafl.GENE_METADATA[node]['initial_ratio']}
        for method in AGGREGATION_METHODS:
            row[f'{method}_NASH'] = aggregate_node(log2fc['NASH'], NODE_GENES[node], method)
            row[f'{method}_NAFL'] = aggregate_node(log2fc['NAFL'], NODE_GENES[node], method)
        node_rows.append(row)
    node_df = pd.DataFrame(node_rows)
    print(node_df[['node', 'n_genes', 'mean_NASH', 'geomean_NASH', 'median_NASH',
                    'max_abs_log2fc_NASH']].head(5).to_string(index=False))

    print('\n[3/3] Re-running the silymarin AUC analysis under each aggregation method...')
    with open(MODEL_NASH_PATH, encoding='utf-8') as f:
        src_nash = f.read()
    with open(MODEL_NAFL_PATH, encoding='utf-8') as f:
        src_nafl = f.read()

    result_rows = []
    print(f'\n{"method":16s} {"output":22s} {"NASH%":>8s} {"NAFL%":>8s} '
          f'{"ratio":>7s} {"diff_pp":>8s}')
    for method in AGGREGATION_METHODS:
        node_ratio_nash = build_node_ratios(log2fc['NASH'], m_nash, method)
        node_ratio_nafl = build_node_ratios(log2fc['NAFL'], m_nafl, method)

        y0_nash = make_y0_from_node_ratios(m_nash, node_ratio_nash)
        y0_nafl = make_y0_from_node_ratios(m_nafl, node_ratio_nafl)

        drug_nash = make_mod_from_source(
            build_ki_source(src_nash, m_nash, SILYMARIN_TARGETS, KI), f'dnash_{method}')
        drug_nafl = make_mod_from_source(
            build_ki_source(src_nafl, m_nafl, SILYMARIN_TARGETS, KI), f'dnafl_{method}')

        tn0, tn1 = run_sim(m_nash, y0_nash), run_sim(drug_nash, y0_nash)
        tf0, tf1 = run_sim(m_nafl, y0_nafl), run_sim(drug_nafl, y0_nafl)

        for pw in CORE_OUTPUTS:
            i = output_active_index(m_nash, pw)
            nash_red = pct_reduction(tn1[:, i], tn0[:, i])
            nafl_red = pct_reduction(tf1[:, i], tf0[:, i])
            ratio = nafl_red / nash_red if nash_red > 0 else np.nan
            result_rows.append({
                'aggregation_method': method,
                'is_manuscript_default': (method == 'mean'),
                'output': pw, 'NASH_red_pct': nash_red, 'NAFL_red_pct': nafl_red,
                'advantage_ratio': ratio, 'diff_pp': nafl_red - nash_red,
            })
            print(f'{method:16s} {pw:22s} {nash_red:8.3f} {nafl_red:8.3f} '
                  f'{ratio:7.3f} {nafl_red-nash_red:8.3f}')

    result_df = pd.DataFrame(result_rows)
    xlsx_path = os.path.join(OUT_DIR, 'alt_node_aggregation_summary.xlsx')
    with pd.ExcelWriter(xlsx_path, engine='openpyxl') as writer:
        node_df.to_excel(writer, sheet_name='node_ratios_by_method', index=False)
        result_df.to_excel(writer, sheet_name='AUC_results_by_method', index=False)
        pivot = result_df.pivot_table(index='aggregation_method', columns='output',
                                       values='advantage_ratio')
        pivot.to_excel(writer, sheet_name='ratio_pivot')
    print(f'\nExcel summary saved: {xlsx_path}')

    print('\n-- Summary: does advantage_ratio > 1 (NAFL > NASH) hold across ALL aggregation methods? --')
    for pw in CORE_OUTPUTS:
        sub = result_df[result_df['output'] == pw]
        all_hold = (sub['advantage_ratio'] > 1).all()
        rng = (sub['advantage_ratio'].min(), sub['advantage_ratio'].max())
        print(f'  {pw:22s} all_hold={all_hold}   ratio range=[{rng[0]:.3f}, {rng[1]:.3f}]')


if __name__ == '__main__':
    main()
