#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
step5_model_output_validation.py  (GSE48452)  [adaptive dynamic indicators]
================================================================
TRUE model-level validation: feed GSE48452 expression into the existing v10
ODE model and check whether disease-progression conclusions are reproduced.

METHODOLOGY (important):
  P_* outputs are modeled as conserved active/inactive pools (initial_level=100,
  following Tseng 2024). Outputs with several upstream activating inputs saturate
  near 100 even at the Normal state, so their ABSOLUTE steady-state level cannot
  discriminate disease stage. However, disease severity is encoded in the RATE of
  approach to saturation (see Fig 5/6: more severe states saturate faster).

  We therefore use ADAPTIVE indicators:
    - Saturating outputs (Normal-driven steady state > 50): time-to-half-maximal
      activation (t-half). Earlier t-half = more severe.
    - Non-saturating outputs (Normal-driven steady state <= 50): steady-state
      level. Higher = more severe.

  This mirrors how Fig 5/6 are read (dynamics, not terminal absolute values) and
  lets all 8 disease outputs be evaluated without exclusion.

Cross-platform note:
  GSE126848 ratios came from pyDESeq2 (2^log2FC, shrinkage); GSE48452 is
  microarray (no shrinkage). Validation is QUALITATIVE (disease-progression
  direction/dynamics), not quantitative replication of absolute P_* ratios.
  Absolute P_* fold-changes are NOT compared to literature here; the Discussion
  separately bounds the model's quantitative limits (mRNA-protein discordance;
  pool saturation) using literature anchors (Feldstein 2003: hepatocyte apoptosis
  ~5x NASH/Normal; TNF/IL mRNA ~15x).

GSE48452 specifics:
  groups Control/Healthy obese/Steatosis/Nash; 'after surgery' excluded via
  'usable'. NAFL=Steatosis, NASH=Nash, reference=Control.

Inputs:
  GSE48452_expression_probe.csv, GSE48452_sample_metadata.csv,
  GPL11532_probe_to_gene.csv, node_name_table_hsa04932.xlsx
Model files (copy in, or set MODEL_DIR):
  ode_model_pydeseq2_NASH_vs_Normal_v10_mean.py (+ NAFL version for reference)

Outputs:
  GSE48452_node_initial_ratios.xlsx
  GSE48452_model_dynamics.xlsx          (t-half + steady state + early AUC per output)
  Validation_model_progression.xlsx     (adaptive-indicator progression table)
  Fig_model_validation_dynamics.png      (P_* trajectories: Normal/NAFL/NASH)
================================================================
"""
import os
import importlib.util
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.integrate import odeint

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
MODEL_DIR = SCRIPT_DIR
MODEL_TEMPLATE = os.path.join(MODEL_DIR, 'ode_model_pydeseq2_NASH_vs_Normal_v10_mean.py')

EXPR_CSV   = os.path.join(SCRIPT_DIR, 'GSE48452_expression_probe.csv')
META_CSV   = os.path.join(SCRIPT_DIR, 'GSE48452_sample_metadata.csv')
ANNOT_CSV  = os.path.join(SCRIPT_DIR, 'GPL11532_probe_to_gene.csv')
NODE_TABLE = os.path.join(SCRIPT_DIR, 'node_name_table_hsa04932.xlsx')

OUT_RATIOS = os.path.join(SCRIPT_DIR, 'GSE48452_node_initial_ratios.xlsx')
OUT_DYN    = os.path.join(SCRIPT_DIR, 'GSE48452_model_dynamics.xlsx')
OUT_PROG   = os.path.join(SCRIPT_DIR, 'Validation_model_progression.xlsx')
OUT_FIG    = os.path.join(SCRIPT_DIR, 'Fig_model_validation_dynamics.png')

DISEASE_GROUPS = {'Obese': 'Healthy obese', 'NAFL': 'Steatosis', 'NASH': 'Nash'}
REFERENCE_GROUP = 'Control'
USE_USABLE_FLAG = True
# Disease-severity ordering used for progression checks and plotting
STAGE_ORDER = ['Normal', 'Obese', 'NAFL', 'NASH']

T_MAX, N_POINTS = 300.0, 3001
HALF_TARGET = 50.0        # half-maximal activation level (pool max = 100)
SATURATION_BASELINE = 50.0  # Normal-driven steady state above this = saturating
AUC_CUT = 75.0

DISEASE_OUTPUTS = [
    'P_Hepatocyte_injury', 'P_Cell_death', 'P_Inflammation',
    'P_Fibrosis', 'P_Apoptosis', 'P_HCC_proliferation',
    'P_Development_of_steatohepatitis', 'P_Development_of_NAFLD',
]


def load_model(path):
    name = 'm_' + os.path.basename(path).replace('.py', '').replace('-', '_')
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec); spec.loader.exec_module(mod)
    return mod


def load_node_gene_map():
    df = pd.read_excel(NODE_TABLE, 'NodeNameTable')
    node_map = {}
    for _, r in df.iterrows():
        genes = []
        for g in str(r['gene_symbols']).split(','):
            g = g.strip().rstrip('*').strip()
            if '(' in g:
                g = g.split('(')[0].strip()
            if g and g.lower() != 'nan' and g not in genes:
                genes.append(g)
        node_map[r['ode_variable_name']] = genes
    return node_map


def build_gene_group_means():
    expr = pd.read_csv(EXPR_CSV, index_col=0); expr.index = expr.index.astype(str)
    meta = pd.read_csv(META_CSV)
    if USE_USABLE_FLAG and 'usable' in meta.columns:
        meta = meta[meta['usable']]
    annot = pd.read_csv(ANNOT_CSV); annot['probe_id'] = annot['probe_id'].astype(str)
    e2 = expr.merge(annot, left_index=True, right_on='probe_id', how='inner')
    ge = e2.groupby('gene_symbol')[list(expr.columns)].mean()
    gm = {}
    for grp in meta['group'].unique():
        s = [x for x in meta.loc[meta['group'] == grp, 'GSM'] if x in ge.columns]
        gm[grp] = ge[s].mean(axis=1)
    return pd.DataFrame(gm)


def compute_node_ratios(model, node_gene_map, gm, disease_grp, ref_grp):
    rows = []
    for node in model.GENE_METADATA.keys():
        genes = node_gene_map.get(node, [])
        ratios, found = [], []
        for g in genes:
            if g in gm.index and disease_grp in gm.columns and ref_grp in gm.columns:
                ratios.append(float(2 ** (gm.loc[g, disease_grp] - gm.loc[g, ref_grp])))
                found.append(g)
        rows.append({'node': node, 'n_genes_total': len(genes),
                     'n_genes_found': len(found), 'genes_found': ','.join(found),
                     'ratio': float(np.mean(ratios)) if ratios else 1.0,
                     'mapped': bool(ratios)})
    return pd.DataFrame(rows)


def simulate_trajectories(model, ratio_map):
    y0 = np.zeros(len(model.STATE_VARS))
    for i, sv in enumerate(model.STATE_VARS):
        if sv.endswith('_inactive'):
            node = sv[:-len('_inactive')]
            if node in ratio_map:
                y0[i] = ratio_map[node]
            elif node in model.PATHWAY_METADATA:
                y0[i] = model.PATHWAY_METADATA[node].get('initial_level', 100.0)
            else:
                y0[i] = 1.0
    t = np.linspace(0, T_MAX, N_POINTS)
    y = odeint(model.ode_system, y0, t)
    return t, y


def t_half(t, traj, target=HALF_TARGET):
    above = np.where(traj >= target)[0]
    return float(t[above[0]]) if len(above) > 0 else np.inf


def steady(traj):
    return float(np.mean(traj[-int(len(traj) * 0.1):]))


def early_auc(t, traj, t_cut=AUC_CUT):
    mask = t <= t_cut
    return float(np.trapezoid(traj[mask], t[mask]))


def main():
    print('=' * 64)
    print('GSE48452 model-level validation (adaptive dynamic indicators)')
    print('=' * 64)
    for f in (MODEL_TEMPLATE, EXPR_CSV, META_CSV, ANNOT_CSV, NODE_TABLE):
        if not os.path.exists(f):
            raise FileNotFoundError(f'Required file not found: {f}')

    model = load_model(MODEL_TEMPLATE)
    p_idx = {pw: model.STATE_VARS.index(pw + '_active') for pw in model.PATHWAY_METADATA}
    node_gene_map = load_node_gene_map()
    n_multi = sum(1 for v in node_gene_map.values() if len(v) > 1)
    print(f'Model: {len(model.GENE_METADATA)} gene nodes, '
          f'{len(model.PATHWAY_METADATA)} P_* outputs '
          f'({n_multi} multi-gene nodes averaged)')

    gm = build_gene_group_means()
    print(f'GSE48452 gene means built. Groups: {list(gm.columns)}')

    # Node ratios per stage (+ Normal = all gene ratios 1.0)
    ratios = {'Normal': {n: 1.0 for n in model.GENE_METADATA}}
    ratio_dfs = {}
    for stage, grp in DISEASE_GROUPS.items():
        df = compute_node_ratios(model, node_gene_map, gm, grp, REFERENCE_GROUP)
        ratio_dfs[stage] = df
        ratios[stage] = dict(zip(df['node'], df['ratio']))
        full = (df['n_genes_found'] == df['n_genes_total']).sum()
        print(f'  {stage} ({grp} vs {REFERENCE_GROUP}): {full}/{len(df)} nodes full-mapped')

    with pd.ExcelWriter(OUT_RATIOS, engine='openpyxl') as xw:
        for stage, df in ratio_dfs.items():
            df.to_excel(xw, sheet_name=f'{stage}_ratios', index=False)
    print(f'Node ratios saved: {OUT_RATIOS}')

    # Stages actually available (Normal always; others if their group exists)
    sim_stages = [s for s in STAGE_ORDER if s == 'Normal' or s in ratios]

    # Simulate trajectories
    print(f'\nSimulating trajectories ({" / ".join(sim_stages)})...')
    traj = {}
    for stage in sim_stages:
        t, y = simulate_trajectories(model, ratios[stage])
        traj[stage] = (t, y)
    print('  done.')

    # Compute indicators per output
    rows = []
    for pw in model.PATHWAY_METADATA:
        i = p_idx[pw]
        base = steady(traj['Normal'][1][:, i])
        saturating = base > SATURATION_BASELINE
        rec = {'P_output': pw, 'Normal_baseline': base,
               'type': 'saturating' if saturating else 'non-saturating'}
        for stage in sim_stages:
            t, y = traj[stage]
            rec[f'{stage}_steady'] = steady(y[:, i])
            rec[f'{stage}_thalf'] = t_half(t, y[:, i])
            rec[f'{stage}_AUC75'] = early_auc(t, y[:, i])
        rows.append(rec)
    dyn_df = pd.DataFrame(rows)
    dyn_df.to_excel(OUT_DYN, index=False)
    print(f'Dynamics table saved: {OUT_DYN}')

    # ----------------------------------------------------------------
    # OBJECTIVE REPORTING (no pass/fail threshold imposed here).
    # We list, for every disease output, the adaptive indicator value at
    # each simulated stage, plus two reference directional summaries:
    #   - NASH vs NAFL  (primary disease-progression contrast)
    #   - full gradient Normal->...->NASH monotonicity (info only)
    # Interpretation is left to manual review of these values + the figure.
    # ----------------------------------------------------------------
    grad_stages = [s for s in sim_stages]   # severity order
    print('\n' + '=' * 78)
    print('Model output by stage (objective; adaptive indicator per output)')
    print('  saturating  -> t-half  (smaller = reaches saturation earlier)')
    print('  non-satur.  -> steady  (larger  = higher activation)')
    print('=' * 78)
    hdr = f'{"P_output":<33}{"type":>9}{"ind":>8}'
    for s in grad_stages:
        hdr += f'{s:>8}'
    hdr += f'{"NASHvNAFL":>11}{"fullMono":>9}'
    print(hdr)
    prog_rows = []
    for pw in DISEASE_OUTPUTS:
        row = dyn_df[dyn_df.P_output == pw].iloc[0]
        if row['type'] == 'saturating':
            vals = [row[f'{s}_thalf'] for s in grad_stages]
            ind = 't-half'
            disp = ['inf' if not np.isfinite(v) else f'{v:.0f}' for v in vals]
            # NASH faster (smaller t-half) than NAFL = progression direction
            v_nafl = row.get('NAFL_thalf', np.nan)
            v_nash = row.get('NASH_thalf', np.nan)
            nash_vs_nafl = ('NASH faster' if (np.isfinite(v_nash) and
                            np.isfinite(v_nafl) and v_nash < v_nafl)
                            else ('equal' if v_nash == v_nafl else 'NAFL faster'))
            full_mono = all(vals[k] >= vals[k+1] for k in range(len(vals)-1)
                            if np.isfinite(vals[k]) and np.isfinite(vals[k+1]))
        else:
            vals = [row[f'{s}_steady'] for s in grad_stages]
            ind = 'steady'
            disp = [f'{v:.2f}' for v in vals]
            v_nafl = row.get('NAFL_steady', np.nan)
            v_nash = row.get('NASH_steady', np.nan)
            nash_vs_nafl = ('NASH higher' if v_nash > v_nafl
                            else ('equal' if v_nash == v_nafl else 'NAFL higher'))
            full_mono = all(vals[k] <= vals[k+1] for k in range(len(vals)-1))
        line = f'{pw:<33}{row["type"][:8]:>9}{ind:>8}'
        for d in disp:
            line += f'{d:>8}'
        line += f'{nash_vs_nafl:>11}{("mono" if full_mono else "-"):>9}'
        print(line)
        rec = {'P_output': pw, 'type': row['type'], 'indicator': ind,
               'NASH_vs_NAFL': nash_vs_nafl, 'full_gradient_monotone': full_mono}
        for s, v in zip(grad_stages, vals):
            rec[s] = v
        prog_rows.append(rec)
    print('\n  (Values reported objectively; no pass/fail threshold imposed.')
    print('   Review the figure + table to judge progression.)')
    pd.DataFrame(prog_rows).to_excel(OUT_PROG, index=False)
    print(f'  Stage-value table saved: {OUT_PROG}')

    # Figure: trajectories of the 6 key disease outputs (matches Fig 5: 4 lines)
    key_outs = ['P_Hepatocyte_injury', 'P_Cell_death', 'P_Inflammation',
                'P_Fibrosis', 'P_Apoptosis', 'P_Development_of_steatohepatitis']
    key_outs = [p for p in key_outs if p in model.PATHWAY_METADATA]
    fig, axes = plt.subplots(2, 3, figsize=(14, 8), constrained_layout=True)
    colors = {'Normal': '#1f3a5f', 'Obese': '#4a90d9',
              'NAFL': '#e09f3e', 'NASH': '#c1121f'}
    for ax, pw in zip(axes.flat, key_outs):
        i = p_idx[pw]
        for stage in sim_stages:
            t, y = traj[stage]
            ax.plot(t, y[:, i], label=stage,
                    color=colors.get(stage, 'gray'), lw=2)
        ax.set_title(pw.replace('P_', ''), fontsize=10)
        ax.set_xlabel('Time (h)', fontsize=8)
        ax.set_ylabel('P_active (a.u.)', fontsize=8)
        ax.legend(fontsize=7)
        ax.grid(alpha=0.3)
    plt.suptitle('GSE48452-driven model: disease-stage P_* dynamics',
                 fontsize=12, fontweight='bold')
    plt.savefig(OUT_FIG, dpi=300, bbox_inches='tight')
    plt.close()
    print(f'\nFigure saved: {OUT_FIG}')
    print('\nDone.')


if __name__ == '__main__':
    main()
