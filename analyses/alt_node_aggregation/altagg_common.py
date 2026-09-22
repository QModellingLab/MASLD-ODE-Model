#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Shared utilities for alt_node_aggregation (R1-3).

Reviewer 1, Major Comment 3:
  "Multi-gene nodes ... are represented by unweighted arithmetic means,
   although their constituent genes have different functions and
   stoichiometry ... The authors should ... test whether alternative
   aggregation rules affect the results."

This reuses the SAME real per-gene DESeq2 pipeline already built and
validated for patient_bootstrap_resampling (single fit here, NOT resampled
-- this recovers the manuscript's actual point-estimate gene-level log2FC
values), then re-aggregates the 16 multi-gene nodes under 3 alternative
rules and compares against the manuscript's default (arithmetic mean of
linear-space ratios, i.e. "v10 Mean").

Requires the SAME two data files already used for patient_bootstrap_
resampling, placed in ./data/:
    GSE126848_Gene_counts_raw.txt
    GSE126848_Gene_counts_mapped.xlsx
"""
import os, sys, types, importlib.util
import numpy as np
import pandas as pd
from scipy.integrate import odeint

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from node_gene_map import NODE_GENES

HERE = os.path.dirname(os.path.abspath(__file__))
# Shared GEO inputs live at <repo>/data/GSE126848/ (single copy for all analyses).
_SHARED = os.path.join(os.path.dirname(os.path.dirname(HERE)), 'data', 'GSE126848')
_LOCAL  = os.path.join(HERE, 'data')
def _find(name):
    for base in (_LOCAL, _SHARED):
        cand = os.path.join(base, name)
        if os.path.exists(cand):
            return cand
    return os.path.join(_SHARED, name)
RAW_COUNTS = _find('GSE126848_Gene_counts_raw.txt')
MAPPED_XLSX = _find('GSE126848_Gene_counts_mapped.xlsx')
MODEL_NASH_PATH = os.path.join(HERE, 'models', 'ode_model_pydeseq2_NASH_vs_Normal_v10_mean.py')
MODEL_NAFL_PATH = os.path.join(HERE, 'models', 'ode_model_pydeseq2_NAFL_vs_Normal_v10_mean.py')
OUT_DIR = os.path.join(HERE, 'outputs')
os.makedirs(OUT_DIR, exist_ok=True)

T_MAX, N_POINTS = 300.0, 3001
T = np.linspace(0, T_MAX, N_POINTS)
CORE_OUTPUTS = ['P_Cell_death', 'P_Hepatocyte_injury', 'P_Inflammation']
KI = 0.3
SILYMARIN_TARGETS = ['CASP3', 'CASP7', 'CASP8', 'CYP2E1', 'IL_8', 'TNFa', 'NF_kB', 'TGF_b1']

MULTI_GENE_NODES = [n for n, genes in NODE_GENES.items() if len(genes) > 1]

AGGREGATION_METHODS = ['mean', 'geomean', 'median', 'max_abs_log2fc']


# ---------------- real per-gene DESeq2 (single fit, NOT resampled) ----------------
def load_symbol_counts():
    raw = pd.read_csv(RAW_COUNTS, sep='\t', index_col=0)
    raw.columns = [str(c) for c in raw.columns]
    df_map = pd.read_excel(MAPPED_XLSX, sheet_name='Sheet1')
    sym = df_map[['Unnamed: 0', 'Gene_Symbol']].iloc[1:]
    sym.columns = ['ENSG', 'Symbol']
    sym = sym.dropna(subset=['Symbol']).drop_duplicates(subset=['ENSG'])
    ensg2sym = dict(zip(sym['ENSG'], sym['Symbol']))
    raw = raw.loc[raw.index.isin(ensg2sym)]
    raw['Symbol'] = raw.index.map(ensg2sym)
    return raw.groupby('Symbol').sum()


def load_condition_map():
    cols = list(pd.read_excel(MAPPED_XLSX, sheet_name='Sheet1', nrows=1).columns)[4:]
    key_row = pd.read_excel(MAPPED_XLSX, sheet_name='Sheet1').iloc[0]

    def cond(c):
        for k in ['Normal', 'Obese', 'NAFL', 'NASH']:
            if c.startswith(k):
                return k
        raise ValueError(c)

    ids = [str(int(key_row[c])).zfill(4) for c in cols]
    conds = [cond(c) for c in cols]
    return pd.DataFrame({'condition': conds}, index=ids)


def run_deseq2_single():
    """Single (non-resampled) PyDESeq2 fit on the full 57-sample cohort --
    recovers the manuscript's actual point-estimate log2FC values. Returns
    dict {'NASH': log2fc_series, 'NAFL': log2fc_series} indexed by gene
    symbol."""
    from pydeseq2.dds import DeseqDataSet
    from pydeseq2.ds import DeseqStats

    sym_counts = load_symbol_counts()
    meta = load_condition_map()
    counts_t = sym_counts[meta.index].T
    keep = counts_t.sum(axis=0) > 0
    counts_t = counts_t.loc[:, keep]

    dds = DeseqDataSet(counts=counts_t.astype(int), metadata=meta,
                        design="~condition", refit_cooks=True, quiet=False)
    dds.deseq2()

    out = {}
    for level in ['NASH', 'NAFL']:
        stat = DeseqStats(dds, contrast=["condition", level, "Normal"], quiet=False)
        stat.summary()
        out[level] = stat.results_df['log2FoldChange']
    return out


def aggregate_node(log2fc_series, genes, method):
    """Aggregate a multi-gene node's ratio (2^log2FC) under one of the
    AGGREGATION_METHODS. Falls back gracefully if some genes are missing
    from the DEG results (uses whichever genes ARE present)."""
    vals_log2fc = [log2fc_series[g] for g in genes if g in log2fc_series.index
                   and np.isfinite(log2fc_series[g])]
    if not vals_log2fc:
        return np.nan
    vals_log2fc = np.array(vals_log2fc)
    vals_linear = 2.0 ** vals_log2fc

    if method == 'mean':
        return float(np.mean(vals_linear))              # manuscript default (v10 Mean)
    elif method == 'geomean':
        return float(2.0 ** np.mean(vals_log2fc))        # geometric mean in log2FC space
    elif method == 'median':
        return float(np.median(vals_linear))
    elif method == 'max_abs_log2fc':
        idx = np.argmax(np.abs(vals_log2fc))
        return float(vals_linear[idx])                   # most-differentially-expressed gene
    else:
        raise ValueError(method)


def build_node_ratios(log2fc_series, model, method):
    """Returns {node_name: ratio} for ALL 58 molecular nodes: multi-gene
    nodes re-aggregated under `method`; single-gene nodes always use the
    model's own baked-in GENE_METADATA value (identical to manuscript,
    since there is only one gene -- no aggregation choice to make)."""
    node_ratio = {}
    for node in model.GENE_METADATA.keys():
        if node in MULTI_GENE_NODES:
            node_ratio[node] = aggregate_node(log2fc_series, NODE_GENES[node], method)
            if not np.isfinite(node_ratio[node]):
                node_ratio[node] = model.GENE_METADATA[node]['initial_ratio']  # fallback
        else:
            node_ratio[node] = model.GENE_METADATA[node]['initial_ratio']
    return node_ratio


# ---------------- ODE plumbing (same pattern as other control folders) ----------------
def load_model(path):
    spec = importlib.util.spec_from_file_location(
        os.path.splitext(os.path.basename(path))[0] + '_altagg', path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def make_y0_from_node_ratios(model, node_ratio):
    y0 = np.zeros(len(model.STATE_VARS))
    for i, sv in enumerate(model.STATE_VARS):
        if sv.endswith('_inactive'):
            node = sv[:-len('_inactive')]
            if node in node_ratio:
                y0[i] = node_ratio[node]
            elif node in model.PATHWAY_METADATA:
                y0[i] = model.PATHWAY_METADATA[node]['initial_level']
            else:
                raise KeyError(node)
    return y0


def build_ki_source(base_src, model, target_nodes, ki):
    idx_pairs = [(model.STATE_VARS.index(n + '_active'),
                  model.STATE_VARS.index(n + '_inactive')) for n in target_nodes]
    lines = base_src.split('\n')
    for ai, ii in idx_pairs:
        for idx in (ii, ai):
            for li, line in enumerate(lines):
                if line.lstrip().startswith(f'dydt[{idx}]'):
                    eq = line.find('=')
                    rhs = line[eq + 1:].strip()
                    if rhs.startswith('-('):
                        sign, bs = '-', 1
                    elif rhs.startswith('('):
                        sign, bs = '', 0
                    else:
                        continue
                    depth, ep = 0, None
                    for p, ch in enumerate(rhs[bs:], start=bs):
                        if ch == '(':
                            depth += 1
                        elif ch == ')':
                            depth -= 1
                            if depth == 0:
                                ep = p
                                break
                    if ep is None:
                        continue
                    act = rhs[bs:ep + 1]
                    rem = rhs[ep + 1:]
                    ind = line[:len(line) - len(line.lstrip())]
                    lines[li] = f"{ind}dydt[{idx}] = {sign}(({act})*{ki:.4f}){rem}"
                    break
    return '\n'.join(lines)


def make_mod_from_source(src, name):
    mod = types.ModuleType(name)
    mod.__file__ = name
    exec(compile(src, name, 'exec'), mod.__dict__)
    return mod


def run_sim(model, y0):
    return odeint(model.ode_system, y0, T, args=(None,), mxstep=10000)


def pct_reduction(traj_drug_col, traj_nodrug_col):
    f = getattr(np, 'trapezoid', None) or np.trapz
    a0 = float(f(traj_nodrug_col, T))
    a1 = float(f(traj_drug_col, T))
    return (a0 - a1) / a0 * 100 if a0 > 0 else np.nan


def output_active_index(model, pw):
    return model.STATE_VARS.index(pw + '_active')
