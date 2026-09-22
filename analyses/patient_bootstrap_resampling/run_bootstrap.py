#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
R1-2 Bootstrap / patient-level robustness (standalone, for local execution).

For each bootstrap iteration:
  1. Within each of the 4 conditions (Normal/Obese/NAFL/NASH), resample patients
     WITH replacement (same group sizes as original: 14/12/15/16).
  2. Re-run PyDESeq2 (single ~condition design, no batch correction) on the
     resampled count matrix -> gene-level log2FC for NASH_vs_Normal and NAFL_vs_Normal.
  3. Aggregate gene-level ratios (2^log2FC) to the 58 ODE nodes via arithmetic
     mean (v10-Mean), using NODE_GENES mapping transcribed from Supplementary Table S1.
  4. Build Y0 for NASH and NAFL baselines, run ODE (no-drug and ki=0.3 Silymarin)
     using the primary NASH model file (STATE_VARS / ode_system, with Silymarin
     ki=0.3 inhibition dynamically injected into the source at 8 target nodes).
  5. Compute AUC-based % reduction for the 3 core outputs and the NAFL/NASH
     early-intervention advantage ratio.

Resumable: results appended to bootstrap_results.jsonl; already-completed
iteration indices are skipped on restart. Fixed master seed -> per-iteration
seed = BASE_SEED + iteration, so a re-run reproduces identical resamples.

Usage:
    python run_bootstrap.py [N_TOTAL] [MAX_THIS_CALL]

    N_TOTAL       total number of bootstrap iterations wanted (default 200)
    MAX_THIS_CALL if given, stop after running this many NEW iterations in
                  this call (useful to checkpoint progress); omit to run
                  until N_TOTAL is reached in one go.

Requires (place in ./data/ next to this script):
    GSE126848_Gene_counts_raw.txt      (raw ENSG counts, 57 samples x ~19786 genes)
    GSE126848_Gene_counts_mapped.xlsx  (Sheet1: sample_id/condition key row + Gene_Symbol column)

Requires (already included in ./models/):
    ode_model_pydeseq2_NASH_vs_Normal_v10_mean.py

Python deps: numpy, pandas, scipy, pydeseq2, openpyxl
"""
import os, sys, json, time, importlib.util, types
import numpy as np
import pandas as pd
from scipy.integrate import odeint

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from node_gene_map import NODE_GENES

# ---------------- paths (all relative to this script's folder) ----------------
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
OUT_JSONL = os.path.join(HERE, 'bootstrap_results.jsonl')

N_BOOTSTRAP = 200
BASE_SEED = 20260908
T_MAX, N_POINTS = 300.0, 3001
t = np.linspace(0, T_MAX, N_POINTS)
CORE_OUTPUTS = ['P_Cell_death', 'P_Hepatocyte_injury', 'P_Inflammation']

SILYMARIN_TARGETS = [
    ('CASP3', 18, 19), ('CASP7', 20, 21), ('CASP8', 22, 23),
    ('CYP2E1', 26, 27), ('IL_8', 58, 59), ('TNFa', 104, 105),
    ('NF_kB', 80, 81), ('TGF_b1', 100, 101),
]


# ---------------- data prep (done once) ----------------
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
    sym_counts = raw.groupby('Symbol').sum()
    return sym_counts  # genes(symbol) x 57 samples, unresampled


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
    df = pd.DataFrame({'sample_id': ids, 'condition': conds})
    return df


def load_model(path):
    spec = importlib.util.spec_from_file_location(
        os.path.splitext(os.path.basename(path))[0], path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def build_sily_source(src, ki):
    lines = src.split('\n')
    for _, ii, ai in SILYMARIN_TARGETS:
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
                    depth = 0
                    ep = None
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


def make_y0_from_node_ratios(model, node_ratio_map):
    y0 = np.zeros(len(model.STATE_VARS))
    for i, sv in enumerate(model.STATE_VARS):
        if sv.endswith('_inactive'):
            node = sv[:-len('_inactive')]
            if node == 'INS':
                y0[i] = 1.0
            elif node in node_ratio_map and np.isfinite(node_ratio_map[node]):
                y0[i] = node_ratio_map[node]
            elif node in model.PATHWAY_METADATA:
                y0[i] = model.PATHWAY_METADATA[node].get('initial_level', 100.)
            else:
                y0[i] = 1.0
    return y0


def run_sim(model, y0):
    return odeint(model.ode_system, y0, t, args=(None,), mxstep=10000)


def auc_of(traj):
    _trapz = getattr(np, 'trapezoid', None) or np.trapz
    return float(_trapz(traj, t))


def auc_pct_reduction(traj_drug, traj_nodrug):
    a0 = auc_of(traj_nodrug)
    a1 = auc_of(traj_drug)
    return (a0 - a1) / a0 * 100 if a0 > 0 else np.nan


def resample_counts(sym_counts, cond_df, rng):
    """Stratified with-replacement resample of patients; returns a genes x N_resampled
    count matrix with unique synthetic column labels."""
    picked_cols = []
    for cond in ['Normal', 'Obese', 'NAFL', 'NASH']:
        pool = cond_df.loc[cond_df['condition'] == cond, 'sample_id'].values
        picks = rng.choice(pool, size=len(pool), replace=True)
        for k, sid in enumerate(picks):
            picked_cols.append((f'{cond}_{k}', sid, cond))
    new_cols = [p[0] for p in picked_cols]
    src_ids = [p[1] for p in picked_cols]
    conditions = [p[2] for p in picked_cols]
    mat = sym_counts[src_ids].copy()
    mat.columns = new_cols
    meta = pd.DataFrame({'condition': conditions}, index=new_cols)
    return mat, meta


def run_deseq2_both_contrasts(counts_df, meta_df):
    """counts_df: genes x samples (raw counts); meta_df: samples x [condition].
    Fits DeseqDataSet ONCE, then extracts both NASH-vs-Normal and NAFL-vs-Normal
    contrasts from the same fitted dispersions (much faster than fitting twice).
    Returns dict: {'NASH': ratio_series, 'NAFL': ratio_series}."""
    from pydeseq2.dds import DeseqDataSet
    from pydeseq2.ds import DeseqStats

    counts_t = counts_df.T  # samples x genes, required by pydeseq2
    keep = counts_t.sum(axis=0) > 0
    counts_t = counts_t.loc[:, keep]

    dds = DeseqDataSet(
        counts=counts_t.astype(int),
        metadata=meta_df,
        design="~condition",
        refit_cooks=False,
        quiet=True,
    )
    dds.deseq2()

    out = {}
    for level in ['NASH', 'NAFL']:
        stat = DeseqStats(dds, contrast=["condition", level, "Normal"], quiet=True)
        stat.summary()
        out[level] = 2.0 ** stat.results_df['log2FoldChange']
    return out


def aggregate_nodes(ratio_series):
    node_ratio = {}
    for node, genes in NODE_GENES.items():
        vals = [ratio_series[g] for g in genes if g in ratio_series.index and np.isfinite(ratio_series[g])]
        node_ratio[node] = float(np.mean(vals)) if vals else np.nan
    return node_ratio


def already_done():
    done = set()
    if os.path.exists(OUT_JSONL):
        with open(OUT_JSONL) as fh:
            for line in fh:
                try:
                    rec = json.loads(line)
                    done.add(rec['iter'])
                except Exception:
                    pass
    return done


def main(n_iter=N_BOOTSTRAP, max_this_run=None):
    for req, label in [(RAW_COUNTS, 'raw counts'), (MAPPED_XLSX, 'mapped xlsx')]:
        if not os.path.exists(req):
            print(f'[ERROR] Missing {label} file: {req}')
            print('        Place your GSE126848 raw-count and mapped-xlsx files in the data/ folder.')
            sys.exit(1)

    print('[Setup] loading symbol-level counts and condition map...')
    sym_counts = load_symbol_counts()
    cond_df = load_condition_map()
    print(f'  sym_counts shape: {sym_counts.shape}')
    print(f'  condition counts:\n{cond_df["condition"].value_counts()}')

    print('[Setup] loading NASH model (base ODE structure + Silymarin ki=0.3 variant)...')
    m_nash = load_model(MODEL_NASH_PATH)
    p_idx = {pw: m_nash.STATE_VARS.index(pw + '_active') for pw in CORE_OUTPUTS}
    with open(MODEL_NASH_PATH, encoding='utf-8') as fh:
        base_src = fh.read()
    sily_mod = make_mod_from_source(build_sily_source(base_src, 0.3), 'sily_03')

    done = already_done()
    print(f'[Resume] {len(done)} iterations already completed.')

    todo = [i for i in range(n_iter) if i not in done]
    if max_this_run is not None:
        todo = todo[:max_this_run]
    print(f'[Plan] running {len(todo)} iterations this call.')

    fh_out = open(OUT_JSONL, 'a')
    for it in todo:
        t0 = time.time()
        rng = np.random.default_rng(BASE_SEED + it)
        counts_rs, meta_rs = resample_counts(sym_counts, cond_df, rng)

        try:
            ratios = run_deseq2_both_contrasts(counts_rs, meta_rs)
            ratio_nash, ratio_nafl = ratios['NASH'], ratios['NAFL']
        except Exception as e:
            rec = {'iter': it, 'seed': BASE_SEED + it, 'error': str(e)}
            fh_out.write(json.dumps(rec) + '\n')
            fh_out.flush()
            print(f'  iter {it}: FAILED ({e}) [{time.time()-t0:.1f}s]')
            continue

        node_ratio_nash = aggregate_nodes(ratio_nash)
        node_ratio_nafl = aggregate_nodes(ratio_nafl)

        y0_nash = make_y0_from_node_ratios(m_nash, node_ratio_nash)
        y0_nafl = make_y0_from_node_ratios(m_nash, node_ratio_nafl)

        traj_nash_nodrug = run_sim(m_nash, y0_nash)
        traj_nash_drug = run_sim(sily_mod, y0_nash)
        traj_nafl_nodrug = run_sim(m_nash, y0_nafl)
        traj_nafl_drug = run_sim(sily_mod, y0_nafl)

        rec = {'iter': it, 'seed': BASE_SEED + it}
        for pw in CORE_OUTPUTS:
            i = p_idx[pw]
            nash_red = auc_pct_reduction(traj_nash_drug[:, i], traj_nash_nodrug[:, i])
            nafl_red = auc_pct_reduction(traj_nafl_drug[:, i], traj_nafl_nodrug[:, i])
            ratio = nafl_red / nash_red if (np.isfinite(nafl_red) and np.isfinite(nash_red) and nash_red > 0) else None
            rec[f'{pw}_NASH_red'] = round(nash_red, 4) if np.isfinite(nash_red) else None
            rec[f'{pw}_NAFL_red'] = round(nafl_red, 4) if np.isfinite(nafl_red) else None
            rec[f'{pw}_ratio'] = round(ratio, 4) if ratio is not None else None

        fh_out.write(json.dumps(rec) + '\n')
        fh_out.flush()
        print(f'  iter {it}: done [{time.time()-t0:.1f}s]  '
              f'ratios: ' + ', '.join(f"{pw}={rec[pw+'_ratio']}" for pw in CORE_OUTPUTS))

    fh_out.close()
    print('\n[Done this call]')


if __name__ == '__main__':
    n = int(sys.argv[1]) if len(sys.argv) > 1 else N_BOOTSTRAP
    mx = int(sys.argv[2]) if len(sys.argv) > 2 else None
    main(n_iter=n, max_this_run=mx)
