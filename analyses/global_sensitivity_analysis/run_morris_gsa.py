#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Global sensitivity analysis (Morris elementary-effects screening) for the
three core pathway outputs, at the steatohepatitis (NASH) baseline -- same
condition as the manuscript's local sensitivity analysis (Fig. 6) -- so the
two are directly comparable.

Addresses npj SBA Reviewer 1, Major Comment 6:
  "A global sensitivity analysis is also needed, e.g. the ranking of these
   genes among all genes and p-values to define significance or confidence."

Method: Morris elementary-effects screening (SALib), chosen over
variance-based Sobol indices because Sobol would need ~50-100x more model
evaluations for the same number of parameters (58) to converge, which is
not a good use of compute for a "screening"-level global sensitivity
question. Morris still (a) varies ALL 58 initial-condition parameters
SIMULTANEOUSLY per trajectory (unlike the manuscript's local one-at-a-time
+/-1% perturbation), (b) explores a much wider range per parameter
(0.3x-3.0x the true NASH-baseline ratio, vs +/-1%), and (c) provides a
bootstrap 95% CI on each parameter's importance (mu*) via SALib's built-in
resampling, giving the "confidence" the reviewer asked for.

Output per core output: mu* (overall importance, robust to sign),
mu (signed importance), sigma (interaction/non-linearity indicator), and
a bootstrap 95% CI on mu*. Nodes are ranked by mu*, and the rank of the
biology-of-interest nodes (CASP3, CASP7, TNFa, FasL, IL_8, TGF_b1, Cytc,
Bax) among all 58 is reported explicitly.

------------------------------------------------------------------
COMPUTE / RUNTIME
------------------------------------------------------------------
Number of model evaluations = R_TRAJ * (58 + 1). Each evaluation is one
300h ODE solve (~0.5-1s on a typical single core). Evaluations are
EMBARRASSINGLY PARALLEL (independent parameter sets), so this script uses
a multiprocessing.Pool sized to N_WORKERS (default: all logical CPU cores)
to parallelize across your machine's cores -- this typically gives a
near-linear speedup with core count. GPU acceleration is NOT applicable
here: the bottleneck is scipy's LSODA ODE solver, a small-dimension
(140-state), inherently sequential CPU routine with no practical GPU
implementation available for this network without a full rewrite.

Approximate wall-clock time for R_TRAJ=100 (5900 evaluations):
  1 core   : 5900 * ~0.6s / 1  ~= ~60 minutes
  4 cores  : ~15-18 minutes
  8 cores  : ~8-10 minutes
Lower R_TRAJ (e.g. 50) roughly halves both evaluation count and time, at
the cost of wider bootstrap CIs on mu*.

Resumable: raw (parameter_row, AUC_Cell_death, AUC_Hepatocyte_injury,
AUC_Inflammation) results are appended to outputs/gsa_raw_results.csv as
they complete; re-running with the same R_TRAJ/seed skips rows already
present in that file.

Usage:
    python run_morris_gsa.py
(edit R_TRAJ / N_WORKERS below, or just press F5 in Spyder)
"""
import os
import sys
import time
import numpy as np
import pandas as pd
from concurrent.futures import ProcessPoolExecutor, as_completed

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from gsa_common import (
    OUT_DIR, CORE_OUTPUTS, load_model, node_names, baseline_vector,
    build_bounds, eval_sample, _init_worker,
)

# ============================================================
# EDIT THESE (defaults chosen for a reasonable desktop run).
# ============================================================
R_TRAJ = 100          # number of Morris trajectories (evaluations = R_TRAJ*(58+1))
N_WORKERS = None      # None -> use os.cpu_count(); set an int to override
RNG_SEED = 20260908   # fixed seed -> reproducible sample matrix
NUM_LEVELS = 4        # Morris grid levels (SALib default-recommended value)

RAW_CSV = os.path.join(OUT_DIR, 'gsa_raw_results.csv')
XLSX_OUT = os.path.join(OUT_DIR, 'global_sensitivity_summary.xlsx')

# Nodes of particular interest (manuscript's local-sensitivity dominant
# drivers, Fig. 6 / Discussion) -- their rank among all 58 is reported.
NODES_OF_INTEREST = ['CASP7', 'CASP3', 'Cytc', 'Bax', 'TNFa', 'FasL', 'IL_8', 'TGF_b1']


def build_problem_and_samples():
    from SALib.sample import morris as morris_sample

    model = load_model()
    names = node_names(model)
    baseline = baseline_vector(model, names)
    bounds = build_bounds(baseline)

    problem = {'num_vars': len(names), 'names': names, 'bounds': bounds}
    np.random.seed(RNG_SEED)
    param_values = morris_sample.sample(problem, N=R_TRAJ, num_levels=NUM_LEVELS,
                                         seed=RNG_SEED)
    return problem, names, baseline, param_values


def run_evaluations(param_values):
    n_workers = N_WORKERS or os.cpu_count() or 1
    n_total = len(param_values)

    done_rows = {}
    if os.path.exists(RAW_CSV):
        prev = pd.read_csv(RAW_CSV)
        for _, row in prev.iterrows():
            done_rows[int(row['row_idx'])] = (
                row['AUC_P_Cell_death'], row['AUC_P_Hepatocyte_injury'], row['AUC_P_Inflammation'])
        print(f'[Resume] {len(done_rows)} / {n_total} evaluations already present in {RAW_CSV}')

    todo_idx = [i for i in range(n_total) if i not in done_rows]
    print(f'[Plan] {len(todo_idx)} new evaluations to run, using {n_workers} worker process(es).')

    if todo_idx:
        write_header = not os.path.exists(RAW_CSV)
        fh = open(RAW_CSV, 'a', newline='')
        if write_header:
            fh.write('row_idx,AUC_P_Cell_death,AUC_P_Hepatocyte_injury,AUC_P_Inflammation\n')
        fh.flush()

        t0 = time.time()
        n_done_this_run = 0
        with ProcessPoolExecutor(max_workers=n_workers, initializer=_init_worker) as ex:
            futures = {ex.submit(eval_sample, param_values[i]): i for i in todo_idx}
            for fut in as_completed(futures):
                i = futures[fut]
                auc_cd, auc_hi, auc_infl = fut.result()
                fh.write(f'{i},{auc_cd},{auc_hi},{auc_infl}\n')
                fh.flush()
                done_rows[i] = (auc_cd, auc_hi, auc_infl)
                n_done_this_run += 1
                if n_done_this_run % 50 == 0 or n_done_this_run == len(todo_idx):
                    elapsed = time.time() - t0
                    rate = elapsed / n_done_this_run
                    eta = rate * (len(todo_idx) - n_done_this_run)
                    print(f'  {n_done_this_run}/{len(todo_idx)} new evaluations done, '
                          f'elapsed={elapsed:.0f}s, eta={eta:.0f}s')
        fh.close()

    Y = np.zeros((n_total, 3))
    for i in range(n_total):
        Y[i] = done_rows[i]
    return Y


def analyze_and_save(problem, names, Y):
    from SALib.analyze import morris as morris_analyze
    from SALib.sample import morris as morris_sample

    # Rebuild the identical (seeded) parameter matrix X that produced Y,
    # since morris.analyze needs both X and Y together.
    np.random.seed(RNG_SEED)
    X = morris_sample.sample(problem, N=R_TRAJ, num_levels=NUM_LEVELS, seed=RNG_SEED)

    all_rows = []
    per_output_frames = {}
    for j, pw in enumerate(CORE_OUTPUTS):
        res = morris_analyze.analyze(
            problem, X, Y[:, j], num_resamples=1000, conf_level=0.95,
            print_to_console=False, seed=RNG_SEED)
        df = pd.DataFrame({
            'node': names,
            'mu_star': res['mu_star'],
            'mu_star_conf': res['mu_star_conf'],
            'mu': res['mu'],
            'sigma': res['sigma'],
        })
        df = df.sort_values('mu_star', ascending=False).reset_index(drop=True)
        df['rank'] = np.arange(1, len(df) + 1)
        per_output_frames[pw] = df

        print(f'\n=== {pw}: top 15 nodes by mu* (global sensitivity) ===')
        print(df.head(15).to_string(index=False))

        print(f'\n  Rank of nodes of interest for {pw}:')
        for node in NODES_OF_INTEREST:
            row = df[df['node'] == node]
            if len(row):
                r = int(row['rank'].values[0])
                ms = row['mu_star'].values[0]
                ci = row['mu_star_conf'].values[0]
                print(f'    {node:8s} rank {r:2d}/{len(df)}   mu*={ms:.4f} (CI ~+/-{ci:.4f})')

        all_rows.append(df.assign(output=pw))

    with pd.ExcelWriter(XLSX_OUT, engine='openpyxl') as writer:
        for pw, df in per_output_frames.items():
            df.to_excel(writer, sheet_name=pw[:31], index=False)
        pd.concat(all_rows, ignore_index=True).to_excel(writer, sheet_name='all_combined', index=False)
    print(f'\nExcel summary saved: {XLSX_OUT}')

    return per_output_frames


if __name__ == '__main__':
    print(f'R_TRAJ={R_TRAJ}  ->  {R_TRAJ * 59} model evaluations planned '
          f'(58 params + 1 per trajectory)')
    t_all = time.time()

    problem, names, baseline, param_values = build_problem_and_samples()
    Y = run_evaluations(param_values)
    per_output_frames = analyze_and_save(problem, names, Y)

    print(f'\nTotal wall-clock time: {(time.time() - t_all) / 60:.1f} minutes')
