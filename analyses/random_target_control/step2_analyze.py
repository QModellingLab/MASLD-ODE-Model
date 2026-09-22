#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Analyze the random_target / matched_control null distributions against the
true 8-silymarin-target result. Run after step1_run_batch.py has produced
outputs/true_result.json and at least one of outputs/random_target_results.jsonl
/ outputs/matched_control_results.jsonl.

For each null model and each of the 3 core outputs, reports:
  - true diff (NAFL_red - NASH_red, percentage points) vs null median
  - "direction" = % of null replicates with diff > 0 (NAFL>NASH direction)
  - P(>=true) = one-sided empirical P-value (proportion of null replicates
    whose diff meets or exceeds the true diff) -- this is the R1-7-style
    "FDR/percentile significance" of the real silymarin target set relative
    to a random-8-node background, with a bootstrap 95% CI on that
    proportion.
  - the analogous statistic for the NAFL/NASH ratio itself.
Also reports a joint (all-three-outputs-simultaneously) test and a
min-diff-across-3 summary statistic, exactly parallel to the edge_scramble_
control analysis, so the two control analyses can be reported side by side
in the Response to Reviewers.

Usage:
    python step2_analyze.py

NOTE ON FILENAME (2026-09-11): reads outputs/random_target_control_true_result.json
(previously true_result.json, renamed to avoid a filename collision with
edge_scramble_control's own true_result.json when both are uploaded into
the same flat Project knowledge folder).
"""
import json, os
import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
OUT_DIR = os.path.join(HERE, 'outputs')
CORE_OUTPUTS = ['P_Cell_death', 'P_Hepatocyte_injury', 'P_Inflammation']
N_BOOT = 2000
RNG = np.random.default_rng(12345)


def load_null(mode):
    path = os.path.join(OUT_DIR, f'{mode}_results.jsonl')
    if not os.path.exists(path):
        return None
    recs = []
    with open(path) as fh:
        for line in fh:
            rec = json.loads(line)
            if 'error' not in rec and not rec.get('nan_flag', False):
                recs.append(rec)
    return recs


def boot_ci_prop(vals, true_val, n_boot=N_BOOT):
    vals = np.asarray(vals)
    n = len(vals)
    obs = np.mean(vals >= true_val)
    boots = np.empty(n_boot)
    for b in range(n_boot):
        sample = vals[RNG.integers(0, n, n)]
        boots[b] = np.mean(sample >= true_val)
    lo, hi = np.percentile(boots, [2.5, 97.5])
    return obs, lo, hi


def boot_ci_mean(vals, n_boot=N_BOOT):
    vals = np.asarray(vals)
    n = len(vals)
    obs = np.mean(vals)
    boots = np.empty(n_boot)
    for b in range(n_boot):
        sample = vals[RNG.integers(0, n, n)]
        boots[b] = np.mean(sample)
    lo, hi = np.percentile(boots, [2.5, 97.5])
    return obs, lo, hi


def analyze_mode(mode, true_result, records):
    n_total = len(records)
    n_excluded = sum(1 for r in records if any(
        r.get(f'{pw}_diff_pp') is None for pw in CORE_OUTPUTS))
    clean = [r for r in records if all(
        r.get(f'{pw}_diff_pp') is not None for pw in CORE_OUTPUTS)]
    print(f'\n{"="*78}\nNULL MODEL: {mode}  (N={len(clean)} clean replicates, '
          f'{n_excluded} excluded for NaN/Inf)\n{"="*78}')

    summary_rows = []
    diff_mat = np.array([[r[f'{pw}_diff_pp'] for pw in CORE_OUTPUTS] for r in clean])

    print('\n-- Marginal (per-output) tests, difference-based (NAFL_red - NASH_red, pp) --')
    for j, pw in enumerate(CORE_OUTPUTS):
        true_diff = true_result[pw]['diff_pp']
        null_vals = diff_mat[:, j]
        direction_pct, dlo, dhi = boot_ci_prop(null_vals, 1e-9)  # % with diff>0 (approx via >=~0)
        direction_pct = np.mean(null_vals > 0) * 100
        # bootstrap CI for direction%
        boots_dir = np.array([np.mean(null_vals[RNG.integers(0, len(null_vals), len(null_vals))] > 0) * 100
                               for _ in range(N_BOOT)])
        dlo, dhi = np.percentile(boots_dir, [2.5, 97.5])
        p_obs, p_lo, p_hi = boot_ci_prop(null_vals, true_diff)
        null_median = np.median(null_vals)
        print(f'  {pw:22s} true diff={true_diff:7.2f}pp  null median={null_median:6.2f}pp  '
              f'direction={direction_pct:5.1f}% (CI {dlo:.1f}-{dhi:.1f}%)  '
              f'P(>=true)={p_obs:.4f} (CI {p_lo:.4f}-{p_hi:.4f})')
        summary_rows.append({
            'null_model': mode, 'output': pw, 'metric': 'diff_pp',
            'true_value': true_diff, 'null_median': null_median,
            'direction_pct': direction_pct, 'direction_CI_lo': dlo, 'direction_CI_hi': dhi,
            'P_ge_true': p_obs, 'P_ge_true_CI_lo': p_lo, 'P_ge_true_CI_hi': p_hi,
        })

    print('\n-- Marginal (per-output) tests, ratio-based (NAFL_red / NASH_red) --')
    ratio_mat = np.array([[r[f'{pw}_ratio'] if r.get(f'{pw}_ratio') is not None else np.nan
                            for pw in CORE_OUTPUTS] for r in clean])
    for j, pw in enumerate(CORE_OUTPUTS):
        true_ratio = true_result[pw]['ratio']
        vals = ratio_mat[:, j]
        vals = vals[np.isfinite(vals)]
        p_obs, p_lo, p_hi = boot_ci_prop(vals, true_ratio)
        null_median = np.median(vals)
        print(f'  {pw:22s} true ratio={true_ratio:7.3f}  null median ratio={null_median:6.3f}  '
              f'P(>=true)={p_obs:.4f} (CI {p_lo:.4f}-{p_hi:.4f})  n_valid={len(vals)}')
        summary_rows.append({
            'null_model': mode, 'output': pw, 'metric': 'ratio',
            'true_value': true_ratio, 'null_median': null_median,
            'direction_pct': None, 'direction_CI_lo': None, 'direction_CI_hi': None,
            'P_ge_true': p_obs, 'P_ge_true_CI_lo': p_lo, 'P_ge_true_CI_hi': p_hi,
        })

    print('\n-- Joint test (all three outputs simultaneously) --')
    true_diffs = np.array([true_result[pw]['diff_pp'] for pw in CORE_OUTPUTS])
    all_three_pos = np.mean(np.all(diff_mat > 0, axis=1)) * 100
    all_three_ge_true = np.mean(np.all(diff_mat >= true_diffs[None, :], axis=1))
    obs_ge, lo_ge, hi_ge = boot_ci_prop(np.all(diff_mat > 0, axis=1).astype(float), 1.0) if False else (None, None, None)
    boots_all3 = np.array([np.mean(np.all(diff_mat[RNG.integers(0, len(diff_mat), len(diff_mat))] > 0, axis=1)) * 100
                            for _ in range(N_BOOT)])
    a3_lo, a3_hi = np.percentile(boots_all3, [2.5, 97.5])
    boots_all3ge = np.array([np.mean(np.all(diff_mat[RNG.integers(0, len(diff_mat), len(diff_mat))] >= true_diffs[None, :], axis=1))
                              for _ in range(N_BOOT)])
    a3ge_lo, a3ge_hi = np.percentile(boots_all3ge, [2.5, 97.5])
    print(f'  All-three direction (diff>0 for all 3): {all_three_pos:.1f}% (CI {a3_lo:.1f}-{a3_hi:.1f}%)')
    print(f'  All-three diff >= true (joint): {all_three_ge_true:.4f} (CI {a3ge_lo:.4f}-{a3ge_hi:.4f})')

    min_diff_true = true_diffs.min()
    min_diff_null = diff_mat.min(axis=1)
    p_min, p_min_lo, p_min_hi = boot_ci_prop(min_diff_null, min_diff_true)
    print(f'  Min-diff-across-3 summary stat: true={min_diff_true:.2f}pp, '
          f'null median={np.median(min_diff_null):.2f}pp, '
          f'P(>=true)={p_min:.4f} (CI {p_min_lo:.4f}-{p_min_hi:.4f})')

    summary_rows.append({
        'null_model': mode, 'output': 'ALL_THREE_joint', 'metric': 'diff_pp_direction_pct',
        'true_value': None, 'null_median': None,
        'direction_pct': all_three_pos, 'direction_CI_lo': a3_lo, 'direction_CI_hi': a3_hi,
        'P_ge_true': all_three_ge_true, 'P_ge_true_CI_lo': a3ge_lo, 'P_ge_true_CI_hi': a3ge_hi,
    })
    summary_rows.append({
        'null_model': mode, 'output': 'MIN_diff_across_3', 'metric': 'diff_pp',
        'true_value': min_diff_true, 'null_median': np.median(min_diff_null),
        'direction_pct': None, 'direction_CI_lo': None, 'direction_CI_hi': None,
        'P_ge_true': p_min, 'P_ge_true_CI_lo': p_min_lo, 'P_ge_true_CI_hi': p_min_hi,
    })

    return summary_rows, len(clean), n_excluded


def main():
    true_path = os.path.join(OUT_DIR, 'random_target_control_true_result.json')
    if not os.path.exists(true_path):
        print('[ERROR] outputs/random_target_control_true_result.json not found -- run '
              '"python step1_run_batch.py true" first.')
        return
    with open(true_path) as fh:
        true_result = json.load(fh)

    print('TRUE (real hsa04932 topology, real 8 silymarin targets) reference:')
    for pw in CORE_OUTPUTS:
        r = true_result[pw]
        print(f'  {pw:22s} NAFL={r["NAFL_red"]:.2f}%  NASH={r["NASH_red"]:.2f}%  '
              f'ratio={r["ratio"]:.2f}  diff={r["diff_pp"]:.2f}pp')

    all_rows = []
    for mode in ['random_target', 'matched_control']:
        records = load_null(mode)
        if records is None:
            print(f'\n[skip] outputs/{mode}_results.jsonl not found yet.')
            continue
        rows, n_clean, n_exc = analyze_mode(mode, true_result, records)
        all_rows.extend(rows)

    if not all_rows:
        print('\nNo null-model results found yet -- nothing to summarize.')
        return

    df = pd.DataFrame(all_rows)
    xlsx_path = os.path.join(OUT_DIR, 'random_target_control_summary.xlsx')
    df.to_excel(xlsx_path, index=False)
    print(f'\nExcel summary saved: {xlsx_path}')


if __name__ == '__main__':
    main()
