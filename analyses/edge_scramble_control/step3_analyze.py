#!/usr/bin/env python3
"""
step3_analyze.py
=================
Reads outputs/edge_scramble_true_result.json and outputs/{scramble,er}_results.jsonl,
computes marginal (per-output) and joint (all-three-outputs-simultaneously)
empirical significance tests against both null models, with Wilson
confidence intervals, and writes a summary Excel file plus a printed
report.

Two test statistics are reported for the marginal test (see README for
why): the percentage-point DIFFERENCE (NAFL_red - NASH_red), which is
robust to near-zero-denominator blow-ups, and is used as primary; the
RATIO (NAFL_red / NASH_red), matching the manuscript's headline metric,
reported for reference only.

The joint test (all three core outputs simultaneously exceeding the true
topology's per-output difference) is pre-specified by reference to the
manuscript's own Supplementary Table S4 "all_outputs_hold" column design
-- it is not a metric chosen post hoc to rescue a null marginal result,
though we report the marginal result plainly alongside it either way.

Usage:
    python step3_analyze.py

NOTE ON FILENAME (2026-09-11): reads outputs/edge_scramble_true_result.json
(previously outputs/true_result.json, renamed to avoid a filename
collision with random_target_control's own true_result.json when both
are uploaded into the same flat Project knowledge folder).
"""
import json
import math
import statistics

CORE_OUTPUTS = ['P_Hepatocyte_injury', 'P_Cell_death', 'P_Inflammation']


def wilson_ci(k, n, z=1.96):
    if n == 0:
        return (float('nan'), float('nan'))
    phat = k / n
    denom = 1 + z**2 / n
    center = (phat + z**2 / (2 * n)) / denom
    half = z * math.sqrt(phat * (1 - phat) / n + z**2 / (4 * n**2)) / denom
    return (max(0, center - half), min(1, center + half))


def load_clean(path):
    data = [json.loads(l) for l in open(path)]
    clean = [d for d in data if not any(math.isnan(r['ratio']) for r in d['result'].values())]
    return clean, len(data) - len(clean)


def analyze_one_null(label, path, true_result):
    print("\n" + "=" * 74)
    print(f"NULL MODEL: {label}")
    print("=" * 74)
    clean, n_excluded = load_clean(path)
    n = len(clean)
    print(f"Replicates: {n} clean ({n_excluded} excluded for numerical NaN/Inf)")

    rows = []
    print("\n-- Marginal (per-output) tests, difference-based --")
    for node in CORE_OUTPUTS:
        diffs = [d['result'][node]['nafl_red'] - d['result'][node]['nash_red'] for d in clean]
        true_diff = true_result[node]['nafl_red'] - true_result[node]['nash_red']
        k_ge = sum(1 for x in diffs if x >= true_diff)
        lo, hi = wilson_ci(k_ge, n)
        k_dir = sum(1 for x in diffs if x > 0)
        dlo, dhi = wilson_ci(k_dir, n)
        print(f"  {node:<22} true diff={true_diff:6.2f}pp  "
              f"null median={statistics.median(diffs):6.2f}pp  "
              f"direction={100*k_dir/n:5.1f}% (CI {100*dlo:.1f}-{100*dhi:.1f}%)  "
              f"P(>=true)={k_ge/n:.4f} (CI {lo:.4f}-{hi:.4f})")
        rows.append({
            'null_model': label, 'test': 'marginal', 'output': node,
            'true_diff_pp': round(true_diff, 3),
            'null_median_diff_pp': round(statistics.median(diffs), 3),
            'direction_pct': round(100 * k_dir / n, 2),
            'direction_ci_lo': round(100 * dlo, 2), 'direction_ci_hi': round(100 * dhi, 2),
            'p_value': round(k_ge / n, 5), 'p_ci_lo': round(lo, 5), 'p_ci_hi': round(hi, 5),
            'n_replicates': n,
        })

    print("\n-- Joint test (all three outputs simultaneously) --")
    joint_dir_k = sum(1 for d in clean if all(
        d['result'][node]['nafl_red'] - d['result'][node]['nash_red'] > 0 for node in CORE_OUTPUTS))
    jlo, jhi = wilson_ci(joint_dir_k, n)
    print(f"  All-three direction (diff>0 for all 3): {joint_dir_k}/{n} = {100*joint_dir_k/n:.1f}% "
          f"(CI {100*jlo:.1f}-{100*jhi:.1f}%)")

    true_diffs = {node: true_result[node]['nafl_red'] - true_result[node]['nash_red'] for node in CORE_OUTPUTS}
    joint_ge_k = sum(1 for d in clean if all(
        d['result'][node]['nafl_red'] - d['result'][node]['nash_red'] >= true_diffs[node]
        for node in CORE_OUTPUTS))
    jglo, jghi = wilson_ci(joint_ge_k, n)
    print(f"  All-three diff >= true (joint): {joint_ge_k}/{n} = {joint_ge_k/n:.4f} "
          f"(CI {jglo:.4f}-{jghi:.4f})")

    true_min_diff = min(true_diffs.values())
    null_min_diffs = [min(d['result'][node]['nafl_red'] - d['result'][node]['nash_red']
                           for node in CORE_OUTPUTS) for d in clean]
    k_min = sum(1 for m in null_min_diffs if m >= true_min_diff)
    mlo, mhi = wilson_ci(k_min, n)
    print(f"  Min-diff-across-3 summary stat: true={true_min_diff:.2f}pp, "
          f"null median={statistics.median(null_min_diffs):.2f}pp, "
          f"P(>=true)={k_min/n:.4f} (CI {mlo:.4f}-{mhi:.4f})")

    rows.append({
        'null_model': label, 'test': 'joint_direction', 'output': 'all_three',
        'direction_pct': round(100 * joint_dir_k / n, 2),
        'direction_ci_lo': round(100 * jlo, 2), 'direction_ci_hi': round(100 * jhi, 2),
        'n_replicates': n,
    })
    rows.append({
        'null_model': label, 'test': 'joint_diff_ge_true', 'output': 'all_three',
        'p_value': round(joint_ge_k / n, 5), 'p_ci_lo': round(jglo, 5), 'p_ci_hi': round(jghi, 5),
        'n_replicates': n,
    })
    rows.append({
        'null_model': label, 'test': 'joint_min_diff', 'output': 'all_three',
        'true_diff_pp': round(true_min_diff, 3),
        'null_median_diff_pp': round(statistics.median(null_min_diffs), 3),
        'p_value': round(k_min / n, 5), 'p_ci_lo': round(mlo, 5), 'p_ci_hi': round(mhi, 5),
        'n_replicates': n,
    })
    return rows


def main():
    true_result = json.load(open('outputs/edge_scramble_true_result.json'))
    print("TRUE (real hsa04932) topology reference:")
    for node, r in true_result.items():
        print(f"  {node:<22} NAFL={r['nafl_red']:.2f}%  NASH={r['nash_red']:.2f}%  "
              f"ratio={r['ratio']:.2f}  diff={r['nafl_red']-r['nash_red']:.2f}pp")

    all_rows = []
    import os
    if os.path.exists('outputs/scramble_results.jsonl'):
        all_rows += analyze_one_null('degree-preserving scramble (Maslov-Sneppen)',
                                      'outputs/scramble_results.jsonl', true_result)
    else:
        print("\n(outputs/scramble_results.jsonl not found -- run step2_run_batch.py scramble first)")

    if os.path.exists('outputs/er_results.jsonl'):
        all_rows += analyze_one_null('Erdos-Renyi (degree-unconstrained)',
                                      'outputs/er_results.jsonl', true_result)
    else:
        print("\n(outputs/er_results.jsonl not found -- run step2_run_batch.py er first)")

    try:
        import pandas as pd
        df = pd.DataFrame(all_rows)
        with pd.ExcelWriter('outputs/edge_scramble_summary.xlsx', engine='openpyxl') as w:
            df.to_excel(w, sheet_name='Summary', index=False)
            pd.DataFrame([{'output': k, **v} for k, v in true_result.items()]).to_excel(
                w, sheet_name='True_topology_reference', index=False)
        print("\nExcel summary saved: outputs/edge_scramble_summary.xlsx")
    except ImportError:
        print("\n(pandas/openpyxl not installed -- skipping Excel export; "
              "printed report above has all the numbers)")


if __name__ == '__main__':
    main()
