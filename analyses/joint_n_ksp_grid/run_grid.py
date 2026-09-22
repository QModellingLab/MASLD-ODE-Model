#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Joint n x ksp parameter grid (R3, minor point 3).

The manuscript's Supplementary Table S4 varies the Hill coefficient (n) and
half-saturation constant (ksp) ONE AT A TIME from the primary n=ksp=2 model
(n=1,2,4 at ksp=2; ksp=1,4 at n=2 -- 5 combinations total including the
primary). Reviewer 3 asks to vary them JOINTLY (a full grid) and to report
the EFFECT MAGNITUDE, not just its sign/direction.

This script runs the full 3x3 grid (n in {1,2,4}, ksp in {1,2,4} = 9
combinations, including all 5 already in Table S4 plus 4 new joint
combinations: (1,1),(1,4),(4,1),(4,4)) and reports, for each of the 3 core
outputs: NASH %reduction, NAFL %reduction, the NAFL/NASH advantage ratio
(same metric as Table 4 / Table S4), and the raw difference in percentage
points (diff_pp) as an explicit magnitude measure.

Fast, deterministic, single run -- no batching or resampling needed
(9 combinations x 4 ODE solves = 36 solves total, ~20-30s).

Usage: python run_grid.py  (no parameters; just press F5)
"""
import os
import numpy as np
import pandas as pd

from nksp_common import (
    MODEL_NASH_PATH, MODEL_NAFL_PATH, OUT_DIR, CORE_OUTPUTS, KI,
    SILYMARIN_TARGETS, load_model, make_y0, parameterize_nksp,
    build_ki_source, make_mod_from_source, run_sim, pct_reduction,
    output_active_index,
)

N_VALUES = [1, 2, 4]
KSP_VALUES = [1, 2, 4]


def run_combo(base_nash, base_nafl, n, ksp):
    src_nash = parameterize_nksp(base_nash, n, ksp)
    src_nafl = parameterize_nksp(base_nafl, n, ksp)
    m_nash = make_mod_from_source(src_nash, f'nash_{n}_{ksp}')
    m_nafl = make_mod_from_source(src_nafl, f'nafl_{n}_{ksp}')
    y0_nash = make_y0(m_nash)
    y0_nafl = make_y0(m_nafl)
    drug_nash = make_mod_from_source(
        build_ki_source(src_nash, m_nash, SILYMARIN_TARGETS, KI), f'dnash_{n}_{ksp}')
    drug_nafl = make_mod_from_source(
        build_ki_source(src_nafl, m_nafl, SILYMARIN_TARGETS, KI), f'dnafl_{n}_{ksp}')

    tn0, tn1 = run_sim(m_nash, y0_nash), run_sim(drug_nash, y0_nash)
    tf0, tf1 = run_sim(m_nafl, y0_nafl), run_sim(drug_nafl, y0_nafl)

    rows = []
    for pw in CORE_OUTPUTS:
        i = output_active_index(m_nash, pw)
        nash_red = pct_reduction(tn1[:, i], tn0[:, i])
        nafl_red = pct_reduction(tf1[:, i], tf0[:, i])
        ratio = nafl_red / nash_red if nash_red > 0 else np.nan
        rows.append({
            'n': n, 'ksp': ksp, 'primary': (n == 2 and ksp == 2),
            'in_manuscript_table_S4': (n, ksp) in [(1, 2), (2, 2), (4, 2), (2, 1), (2, 4)],
            'output': pw, 'NASH_red_pct': nash_red, 'NAFL_red_pct': nafl_red,
            'advantage_ratio': ratio, 'diff_pp': nafl_red - nash_red,
        })
    return rows


def main():
    with open(MODEL_NASH_PATH, encoding='utf-8') as f:
        base_nash = f.read()
    with open(MODEL_NAFL_PATH, encoding='utf-8') as f:
        base_nafl = f.read()

    all_rows = []
    print(f'{"n":>3} {"ksp":>4}  {"output":22s} {"NASH%":>8s} {"NAFL%":>8s} '
          f'{"ratio":>7s} {"diff_pp":>8s}')
    for n in N_VALUES:
        for ksp in KSP_VALUES:
            rows = run_combo(base_nash, base_nafl, n, ksp)
            for r in rows:
                print(f'{r["n"]:>3} {r["ksp"]:>4}  {r["output"]:22s} '
                      f'{r["NASH_red_pct"]:8.3f} {r["NAFL_red_pct"]:8.3f} '
                      f'{r["advantage_ratio"]:7.3f} {r["diff_pp"]:8.3f}')
            all_rows.extend(rows)

    df = pd.DataFrame(all_rows)
    xlsx_path = os.path.join(OUT_DIR, 'joint_n_ksp_grid_summary.xlsx')
    with pd.ExcelWriter(xlsx_path, engine='openpyxl') as writer:
        df.to_excel(writer, sheet_name='all_9_combos', index=False)
        pivot_ratio = df.pivot_table(index=['n', 'ksp'], columns='output', values='advantage_ratio')
        pivot_ratio.to_excel(writer, sheet_name='ratio_pivot')
        pivot_diff = df.pivot_table(index=['n', 'ksp'], columns='output', values='diff_pp')
        pivot_diff.to_excel(writer, sheet_name='diff_pp_pivot')
    print(f'\nExcel summary saved: {xlsx_path}')

    print('\n-- Summary: does advantage_ratio > 1 (NAFL > NASH) hold for ALL 9 combos? --')
    for pw in CORE_OUTPUTS:
        sub = df[df['output'] == pw]
        all_hold = (sub['advantage_ratio'] > 1).all()
        rng = (sub['advantage_ratio'].min(), sub['advantage_ratio'].max())
        print(f'  {pw:22s} all_hold={all_hold}   ratio range=[{rng[0]:.3f}, {rng[1]:.3f}]')


if __name__ == '__main__':
    main()
