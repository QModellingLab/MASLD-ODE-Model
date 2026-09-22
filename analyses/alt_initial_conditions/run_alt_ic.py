#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Alternative initial-condition assumptions (R1-4).

Reviewer 1, Major Comment 4 asks:
  "Why should mRNA expression represent inactive protein concentration, and
   why should all signaling nodes have no initial activity in established
   NAFL or NASH tissue? Authors also set the inactive value of each pathway
   output to 100 without providing a biological explanation. These
   assumptions may strongly affect the simulated trajectories and should be
   justified and tested using alternative initial settings."

This script tests two of the manuscript's initial-condition conventions by
sweeping them over a grid and re-running the SAME silymarin AUC-reduction
metric used throughout the paper (Table 4-style: NASH/NAFL %reduction,
NAFL/NASH advantage ratio):

  1. active_fraction: instead of ALL active forms starting at exactly 0
     (manuscript default), initialize each molecular node's active form as
     a FRACTION of its own inactive-form initial ratio (0, 0.1, 0.3 tested
     -- i.e. 0%, 10%, 30% baseline activity already present in established
     disease tissue).
  2. p_output_inactive_level: the fixed starting value for the 12 pathway-
     output nodes' inactive forms (manuscript default 100; 50 and 200 also
     tested).

Full 3x3 = 9 combinations, each requiring 4 ODE solves (NASH/NAFL x
drug/no-drug) -- fast, deterministic, single run (~30-40s total).

Usage: python run_alt_ic.py (no parameters; just press F5)
"""
import os
import numpy as np
import pandas as pd

from altic_common import (
    MODEL_NASH_PATH, MODEL_NAFL_PATH, OUT_DIR, CORE_OUTPUTS, KI,
    SILYMARIN_TARGETS, load_model, make_y0, build_ki_source,
    make_mod_from_source, run_sim, pct_reduction, output_active_index,
)

ACTIVE_FRACTIONS = [0.0, 0.1, 0.3]
P_OUTPUT_LEVELS = [50.0, 100.0, 200.0]


def run_combo(m_nash, m_nafl, drug_nash, drug_nafl, active_frac, p_out_level):
    y0_nash = make_y0(m_nash, active_frac, p_out_level)
    y0_nafl = make_y0(m_nafl, active_frac, p_out_level)

    tn0, tn1 = run_sim(m_nash, y0_nash), run_sim(drug_nash, y0_nash)
    tf0, tf1 = run_sim(m_nafl, y0_nafl), run_sim(drug_nafl, y0_nafl)

    rows = []
    for pw in CORE_OUTPUTS:
        i = output_active_index(m_nash, pw)
        nash_red = pct_reduction(tn1[:, i], tn0[:, i])
        nafl_red = pct_reduction(tf1[:, i], tf0[:, i])
        ratio = nafl_red / nash_red if nash_red > 0 else np.nan
        rows.append({
            'active_fraction': active_frac, 'p_output_inactive_level': p_out_level,
            'is_manuscript_default': (active_frac == 0.0 and p_out_level == 100.0),
            'output': pw, 'NASH_red_pct': nash_red, 'NAFL_red_pct': nafl_red,
            'advantage_ratio': ratio, 'diff_pp': nafl_red - nash_red,
        })
    return rows


def main():
    m_nash = load_model(MODEL_NASH_PATH)
    m_nafl = load_model(MODEL_NAFL_PATH)
    with open(MODEL_NASH_PATH, encoding='utf-8') as f:
        src_nash = f.read()
    with open(MODEL_NAFL_PATH, encoding='utf-8') as f:
        src_nafl = f.read()
    drug_nash = make_mod_from_source(build_ki_source(src_nash, m_nash, SILYMARIN_TARGETS, KI), 'drug_nash')
    drug_nafl = make_mod_from_source(build_ki_source(src_nafl, m_nafl, SILYMARIN_TARGETS, KI), 'drug_nafl')

    all_rows = []
    print(f'{"active_frac":>11} {"P_out_level":>11}  {"output":22s} {"NASH%":>8s} '
          f'{"NAFL%":>8s} {"ratio":>7s} {"diff_pp":>8s}')
    for af in ACTIVE_FRACTIONS:
        for pol in P_OUTPUT_LEVELS:
            rows = run_combo(m_nash, m_nafl, drug_nash, drug_nafl, af, pol)
            for r in rows:
                print(f'{r["active_fraction"]:>11.2f} {r["p_output_inactive_level"]:>11.0f}  '
                      f'{r["output"]:22s} {r["NASH_red_pct"]:8.3f} {r["NAFL_red_pct"]:8.3f} '
                      f'{r["advantage_ratio"]:7.3f} {r["diff_pp"]:8.3f}')
            all_rows.extend(rows)

    df = pd.DataFrame(all_rows)
    xlsx_path = os.path.join(OUT_DIR, 'alt_initial_conditions_summary.xlsx')
    with pd.ExcelWriter(xlsx_path, engine='openpyxl') as writer:
        df.to_excel(writer, sheet_name='all_9_combos', index=False)
        pivot = df.pivot_table(index=['active_fraction', 'p_output_inactive_level'],
                                columns='output', values='advantage_ratio')
        pivot.to_excel(writer, sheet_name='ratio_pivot')
    print(f'\nExcel summary saved: {xlsx_path}')

    print('\n-- Summary: does advantage_ratio > 1 (NAFL > NASH) hold for ALL 9 combos? --')
    for pw in CORE_OUTPUTS:
        sub = df[df['output'] == pw]
        all_hold = (sub['advantage_ratio'] > 1).all()
        rng = (sub['advantage_ratio'].min(), sub['advantage_ratio'].max())
        print(f'  {pw:22s} all_hold={all_hold}   ratio range=[{rng[0]:.3f}, {rng[1]:.3f}]')


if __name__ == '__main__':
    main()
