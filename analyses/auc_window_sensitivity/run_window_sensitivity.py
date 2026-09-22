#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
AUC time-window sensitivity (R1-5 / R3 minor point 1).

Reviewer 1, Major Comment 5:
  "the kinetic parameters were not calibrated using time-course data, so it
   is unclear why the time axis is presented in hours. Authors should ...
   test whether the AUC results are robust to different simulation windows
   and outcome measurements."
Reviewer 3, Minor Point 1:
  "test output sensitivity to ... the 0 to 300 h AUC window."

This script simulates each condition (NASH/NAFL, drug/no-drug) ONCE over a
long horizon (0-1000h, same 0.1h resolution as the manuscript's 300h run),
then computes the AUC-based %reduction and NAFL/NASH advantage ratio at
FOUR different window widths (100h, 300h [manuscript], 600h, 1000h) as
sub-slices of that single trajectory -- no re-simulation needed per window,
so this is fast (4 ODE solves total, ~5-10s).

Usage: python run_window_sensitivity.py (no parameters; just press F5)
"""
import os
import numpy as np
import pandas as pd

from aucwin_common import (
    MODEL_NASH_PATH, MODEL_NAFL_PATH, OUT_DIR, CORE_OUTPUTS, KI,
    SILYMARIN_TARGETS, load_model, make_y0, build_ki_source,
    make_mod_from_source, run_sim_long, pct_reduction_window,
    output_active_index,
)

WINDOWS_H = [100.0, 300.0, 600.0, 1000.0]


def main():
    m_nash = load_model(MODEL_NASH_PATH)
    m_nafl = load_model(MODEL_NAFL_PATH)
    y0_nash = make_y0(m_nash)
    y0_nafl = make_y0(m_nafl)

    with open(MODEL_NASH_PATH, encoding='utf-8') as f:
        src_nash = f.read()
    with open(MODEL_NAFL_PATH, encoding='utf-8') as f:
        src_nafl = f.read()
    drug_nash = make_mod_from_source(build_ki_source(src_nash, m_nash, SILYMARIN_TARGETS, KI), 'drug_nash')
    drug_nafl = make_mod_from_source(build_ki_source(src_nafl, m_nafl, SILYMARIN_TARGETS, KI), 'drug_nafl')

    print('Running 4 long-horizon (0-1000h) ODE solves (NASH/NAFL x drug/no-drug)...')
    tn0 = run_sim_long(m_nash, y0_nash)
    tn1 = run_sim_long(drug_nash, y0_nash)
    tf0 = run_sim_long(m_nafl, y0_nafl)
    tf1 = run_sim_long(drug_nafl, y0_nafl)
    print('Done. Slicing AUC at each window...')

    rows = []
    print(f'\n{"window_h":>9} {"output":22s} {"NASH%":>8s} {"NAFL%":>8s} '
          f'{"ratio":>7s} {"diff_pp":>8s}')
    for w in WINDOWS_H:
        for pw in CORE_OUTPUTS:
            i = output_active_index(m_nash, pw)
            nash_red = pct_reduction_window(tn1[:, i], tn0[:, i], w)
            nafl_red = pct_reduction_window(tf1[:, i], tf0[:, i], w)
            ratio = nafl_red / nash_red if nash_red > 0 else np.nan
            rows.append({
                'window_h': w, 'is_manuscript_default': (w == 300.0),
                'output': pw, 'NASH_red_pct': nash_red, 'NAFL_red_pct': nafl_red,
                'advantage_ratio': ratio, 'diff_pp': nafl_red - nash_red,
            })
            print(f'{w:>9.0f} {pw:22s} {nash_red:8.3f} {nafl_red:8.3f} '
                  f'{ratio:7.3f} {nafl_red-nash_red:8.3f}')

    df = pd.DataFrame(rows)
    xlsx_path = os.path.join(OUT_DIR, 'auc_window_sensitivity_summary.xlsx')
    with pd.ExcelWriter(xlsx_path, engine='openpyxl') as writer:
        df.to_excel(writer, sheet_name='all_windows', index=False)
        pivot = df.pivot_table(index='window_h', columns='output', values='advantage_ratio')
        pivot.to_excel(writer, sheet_name='ratio_pivot')
    print(f'\nExcel summary saved: {xlsx_path}')

    print('\n-- Summary: does advantage_ratio > 1 (NAFL > NASH) hold across ALL windows? --')
    for pw in CORE_OUTPUTS:
        sub = df[df['output'] == pw]
        all_hold = (sub['advantage_ratio'] > 1).all()
        rng = (sub['advantage_ratio'].min(), sub['advantage_ratio'].max())
        print(f'  {pw:22s} all_hold={all_hold}   ratio range=[{rng[0]:.3f}, {rng[1]:.3f}]')


if __name__ == '__main__':
    main()
