# -*- coding: utf-8 -*-
"""Leave-one-out over the eight silymarin target nodes (R1-7 / R2-4 follow-up)."""
import os, sys, json
import numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from rtc_common import (MODEL_NASH_PATH, MODEL_NAFL_PATH, CORE_OUTPUTS, KI,
                        SILYMARIN_TARGETS, load_model, make_y0, build_ki_source,
                        make_mod_from_source, run_sim, pct_reduction, output_active_index)

m_nash = load_model(MODEL_NASH_PATH); m_nafl = load_model(MODEL_NAFL_PATH)
src_nash = open(MODEL_NASH_PATH, encoding='utf-8').read()
src_nafl = open(MODEL_NAFL_PATH, encoding='utf-8').read()
y0_nash, y0_nafl = make_y0(m_nash), make_y0(m_nafl)
base_nash, base_nafl = run_sim(m_nash, y0_nash), run_sim(m_nafl, y0_nafl)

def evaluate(targets, tag):
    dn = make_mod_from_source(build_ki_source(src_nash, m_nash, targets, KI), f'n_{tag}')
    df = make_mod_from_source(build_ki_source(src_nafl, m_nafl, targets, KI), f'f_{tag}')
    tn, tf = run_sim(dn, y0_nash), run_sim(df, y0_nafl)
    out = {}
    for pw in CORE_OUTPUTS:
        i = output_active_index(m_nash, pw)
        nash = pct_reduction(tn[:, i], base_nash[:, i])
        nafl = pct_reduction(tf[:, i], base_nafl[:, i])
        out[pw] = {'NASH_red_pct': nash, 'NAFL_red_pct': nafl,
                   'ratio': (nafl / nash) if nash > 0 else np.nan,
                   'diff_pp': nafl - nash}
    return out

rows = []
full = evaluate(SILYMARIN_TARGETS, 'full')
for pw in CORE_OUTPUTS:
    rows.append({'scenario': 'full 8-target set', 'omitted_target': '-', 'n_targets': 8,
                 'output': pw, **{k: round(v, 4) for k, v in full[pw].items()}})
for t in SILYMARIN_TARGETS:
    sub = [x for x in SILYMARIN_TARGETS if x != t]
    r = evaluate(sub, t)
    for pw in CORE_OUTPUTS:
        rows.append({'scenario': f'leave-one-out (-{t})', 'omitted_target': t, 'n_targets': 7,
                     'output': pw, **{k: round(v, 4) for k, v in r[pw].items()}})
df = pd.DataFrame(rows)
os.makedirs('outputs', exist_ok=True)
df.to_excel('outputs/leave_one_out_summary.xlsx', index=False)
piv = df.pivot_table(index=['scenario', 'omitted_target'], columns='output',
                     values='ratio', sort=False)
print(piv.to_string())
print()
bad = df[(df.n_targets == 7) & (df.ratio <= 1)]
print('leave-one-out combinations with advantage ratio <= 1:', len(bad))
if len(bad): print(bad.to_string())
print()
print('ratio range per output across the 8 leave-one-out runs:')
for pw in CORE_OUTPUTS:
    v = df[(df.n_targets == 7) & (df.output == pw)]['ratio']
    print(f'  {pw:22s} [{v.min():.3f}, {v.max():.3f}]   full-set value {full[pw]["ratio"]:.3f}')
