#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
plot_supp_s6.py
================================================================
Independent plotting script for Supplementary Figure S6 (parameter
robustness of the stage-dependent response).

Does NOT re-derive any numbers of its own -- it reuses the exact same
verified computation utilities from joint_n_ksp_grid/nksp_common.py
(the same module already confirmed, via run_grid.py, to reproduce the
manuscript's published numbers exactly, including the (n=1,ksp=1)
P_Hepatocyte_injury reversal to 0.952).

Two outputs:
  1. A 5-point subset plot (n1/ksp2, n2/ksp2, n4/ksp2, n2/ksp1, n2/ksp4)
     styled to match the OLD pre-revision Fig. S6, for direct visual
     comparison against the original image.
  2. The FULL new 9-point joint-grid plot (all n x ksp combinations),
     which is what this revision's Fig. S6 should actually show.
"""
import os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from nksp_common import (
    MODEL_NASH_PATH, MODEL_NAFL_PATH, OUT_DIR, CORE_OUTPUTS, KI,
    SILYMARIN_TARGETS, load_model, make_y0, parameterize_nksp,
    build_ki_source, make_mod_from_source, run_sim, pct_reduction,
    output_active_index,
)

# Match the OLD figure's legend order/labels and default matplotlib
# color cycle (C0=blue, C1=orange, C2=green) -- confirmed by eye against
# the uploaded original image.
LABELS = {'P_Hepatocyte_injury': 'Hepatocyte_injury',
          'P_Cell_death': 'Cell_death',
          'P_Inflammation': 'Inflammation'}
COLORS = {'P_Hepatocyte_injury': 'C0', 'P_Cell_death': 'C1', 'P_Inflammation': 'C2'}
PLOT_ORDER = ['P_Hepatocyte_injury', 'P_Cell_death', 'P_Inflammation']

FIVE_POINT_COMBOS = [(1, 2), (2, 2), (4, 2), (2, 1), (2, 4)]
NINE_POINT_COMBOS = [(n, ksp) for n in [1, 2, 4] for ksp in [1, 2, 4]]


def compute_all():
    with open(MODEL_NASH_PATH, encoding='utf-8') as f:
        base_nash = f.read()
    with open(MODEL_NAFL_PATH, encoding='utf-8') as f:
        base_nafl = f.read()

    results = {}
    for n, ksp in NINE_POINT_COMBOS:
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
        for pw in CORE_OUTPUTS:
            i = output_active_index(m_nash, pw)
            nash_red = pct_reduction(tn1[:, i], tn0[:, i])
            nafl_red = pct_reduction(tf1[:, i], tf0[:, i])
            ratio = nafl_red / nash_red if nash_red > 0 else np.nan
            results[(n, ksp, pw)] = ratio
        print(f'  computed n={n} ksp={ksp}')
    return results


def plot_subset(results, combos, out_path, title, mark_reversal=False):
    fig, ax = plt.subplots(figsize=(10.5, 6))
    x = np.arange(len(combos))
    for pw in PLOT_ORDER:
        y = [results[(n, ksp, pw)] for n, ksp in combos]
        ax.plot(x, y, marker='o', markersize=10, linewidth=2,
                color=COLORS[pw], label=LABELS[pw])
    ax.axhline(1.0, color='gray', ls='--', lw=1.2)
    ax.text(len(combos) - 1, 1.03, 'advantage = 1 (no stage effect)',
            ha='right', va='bottom', fontsize=9, color='gray')
    ax.set_xticks(x)
    ax.set_xticklabels([f'n{n}\nksp{ksp}' for n, ksp in combos])
    ax.set_ylabel('NAFL / NASH stage-dependent response')
    ax.set_title(title)
    # Legend placed OUTSIDE the axes (right margin) so it can never overlap
    # a data point, regardless of where the highest/lowest values fall.
    ax.legend(loc='upper left', bbox_to_anchor=(1.01, 1.0), borderaxespad=0.)
    ax.grid(alpha=0.2)

    if mark_reversal:
        for i, (n, ksp) in enumerate(combos):
            if (n, ksp) == (1, 1):
                y = results[(1, 1, 'P_Hepatocyte_injury')]
                ax.annotate(f'n=ksp=1\nHepatocyte_injury={y:.3f}\n(only reversal)',
                            xy=(i, y), xytext=(i + 0.4, y + 0.8),
                            fontsize=9, color=COLORS['P_Hepatocyte_injury'],
                            arrowprops=dict(arrowstyle='->', color=COLORS['P_Hepatocyte_injury']))

    fig.tight_layout()
    fig.savefig(out_path, dpi=300, bbox_inches='tight')
    plt.close(fig)
    print(f'  Saved: {out_path}')


def main():
    print('Computing all 9 (n, ksp) combinations from the verified model files...')
    results = compute_all()

    print('\n[1/2] 5-point subset (matches OLD pre-revision Fig. S6 exactly) ...')
    plot_subset(results, FIVE_POINT_COMBOS,
                os.path.join(OUT_DIR, 'SuppS6_5point_match_check.png'),
                'Parameter robustness of the stage-dependent advantage\n'
                '(5-point subset -- for comparison against the original figure)')

    print('\n[2/2] Full 9-point joint grid (this revision\'s actual Supp. Fig. S6) ...')
    plot_subset(results, NINE_POINT_COMBOS,
                os.path.join(OUT_DIR, 'SuppS6_9point_full_grid.png'),
                'Parameter robustness of the stage-dependent response\n'
                '(full 3\u00d73 joint n \u00d7 ksp grid)',
                mark_reversal=True)

    print('\nDone.')


if __name__ == '__main__':
    main()
