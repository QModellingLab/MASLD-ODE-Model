#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Static Hill-curve control analysis (R3-1, second control).

Reviewer 3, Major Comment 1 asks for "a static Hill-curve calculation that
predicts the effect without integrating the ODEs" to test whether the
stage-dependent (NAFL > NASH) silymarin advantage is close to algebraically
guaranteed by (a) Hill-saturation kinetics and (b) the disease-stage-
dependent initial conditions alone -- i.e. a receptor-reserve / "law of
initial value" artifact of the model's construction rather than a genuine,
non-trivial network property.

This script provides two independent, non-AUC, non-time-course analyses:

PART A -- Network steady-state (root-finding, not time-integration)
  For NASH and NAFL, no-drug and drug (ki=0.3, real 8 silymarin targets),
  the TRUE equilibrium of the full 70-node network is found by root-
  finding dy/dt = 0 (scipy.optimize.fsolve/root) rather than by choosing
  an integration window. A short ODE run is used only to generate a good
  starting guess for the root-finder (a standard numerical convenience;
  the reported steady state is verified by checking the residual
  ||dy/dt|| at the solution, not by the integration window). The steady-
  state-based %reduction and NAFL/NASH ratio for the 3 core outputs are
  reported alongside the manuscript's AUC-based values, so the two
  metrics -- one window-dependent (AUC, 0-300h) and one window-free
  (steady state) -- can be compared directly.

PART B -- Single-reaction algebraic Hill/receptor-reserve calculation
  For each core output's direct upstream driver(s) (Table S3: CASP3/CASP7
  for P_Hepatocyte_injury; FasL/TNFa for P_Cell_death; IL_8/TGF_b1 for
  P_Inflammation), the Hill-saturation fraction h(x) = x^n/(ksp^n+x^n) is
  evaluated purely algebraically at the driver's own stage-specific
  initial ratio x (from GENE_METADATA; n=ksp=2.0 as in the manuscript).
  The "reserve" 1-h(x) is a direct receptor-reserve-style prediction of
  how much unsaturated capacity that driver's own activation kinetics
  have at each stage. This never touches the ODE at all -- it only uses
  the algebraic Hill formula and the two stage-specific input numbers.

Usage:
    python run_static_analysis.py
"""
import os
import numpy as np
import pandas as pd
from scipy.integrate import odeint
from scipy.optimize import fsolve

from hillcurve_common import (
    MODEL_NASH_PATH, MODEL_NAFL_PATH, OUT_DIR, CORE_OUTPUTS, KI,
    SILYMARIN_TARGETS, load_model, make_y0, build_ki_source,
    make_mod_from_source, output_active_index,
)

N_HILL, KSP_HILL = 2.0, 2.0  # matches manuscript's fixed Vmax=ksp=n=kcat=2.0

# Direct upstream drivers per core output (Table S3)
DIRECT_DRIVERS = {
    'P_Hepatocyte_injury': ['CASP3', 'CASP7'],
    'P_Cell_death': ['FasL', 'TNFa'],
    'P_Inflammation': ['IL_8', 'TGF_b1'],
}


# ---------------------------------------------------------------- PART A ---
def find_steady_state(model, y0_guess_seed, tmax_warmstart=5000.0):
    """Warm-start with a long (but still finite, cheap) ODE run, then root-
    find dy/dt=0 from that point. Returns (y_ss, residual_norm)."""
    t_warm = np.linspace(0, tmax_warmstart, 200)
    warm_traj = odeint(model.ode_system, y0_guess_seed, t_warm, args=(None,), mxstep=20000)
    guess = warm_traj[-1]

    def resid(y):
        return np.asarray(model.ode_system(y, 0.0, None))

    y_ss, info, ier, msg = fsolve(resid, guess, full_output=True, xtol=1e-10, maxfev=20000)
    res_norm = float(np.linalg.norm(resid(y_ss)))
    return y_ss, res_norm, ier


def part_a_steady_state():
    print('=' * 78)
    print('PART A: network steady state (root-finding, window-free)')
    print('=' * 78)

    m_nash = load_model(MODEL_NASH_PATH)
    m_nafl = load_model(MODEL_NAFL_PATH)
    y0_nash = make_y0(m_nash)
    y0_nafl = make_y0(m_nafl)

    with open(MODEL_NASH_PATH, encoding='utf-8') as fh:
        src_nash = fh.read()
    with open(MODEL_NAFL_PATH, encoding='utf-8') as fh:
        src_nafl = fh.read()
    drug_nash_mod = make_mod_from_source(
        build_ki_source(src_nash, m_nash, SILYMARIN_TARGETS, KI), 'ss_drug_nash')
    drug_nafl_mod = make_mod_from_source(
        build_ki_source(src_nafl, m_nafl, SILYMARIN_TARGETS, KI), 'ss_drug_nafl')

    combos = {
        ('NASH', 'nodrug'): (m_nash, y0_nash),
        ('NASH', 'drug'): (drug_nash_mod, y0_nash),
        ('NAFL', 'nodrug'): (m_nafl, y0_nafl),
        ('NAFL', 'drug'): (drug_nafl_mod, y0_nafl),
    }

    ss_active = {}  # (stage, cond) -> {output: active_ss_value}
    residuals = {}
    print('\nSolving for steady state (this may take ~10-30s per combo)...')
    for (stage, cond), (mod, y0) in combos.items():
        y_ss, res_norm, ier = find_steady_state(mod, y0)
        print(f'  {stage:5s} {cond:6s}  fsolve ier={ier}  residual_norm={res_norm:.3e}')
        ss_active[(stage, cond)] = {
            pw: y_ss[output_active_index(mod, pw)] for pw in CORE_OUTPUTS
        }
        residuals[(stage, cond)] = res_norm

    print('\n-- Steady-state vs AUC-based (manuscript Table 4) comparison --')
    print(f'{"Output":22s} {"SS %reduction NAFL":>19s} {"SS %reduction NASH":>19s} '
          f'{"SS ratio":>9s} {"AUC ratio (Table 4)":>20s}')
    auc_table4_ratio = {'P_Cell_death': 2.44, 'P_Hepatocyte_injury': 2.10, 'P_Inflammation': 3.03}
    results = {}
    rows_a = []
    for pw in CORE_OUTPUTS:
        p_nash_0 = ss_active[('NASH', 'nodrug')][pw]
        p_nash_1 = ss_active[('NASH', 'drug')][pw]
        p_nafl_0 = ss_active[('NAFL', 'nodrug')][pw]
        p_nafl_1 = ss_active[('NAFL', 'drug')][pw]
        red_nash = (p_nash_0 - p_nash_1) / p_nash_0 * 100 if p_nash_0 > 0 else np.nan
        red_nafl = (p_nafl_0 - p_nafl_1) / p_nafl_0 * 100 if p_nafl_0 > 0 else np.nan
        ss_ratio = red_nafl / red_nash if red_nash > 0 else np.nan
        results[pw] = dict(red_nash=red_nash, red_nafl=red_nafl, ss_ratio=ss_ratio,
                            p_nash_0=p_nash_0, p_nash_1=p_nash_1,
                            p_nafl_0=p_nafl_0, p_nafl_1=p_nafl_1)
        print(f'{pw:22s} {red_nafl:19.4f} {red_nash:19.4f} {ss_ratio:9.3f} '
              f'{auc_table4_ratio[pw]:20.2f}')
        rows_a.append({
            'output': pw,
            'SS_pct_reduction_NAFL': red_nafl,
            'SS_pct_reduction_NASH': red_nash,
            'SS_ratio_NAFL_over_NASH': ss_ratio,
            'AUC_ratio_manuscript_Table4': auc_table4_ratio[pw],
            'P_active_ss_NASH_nodrug': p_nash_0,
            'P_active_ss_NASH_drug': p_nash_1,
            'P_active_ss_NAFL_nodrug': p_nafl_0,
            'P_active_ss_NAFL_drug': p_nafl_1,
            'fsolve_residual_NASH_nodrug': residuals[('NASH', 'nodrug')],
            'fsolve_residual_NASH_drug': residuals[('NASH', 'drug')],
            'fsolve_residual_NAFL_nodrug': residuals[('NAFL', 'nodrug')],
            'fsolve_residual_NAFL_drug': residuals[('NAFL', 'drug')],
        })

    return results, rows_a


# ---------------------------------------------------------------- PART B ---
def hill(x, n=N_HILL, ksp=KSP_HILL):
    return x ** n / (ksp ** n + x ** n)


def part_b_algebraic_reserve():
    print('\n' + '=' * 78)
    print('PART B: single-reaction algebraic Hill / receptor-reserve check')
    print('=' * 78)
    print(f'(n = ksp = {N_HILL:.1f}, matching the manuscript\'s uniform kinetic constants)\n')

    m_nash = load_model(MODEL_NASH_PATH)
    m_nafl = load_model(MODEL_NAFL_PATH)

    print(f'{"Output":22s} {"Driver":8s} {"x0(NAFL)":>9s} {"x0(NASH)":>9s} '
          f'{"h(NAFL)":>8s} {"h(NASH)":>8s} {"reserve(NAFL)":>14s} {"reserve(NASH)":>14s} '
          f'{"NAFL>NASH reserve?":>19s}')
    rows_b = []
    for pw, drivers in DIRECT_DRIVERS.items():
        for drv in drivers:
            x_nafl = m_nafl.GENE_METADATA[drv]['initial_ratio']
            x_nash = m_nash.GENE_METADATA[drv]['initial_ratio']
            h_nafl, h_nash = hill(x_nafl), hill(x_nash)
            r_nafl, r_nash = 1 - h_nafl, 1 - h_nash
            flag = 'YES' if r_nafl > r_nash else 'no'
            print(f'{pw:22s} {drv:8s} {x_nafl:9.4f} {x_nash:9.4f} '
                  f'{h_nafl:8.4f} {h_nash:8.4f} {r_nafl:14.6f} {r_nash:14.6f} {flag:>19s}')
            rows_b.append({
                'output': pw, 'driver': drv,
                'x0_NAFL': x_nafl, 'x0_NASH': x_nash,
                'h_NAFL': h_nafl, 'h_NASH': h_nash,
                'reserve_NAFL': r_nafl, 'reserve_NASH': r_nash,
                'NAFL_gt_NASH_reserve': flag,
            })

    print('\nInterpretation guide:')
    print('  "reserve" = 1 - x0^n/(ksp^n+x0^n): the fraction of a driver\'s OWN')
    print('  activation Hill-curve that is still unsaturated at its stage-specific')
    print('  initial ratio x0. A larger reserve at NAFL than at NASH is the single-')
    print('  reaction, purely-algebraic form of the receptor-reserve argument the')
    print('  manuscript cites (Buchwald 2020; Wilder 1957). Because x0 for these')
    print('  drivers is only mildly larger in NASH than NAFL (both roughly the same')
    print('  order as ksp=2), the driver-level algebraic reserve difference is small')
    print('  -- if the full-network (Part A / manuscript AUC) NAFL/NASH advantage is')
    print('  much larger than this small single-reaction reserve difference would')
    print('  suggest, that indicates the cascading network dynamics amplify a small')
    print('  initial-condition difference, rather than the effect being a trivial,')
    print('  single-step algebraic consequence of Hill saturation alone.')

    return rows_b


def save_excel(rows_a, rows_b):
    os.makedirs(OUT_DIR, exist_ok=True)
    xlsx_path = os.path.join(OUT_DIR, 'static_hill_curve_summary.xlsx')
    df_a = pd.DataFrame(rows_a)
    df_b = pd.DataFrame(rows_b)
    with pd.ExcelWriter(xlsx_path, engine='openpyxl') as writer:
        df_a.to_excel(writer, sheet_name='PartA_steady_state', index=False)
        df_b.to_excel(writer, sheet_name='PartB_algebraic_reserve', index=False)
    print(f'\nExcel summary saved: {xlsx_path}')


if __name__ == '__main__':
    ss_results, rows_a = part_a_steady_state()
    rows_b = part_b_algebraic_reserve()
    save_excel(rows_a, rows_b)

    print('\n' + '=' * 78)
    print('Done. See the printed tables above, or outputs/static_hill_curve_summary.xlsx,')
    print('for the Response-to-Reviewers write-up.')
    print('=' * 78)
