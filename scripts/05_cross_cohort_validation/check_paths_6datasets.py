#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
check_paths_6datasets.py
================================================================
Checks all paths required before running fig_ValidationCrossCohort_6datasets.py
and whether the files exist, instead of failing midway through the main script.

How to run:
  %run check_paths_6datasets.py
================================================================
"""
import os

# ================================================================
# Fully self-contained: all data paths resolve relative to this
# repository (no external folders required). Must exactly match the
# settings in fig_ValidationCrossCohort_6datasets.py.
# ================================================================
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
REPO_ROOT = os.path.dirname(os.path.dirname(SCRIPT_DIR))
MODEL_DIR = os.path.join(REPO_ROOT, 'data', 'models')
RATIOS_DIR = os.path.join(REPO_ROOT, 'data', 'ratios')
# ================================================================


def check(path, label):
    ok = os.path.exists(path)
    print(f'  {"OK  " if ok else "MISS"} {label}')
    if not ok:
        print(f'        -> {path}')
    return ok


def main():
    print('=' * 70)
    print('[1] Checking the ODE model files inside data/models/')
    print('=' * 70)
    check(MODEL_DIR, f'data/models/ folder: {MODEL_DIR}')
    all_ok = True
    for fn in ['ode_model_pydeseq2_NASH_vs_Normal_v10_mean.py',
               'ode_model_pydeseq2_NAFL_vs_Normal_v10_mean.py',
               'ode_model_pydeseq2_Obese_vs_Normal_v10_mean.py']:
        ok = check(os.path.join(MODEL_DIR, fn), fn)
        all_ok = all_ok and (ok or fn.startswith('ode_model_pydeseq2_Obese'))
    if not all_ok:
        print('\n  [!] Both the NASH and NAFL model files are required (a missing Obese file only skips the Obese curve and does not affect the rest).')

    print('\n' + '=' * 70)
    print('[2] Checking the GSE48452 / GSE89632 ratio files inside data/ratios/')
    print('=' * 70)
    check(RATIOS_DIR, f'data/ratios/ folder: {RATIOS_DIR}')
    check(os.path.join(RATIOS_DIR, 'GSE48452_node_initial_ratios.xlsx'),
          'GSE48452_node_initial_ratios.xlsx')
    check(os.path.join(RATIOS_DIR, 'GSE89632_node_initial_ratios.xlsx'),
          'GSE89632_node_initial_ratios.xlsx')

    print('\n' + '=' * 70)
    print('[3] Checking the GSE130970 / GSE162694 / GSE213621 ratio files inside data/ratios/')
    print('=' * 70)
    for ds in ['GSE130970', 'GSE162694', 'GSE213621']:
        check(os.path.join(RATIOS_DIR, f'{ds}_node_initial_ratios.xlsx'),
              f'{ds}_node_initial_ratios.xlsx')

    print('\n' + '=' * 70)
    print('[4] Checking whether the current working directory is where fig_ValidationCrossCohort_6datasets.py resides')
    print('=' * 70)
    print(f'  Current working directory: {os.getcwd()}')
    main_script = os.path.join(os.getcwd(), 'fig_ValidationCrossCohort_6datasets.py')
    check(main_script, 'fig_ValidationCrossCohort_6datasets.py (should be in the same folder as this check script)')

    print('\nDone. If everything above shows OK (Obese model optional), you can run:')
    print('  %run fig_ValidationCrossCohort_6datasets.py')


if __name__ == '__main__':
    main()
