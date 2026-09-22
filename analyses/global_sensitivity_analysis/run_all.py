#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
One-click runner: open this file in Spyder and press F5.
Runs run_morris_gsa (Morris global sensitivity) then step2_topology_comparison
(Spearman correlation against out-degree / distance-to-output), in sequence.

Edit run_morris_gsa.R_TRAJ / N_WORKERS at the top of run_morris_gsa.py before
running if you want a different sample size or core count.
"""
import time
import run_morris_gsa as gsa
import step2_topology_comparison as topo

if __name__ == '__main__':
    t0 = time.time()

    print('=' * 70)
    print(f'STEP 1/2: Morris global sensitivity (R_TRAJ={gsa.R_TRAJ}, '
          f'{gsa.R_TRAJ * 59} evaluations planned)')
    print('=' * 70)
    problem, names, baseline, param_values = gsa.build_problem_and_samples()
    Y = gsa.run_evaluations(param_values)
    gsa.analyze_and_save(problem, names, Y)

    print('\n' + '=' * 70)
    print('STEP 2/2: topology comparison')
    print('=' * 70)
    topo.main()

    print(f'\nAll done in {(time.time() - t0) / 60:.1f} minutes.')
