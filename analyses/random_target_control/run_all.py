#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
One-click runner: executes the whole random_target_control pipeline in a
single call -- just open this file in Spyder and press F5.

Runs, in order:
  1. step1_run_batch.run_true()                          (~seconds)
  2. step1_run_batch.run_null_batch('random_target', ...) (~N_NULL seconds,
     since each replicate costs ~1s)
  3. step1_run_batch.run_null_batch('matched_control', ...)
  4. step2_analyze.main()                                 (~seconds)

Total wall-clock time for the default N_NULL=2000 is roughly
2000 x 2 x ~1s ~= 65-70 minutes. Progress prints every 10 replicates as
usual, so you can see it's alive while it runs.

To change the null-model size, just edit N_NULL below and press F5 again
-- already-completed replicate indices are skipped automatically (the
underlying jsonl files are resumable), so re-running after an interrupt
or after raising N_NULL only computes what's missing.
"""
import time
import step1_run_batch as s1
import step2_analyze as s2

# ============================================================
# EDIT THIS if you want a different null-model size (default 2000).
# ============================================================
N_NULL = 2000

if __name__ == '__main__':
    t_start = time.time()

    print('\n' + '=' * 70)
    print('STEP 1/4: true silymarin-target result')
    print('=' * 70)
    s1.run_true()

    print('\n' + '=' * 70)
    print(f'STEP 2/4: random_target null (0-{N_NULL})')
    print('=' * 70)
    s1.run_null_batch('random_target', 0, N_NULL)

    print('\n' + '=' * 70)
    print(f'STEP 3/4: matched_control null (0-{N_NULL})')
    print('=' * 70)
    s1.run_null_batch('matched_control', 0, N_NULL)

    print('\n' + '=' * 70)
    print('STEP 4/4: statistical analysis')
    print('=' * 70)
    s2.main()

    print(f'\nAll done in {(time.time() - t_start) / 60:.1f} minutes.')
