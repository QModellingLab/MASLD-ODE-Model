#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Compute the TRUE silymarin-target result, and/or batches of the two null
distributions ("random_target", "matched_control") for the R1-7 / R2-4 /
R3-5 control analysis.

Usage:
    python step1_run_batch.py true
        -> computes the real 8-silymarin-target result once, saves
           outputs/random_target_control_true_result.json (sanity-checked
           against manuscript Table 4 GSE126848 numbers: P_Hepatocyte_injury
           NAFL=34.31% NASH=16.30% ratio=2.10; P_Cell_death NAFL=13.85%
           NASH=5.67% ratio=2.44; P_Inflammation NAFL=9.01% NASH=2.97%
           ratio=3.03).

    python step1_run_batch.py random_target <start> <end>
        -> runs replicates [start, end) of the "any 8 of 58 nodes" null.
           Resumable: replicate indices already present in
           outputs/random_target_results.jsonl are skipped.

    python step1_run_batch.py matched_control <start> <end>
        -> same, but drawing from the pool that excludes both the real
           silymarin targets AND the top Fig.6 sensitivity drivers
           (outputs/matched_control_results.jsonl).

NOTE ON FILENAME (2026-09-11): this script's "true" output is named
random_target_control_true_result.json (previously true_result.json) --
the generic name collided with edge_scramble_control's own true_result.json
when both were uploaded into the same flat Project knowledge folder. If
you have an OLD outputs/true_result.json from a prior run of this script,
rename or delete it; step2_analyze.py only looks for the new name.

Each replicate only costs 2 ODE solves (NASH-drug, NAFL-drug); the two
no-drug baselines are solved once at the start of the script and reused.
Expect roughly 1-3 s/replicate depending on your machine.

Recommended total N per null: 2000 (run in batches, e.g.
"python step1_run_batch.py random_target 0 2000" all at once, or split
into several calls with different <start>/<end> ranges if you'd rather
checkpoint progress).
"""
import sys, os, json, time
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from rtc_common import (
    MODEL_NASH_PATH, MODEL_NAFL_PATH, OUT_DIR, CORE_OUTPUTS, KI, N_TARGETS,
    SILYMARIN_TARGETS, EXCLUDE_FOR_MATCHED_CONTROL,
    load_model, eligible_pool, make_y0, build_ki_source, make_mod_from_source,
    run_sim, pct_reduction, output_active_index, already_done,
)

BASE_SEED = {'random_target': 30260908, 'matched_control': 40260908}


def run_true():
    m_nash = load_model(MODEL_NASH_PATH)
    m_nafl = load_model(MODEL_NAFL_PATH)
    y0_nash = make_y0(m_nash)
    y0_nafl = make_y0(m_nafl)

    with open(MODEL_NASH_PATH, encoding='utf-8') as fh:
        src_nash = fh.read()
    with open(MODEL_NAFL_PATH, encoding='utf-8') as fh:
        src_nafl = fh.read()

    drug_nash_mod = make_mod_from_source(
        build_ki_source(src_nash, m_nash, SILYMARIN_TARGETS, KI), 'true_drug_nash')
    drug_nafl_mod = make_mod_from_source(
        build_ki_source(src_nafl, m_nafl, SILYMARIN_TARGETS, KI), 'true_drug_nafl')

    traj_nash_nodrug = run_sim(m_nash, y0_nash)
    traj_nash_drug = run_sim(drug_nash_mod, y0_nash)
    traj_nafl_nodrug = run_sim(m_nafl, y0_nafl)
    traj_nafl_drug = run_sim(drug_nafl_mod, y0_nafl)

    result = {'targets': SILYMARIN_TARGETS, 'ki': KI}
    print('TRUE silymarin-target result (sanity check vs manuscript Table 4):')
    for pw in CORE_OUTPUTS:
        i_nash = output_active_index(m_nash, pw)
        i_nafl = output_active_index(m_nafl, pw)
        nash_red = pct_reduction(traj_nash_drug[:, i_nash], traj_nash_nodrug[:, i_nash])
        nafl_red = pct_reduction(traj_nafl_drug[:, i_nafl], traj_nafl_nodrug[:, i_nafl])
        ratio = nafl_red / nash_red if nash_red > 0 else np.nan
        result[pw] = {'NASH_red': nash_red, 'NAFL_red': nafl_red, 'ratio': ratio,
                       'diff_pp': nafl_red - nash_red}
        print(f'  {pw:22s} NAFL={nafl_red:.4f}%  NASH={nash_red:.4f}%  '
              f'ratio={ratio:.4f}  diff={nafl_red-nash_red:.4f}pp')

    with open(os.path.join(OUT_DIR, 'random_target_control_true_result.json'), 'w') as fh:
        json.dump(result, fh, indent=2)
    print(f"\nSaved {os.path.join(OUT_DIR, 'random_target_control_true_result.json')}")


def run_null_batch(mode, start, end):
    assert mode in ('random_target', 'matched_control')
    out_path = os.path.join(OUT_DIR, f'{mode}_results.jsonl')

    m_nash = load_model(MODEL_NASH_PATH)
    m_nafl = load_model(MODEL_NAFL_PATH)
    y0_nash = make_y0(m_nash)
    y0_nafl = make_y0(m_nafl)

    with open(MODEL_NASH_PATH, encoding='utf-8') as fh:
        src_nash = fh.read()
    with open(MODEL_NAFL_PATH, encoding='utf-8') as fh:
        src_nafl = fh.read()

    print('[Setup] computing no-drug baselines once (reused for every replicate)...')
    traj_nash_nodrug = run_sim(m_nash, y0_nash)
    traj_nafl_nodrug = run_sim(m_nafl, y0_nafl)
    idx_nash = {pw: output_active_index(m_nash, pw) for pw in CORE_OUTPUTS}
    idx_nafl = {pw: output_active_index(m_nafl, pw) for pw in CORE_OUTPUTS}

    pool = eligible_pool(m_nash)
    if mode == 'matched_control':
        pool = [n for n in pool if n not in EXCLUDE_FOR_MATCHED_CONTROL]
    print(f'[Setup] mode={mode}  eligible pool size={len(pool)}  '
          f'(drawing {N_TARGETS} nodes/replicate)')

    done = already_done(out_path)
    todo = [i for i in range(start, end) if i not in done]
    print(f'[Resume] {len(done & set(range(start, end)))} replicates in [{start},{end}) '
          f'already present; running {len(todo)} new.')

    fh_out = open(out_path, 'a')
    t_start = time.time()
    for k, rep in enumerate(todo):
        t0 = time.time()
        rng = np.random.default_rng(BASE_SEED[mode] + rep)
        targets = sorted(rng.choice(pool, size=N_TARGETS, replace=False).tolist())

        drug_nash_mod = make_mod_from_source(
            build_ki_source(src_nash, m_nash, targets, KI), f'rt_nash_{rep}')
        drug_nafl_mod = make_mod_from_source(
            build_ki_source(src_nafl, m_nafl, targets, KI), f'rt_nafl_{rep}')

        try:
            traj_nash_drug = run_sim(drug_nash_mod, y0_nash)
            traj_nafl_drug = run_sim(drug_nafl_mod, y0_nafl)
        except Exception as e:
            fh_out.write(json.dumps({'replicate': rep, 'targets': targets, 'error': str(e)}) + '\n')
            fh_out.flush()
            continue

        rec = {'replicate': rep, 'targets': targets}
        bad = False
        for pw in CORE_OUTPUTS:
            nash_red = pct_reduction(traj_nash_drug[:, idx_nash[pw]], traj_nash_nodrug[:, idx_nash[pw]])
            nafl_red = pct_reduction(traj_nafl_drug[:, idx_nafl[pw]], traj_nafl_nodrug[:, idx_nafl[pw]])
            if not (np.isfinite(nash_red) and np.isfinite(nafl_red)):
                bad = True
            ratio = nafl_red / nash_red if (np.isfinite(nash_red) and nash_red > 0) else None
            rec[f'{pw}_NASH_red'] = round(float(nash_red), 4) if np.isfinite(nash_red) else None
            rec[f'{pw}_NAFL_red'] = round(float(nafl_red), 4) if np.isfinite(nafl_red) else None
            rec[f'{pw}_ratio'] = round(float(ratio), 4) if ratio is not None else None
            rec[f'{pw}_diff_pp'] = round(float(nafl_red - nash_red), 4) if (np.isfinite(nash_red) and np.isfinite(nafl_red)) else None
        rec['nan_flag'] = bad

        fh_out.write(json.dumps(rec) + '\n')
        fh_out.flush()

        if (k + 1) % 10 == 0 or (k + 1) == len(todo):
            elapsed = time.time() - t_start
            rate = elapsed / (k + 1)
            eta = rate * (len(todo) - (k + 1))
            print(f'  {k+1}/{len(todo)} new replicates done, '
                  f'elapsed={elapsed:.0f}s, eta={eta:.0f}s')

    fh_out.close()
    print(f'\nDone: {mode} range [{start}, {end}) -- '
          f'computed {len(todo)} new replicates, skipped {len(done & set(range(start, end)))} '
          f'already present. Total in {out_path}: '
          f'{sum(1 for _ in open(out_path))}')


# ============================================================
# DEFAULT PARAMETERS -- edit these and just press F5 in Spyder.
# Command-line arguments (if given) always override these defaults,
# so "python step1_run_batch.py random_target 0 500" still works too.
#
#   DEFAULT_MODE  : 'true' | 'random_target' | 'matched_control'
#   DEFAULT_START : starting replicate index (ignored for 'true')
#   DEFAULT_END   : ending replicate index, exclusive (ignored for 'true')
#
# Typical workflow: run once with DEFAULT_MODE='true', then repeatedly
# edit DEFAULT_START/DEFAULT_END (e.g. 0->500, 500->1000, ...) for each
# of 'random_target' and 'matched_control', pressing F5 each time.
# ============================================================
DEFAULT_MODE = 'true'
DEFAULT_START = 0
DEFAULT_END = 500


if __name__ == '__main__':
    if len(sys.argv) >= 2:
        mode = sys.argv[1]
        start = int(sys.argv[2]) if len(sys.argv) > 2 else DEFAULT_START
        end = int(sys.argv[3]) if len(sys.argv) > 3 else DEFAULT_END
    else:
        mode, start, end = DEFAULT_MODE, DEFAULT_START, DEFAULT_END
        print(f'[No command-line args given -- using DEFAULT_MODE={mode!r}, '
              f'DEFAULT_START={start}, DEFAULT_END={end} from the top of this file]')

    if mode == 'true':
        run_true()
    elif mode in ('random_target', 'matched_control'):
        run_null_batch(mode, start, end)
    else:
        print(f'Unknown mode: {mode!r}. Must be one of: true, random_target, matched_control')
        sys.exit(1)
