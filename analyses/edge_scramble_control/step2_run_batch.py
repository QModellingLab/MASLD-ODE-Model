#!/usr/bin/env python3
"""
step2_run_batch.py
===================
Runs replicates of the silymarin intervention (k_i = 0.3) on randomized
network topologies, for either of two null models, and appends results
incrementally to a .jsonl file (safe to interrupt and resume -- each line
is one completed replicate).

Two null models (see README for rationale):
  scramble  -- degree-preserving double edge-swap (Maslov & Sneppen, 2002).
               Preserves each node's in-degree and out-degree exactly.
  er        -- Erdos-Renyi random rewiring. Does NOT preserve hsa04932's
               degree sequence; a genuinely different null-model class,
               used to check whether results depend on matching the real
               degree distribution.

Usage (run from this folder):
    python step2_run_batch.py true              # compute the real-topology reference once (required first)
    python step2_run_batch.py scramble 0 200
    python step2_run_batch.py scramble 0 400        # extend / resume -- already-completed
                                                     # replicates (0-199) are auto-detected and
                                                     # skipped, only 200-399 are computed
    python step2_run_batch.py er 0 200
    python step2_run_batch.py er 0 400

NOTE ON FILENAME (2026-09-11): the "true" mode's output is named
outputs/edge_scramble_true_result.json (previously outputs/true_result.json)
-- the generic name collided with random_target_control's own
true_result.json when both were uploaded into the same flat Project
knowledge folder. If you have an OLD outputs/true_result.json from a
prior run of this script, rename or delete it; step3_analyze.py only
looks for the new name.

Each call is independent and safe to run in a fresh terminal session, and
safe to re-run with the SAME or an OVERLAPPING range as a previous call:
the script reads whatever is already in the .jsonl output file, and skips
any replicate index already present there, so re-running the same command
after an interrupted run (or just to be safe) never creates duplicate
rows. The underlying RNG stream is still fixed-seed and advanced
deterministically through every replicate index from 0 up to n_end
regardless of what's already done, so the *topology* generated for a
given replicate index is always identical no matter how many times, or
in what order, you've called this script.

Runtime: roughly 1.5-2.5 seconds per NEW replicate on a typical laptop (4
ODE integrations per replicate: NASH baseline, NASH+drug, NAFL baseline,
NAFL+drug). 200 new replicates ~= 5-8 minutes; already-completed
replicates are skipped almost instantly.
"""
import sys
import time
import random
import json

import numpy as np
if not hasattr(np, 'trapezoid'):
    np.trapezoid = np.trapz
from scipy.integrate import odeint

from model_loader import load_models
from generic_model import (
    build_generic_ode, INH_EDGES, ACT_EDGES,
    degree_preserving_scramble, erdos_renyi_edges,
    SILYMARIN_TARGETS, CORE_OUTPUTS, ALL_MOLECULAR_NODES,
)

KI = 0.3
KI_MAP = {n: KI for n in SILYMARIN_TARGETS}
N_SWAPS_PER_REPLICATE = len(ACT_EDGES) * 20  # Milo et al. recommend ~10-100x edge count


def run_replicate(nash, nafl, SV, active_map, t, act_edges):
    ode_base = build_generic_ode(SV, act_edges, INH_EDGES)
    ode_drug = build_generic_ode(SV, act_edges, INH_EDGES, ki_multipliers=KI_MAP)
    yb_nash = odeint(ode_base, nash.Y0, t, args=(nash.PARAMS,), mxstep=5000)
    yd_nash = odeint(ode_drug, nash.Y0, t, args=(nash.PARAMS,), mxstep=5000)
    yb_nafl = odeint(ode_base, nafl.Y0, t, args=(nafl.PARAMS,), mxstep=5000)
    yd_nafl = odeint(ode_drug, nafl.Y0, t, args=(nafl.PARAMS,), mxstep=5000)
    result = {}
    for node in CORE_OUTPUTS:
        pi = active_map[node]
        a0n, a1n = np.trapezoid(yb_nash[:, pi], t), np.trapezoid(yd_nash[:, pi], t)
        a0f, a1f = np.trapezoid(yb_nafl[:, pi], t), np.trapezoid(yd_nafl[:, pi], t)
        rn = 100 * (a0n - a1n) / a0n if a0n else float('nan')
        rf = 100 * (a0f - a1f) / a0f if a0f else float('nan')
        result[node] = {'nafl_red': rf, 'nash_red': rn,
                         'ratio': (rf / rn if rn not in (0, float('nan')) else float('nan'))}
    return result


def main():
    if len(sys.argv) == 2 and sys.argv[1] == 'true':
        # Special mode: compute the TRUE (real hsa04932) topology's result
        # once, through the identical pipeline used for null replicates,
        # so downstream comparisons are apples-to-apples (same time grid,
        # same mxstep, same everything except the edge list).
        nash, nafl = load_models()
        SV = nash.STATE_VARS
        active_map = {s[:-7]: i for i, s in enumerate(SV) if s.endswith('_active')}
        t = np.linspace(0, 300, 1501)
        result = run_replicate(nash, nafl, SV, active_map, t, ACT_EDGES)
        import os
        os.makedirs('outputs', exist_ok=True)
        with open('outputs/edge_scramble_true_result.json', 'w') as f:
            json.dump(result, f, indent=2)
        print("TRUE topology result (identical pipeline as null replicates):")
        for node, r in result.items():
            print(f"  {node:<22} NAFL={r['nafl_red']:.4f}%  NASH={r['nash_red']:.4f}%  "
                  f"ratio={r['ratio']:.4f}  diff={r['nafl_red']-r['nash_red']:.4f}pp")
        print("\nSaved outputs/edge_scramble_true_result.json")
        return

    if len(sys.argv) != 4 or sys.argv[1] not in ('scramble', 'er'):
        print(__doc__)
        sys.exit(1)
    mode, n_start, n_end = sys.argv[1], int(sys.argv[2]), int(sys.argv[3])
    outfile = f'outputs/{mode}_results.jsonl'
    seed = 42 if mode == 'scramble' else 123

    import os
    os.makedirs('outputs', exist_ok=True)

    # Auto-detect already-completed replicate indices, so re-running the
    # same (or an overlapping) range never creates duplicate rows.
    completed = set()
    if os.path.exists(outfile):
        with open(outfile) as f:
            for line in f:
                line = line.strip()
                if not line:
                    continue
                try:
                    completed.add(json.loads(line)['rep'])
                except (json.JSONDecodeError, KeyError):
                    pass  # ignore any corrupted trailing line
    n_already = len([r for r in completed if n_start <= r < n_end])
    print(f"Found {len(completed)} completed replicates in {outfile} "
          f"({n_already} within the requested range [{n_start}, {n_end})); "
          f"these will be skipped.")

    nash, nafl = load_models()
    SV = nash.STATE_VARS
    active_map = {s[:-7]: i for i, s in enumerate(SV) if s.endswith('_active')}
    t = np.linspace(0, 300, 1501)

    rng = random.Random(seed)
    t_start = time.time()
    last_report_time = 0.0
    n_computed = 0
    n_to_compute = (n_end - n_start) - n_already
    with open(outfile, 'a') as f:
        for rep in range(0, n_end):
            # ALWAYS advance the RNG for every replicate index from 0
            # onward, regardless of whether it's in range or already
            # done, so the topology for a given replicate index never
            # depends on how the run was split across calls.
            if mode == 'scramble':
                edges, _ = degree_preserving_scramble(ACT_EDGES, N_SWAPS_PER_REPLICATE, rng)
            else:
                target_pool = ALL_MOLECULAR_NODES + sorted(set(t_ for _, t_ in ACT_EDGES if t_.startswith('P_')))
                edges = erdos_renyi_edges(len(ACT_EDGES), ALL_MOLECULAR_NODES, target_pool, rng)

            if rep < n_start or rep >= n_end:
                continue
            if rep in completed:
                continue

            result = run_replicate(nash, nafl, SV, active_map, t, edges)
            f.write(json.dumps({'rep': rep, 'result': result}) + '\n')
            f.flush()
            n_computed += 1

            elapsed = time.time() - t_start
            # Report whenever EITHER: 15 seconds have passed since the last
            # report, OR every 10 replicates, OR we just crossed a new 10%
            # milestone, OR this is the very last one -- whichever comes
            # first. This guarantees a visible update at least every ~15s
            # even if per-replicate timing is slower than expected, and
            # guarantees at least ~10 updates total regardless of n.
            pct = 100 * n_computed / n_to_compute if n_to_compute else 100
            prev_pct = 100 * (n_computed - 1) / n_to_compute if n_to_compute else 100
            crossed_10pct = int(pct // 10) > int(prev_pct // 10)
            time_to_report = (elapsed - last_report_time) >= 15
            if (n_computed % 10 == 0 or n_computed == n_to_compute
                    or crossed_10pct or time_to_report):
                eta = elapsed / n_computed * (n_to_compute - n_computed)
                print(f"  {n_computed}/{n_to_compute} ({pct:.0f}%) new replicates done, "
                      f"elapsed={elapsed:.0f}s, eta={eta:.0f}s", flush=True)
                last_report_time = elapsed

    print(f"Done: {mode} range [{n_start}, {n_end}) -- computed {n_computed} new replicates "
          f"({time.time()-t_start:.0f}s), skipped {n_already} already present. "
          f"Total in {outfile}: {len(completed) + n_computed}")


if __name__ == '__main__':
    main()
