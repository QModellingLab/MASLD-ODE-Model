#!/usr/bin/env python3
"""
step1_validate.py
==================
MANDATORY first step. Confirms the generic, topology-driven ODE model
(generic_model.py, built purely from the 106-edge Supplementary Table S2
list) exactly reproduces the hardcoded NASH/NAFL model files' trajectories.

If this does not pass, do NOT trust any scramble/null-model results --
it would mean the edge list or gene-to-node mapping has an error.

Usage:
    python step1_validate.py
"""
import numpy as np
import numpy as np
if not hasattr(np, 'trapezoid'):
    np.trapezoid = np.trapz
from scipy.integrate import odeint

from model_loader import load_models
from generic_model import build_generic_ode, ACT_EDGES, INH_EDGES, SILYMARIN_TARGETS, CORE_OUTPUTS


def validate_condition(mod, label, ki_check=0.3):
    SV = mod.STATE_VARS
    active_map = {s[:-7]: i for i, s in enumerate(SV) if s.endswith('_active')}
    t = np.linspace(0, 300, 3001)

    generic_undrugged = build_generic_ode(SV, ACT_EDGES, INH_EDGES)
    y_hard = odeint(mod.ode_system, mod.Y0, t, args=(mod.PARAMS,), mxstep=10000)
    y_gen = odeint(generic_undrugged, mod.Y0, t, args=(mod.PARAMS,), mxstep=10000)
    diff = np.abs(y_hard - y_gen)
    rel = diff / (np.abs(y_hard) + 1e-9)
    print(f"\n[{label}] undrugged: max abs diff = {diff.max():.2e}, "
          f"max rel diff = {rel.max():.2e}  ({'PASS' if rel.max() < 1e-4 else 'FAIL -- CHECK EDGE LIST'})")

    ki_map = {n: ki_check for n in SILYMARIN_TARGETS}
    generic_drugged = build_generic_ode(SV, ACT_EDGES, INH_EDGES, ki_multipliers=ki_map)
    y_drug = odeint(generic_drugged, mod.Y0, t, args=(mod.PARAMS,), mxstep=10000)

    print(f"[{label}] silymarin intervention (k_i={ki_check}) AUC reduction (generic model):")
    for node in CORE_OUTPUTS:
        pi = active_map[node]
        auc0 = np.trapezoid(y_gen[:, pi], t)
        auc1 = np.trapezoid(y_drug[:, pi], t)
        red = 100 * (auc0 - auc1) / auc0
        print(f"    {node:<22} {red:6.2f}%   "
              f"(compare to manuscript Table 4, {label}, k_i={ki_check})")


def main():
    print("=" * 70)
    print("STEP 1: Validating generic topology-driven model against hardcoded files")
    print("=" * 70)
    nash, nafl = load_models()
    validate_condition(nash, 'NASH')
    validate_condition(nafl, 'NAFL')
    print("\nIf both show PASS and the AUC-reduction percentages above match your")
    print("manuscript's Table 4 (GSE126848 row), the generic model is validated --")
    print("proceed to step2_run_batch.py.")


if __name__ == '__main__':
    main()
