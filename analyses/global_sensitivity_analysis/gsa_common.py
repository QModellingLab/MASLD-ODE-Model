#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Shared utilities for the global_sensitivity_analysis package (R1-6).
Module name is deliberately unique (gsa_common, not common) to avoid the
cross-folder Python module-name collision that occurred between
random_target_control/common.py and static_hill_curve_control/common.py.
"""
import os
import importlib.util
import numpy as np
from scipy.integrate import odeint

HERE = os.path.dirname(os.path.abspath(__file__))
MODEL_NASH_PATH = os.path.join(HERE, 'models', 'ode_model_pydeseq2_NASH_vs_Normal_v10_mean.py')
OUT_DIR = os.path.join(HERE, 'outputs')
os.makedirs(OUT_DIR, exist_ok=True)

T_MAX, N_POINTS = 300.0, 3001
T = np.linspace(0, T_MAX, N_POINTS)
CORE_OUTPUTS = ['P_Cell_death', 'P_Hepatocyte_injury', 'P_Inflammation']

# Multiplicative bounds around each node's true NASH-baseline initial ratio,
# defining the range explored by the global sensitivity analysis (much wider
# than the manuscript's local +/-1% perturbation, and varying ALL nodes
# simultaneously per Morris trajectory rather than one at a time).
BOUNDS_FACTOR_LO = 0.3
BOUNDS_FACTOR_HI = 3.0


def _trapz(y):
    f = getattr(np, 'trapezoid', None) or np.trapz
    return float(f(y, T))


def load_model(path=MODEL_NASH_PATH):
    spec = importlib.util.spec_from_file_location(
        os.path.splitext(os.path.basename(path))[0] + '_gsa', path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def node_names(model):
    """The 58 molecular node names, in a fixed sorted order (defines the
    parameter order used throughout the GSA)."""
    return sorted(model.GENE_METADATA.keys())


def baseline_vector(model, names):
    return np.array([model.GENE_METADATA[n]['initial_ratio'] for n in names])


def build_bounds(baseline, lo=BOUNDS_FACTOR_LO, hi=BOUNDS_FACTOR_HI):
    return [[max(1e-6, b * lo), b * hi] for b in baseline]


def make_y0_from_params(model, names, params):
    """params: array aligned with `names` (molecular-node initial ratios).
    Active forms stay 0; P_* pathway-output inactive forms stay at their
    PATHWAY_METADATA initial_level (100.0), unaffected by the GSA."""
    y0 = np.zeros(len(model.STATE_VARS))
    param_map = dict(zip(names, params))
    for i, sv in enumerate(model.STATE_VARS):
        if sv.endswith('_inactive'):
            node = sv[:-len('_inactive')]
            if node in param_map:
                y0[i] = param_map[node]
            elif node in model.PATHWAY_METADATA:
                y0[i] = model.PATHWAY_METADATA[node]['initial_level']
            else:
                raise KeyError(f'Unrecognized node {node!r}')
    return y0


def run_sim(model, y0):
    return odeint(model.ode_system, y0, T, args=(None,), mxstep=10000)


def output_active_index(model, pw):
    return model.STATE_VARS.index(pw + '_active')


# ------------------------------------------------------------------
# Worker-process globals (populated once per process by _init_worker,
# so the model is loaded from disk only once per worker, not once per
# evaluated parameter set).
# ------------------------------------------------------------------
_WORKER_MODEL = None
_WORKER_NAMES = None
_WORKER_IDX = None


def _init_worker():
    global _WORKER_MODEL, _WORKER_NAMES, _WORKER_IDX
    _WORKER_MODEL = load_model()
    _WORKER_NAMES = node_names(_WORKER_MODEL)
    _WORKER_IDX = {pw: output_active_index(_WORKER_MODEL, pw) for pw in CORE_OUTPUTS}


def eval_sample(params):
    """Top-level, picklable worker function: params -> (AUC_Cell_death,
    AUC_Hepatocyte_injury, AUC_Inflammation). Requires _init_worker() to
    have been called in this process first (done automatically by
    ProcessPoolExecutor's `initializer=` argument)."""
    global _WORKER_MODEL, _WORKER_NAMES, _WORKER_IDX
    if _WORKER_MODEL is None:
        _init_worker()
    y0 = make_y0_from_params(_WORKER_MODEL, _WORKER_NAMES, params)
    traj = run_sim(_WORKER_MODEL, y0)
    return tuple(_trapz(traj[:, _WORKER_IDX[pw]]) for pw in CORE_OUTPUTS)
