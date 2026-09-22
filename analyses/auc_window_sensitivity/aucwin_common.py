#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Shared utilities for auc_window_sensitivity (R1-5 / R3 minor-1)."""
import os, types, importlib.util
import numpy as np
from scipy.integrate import odeint

HERE = os.path.dirname(os.path.abspath(__file__))
MODEL_NASH_PATH = os.path.join(HERE, 'models', 'ode_model_pydeseq2_NASH_vs_Normal_v10_mean.py')
MODEL_NAFL_PATH = os.path.join(HERE, 'models', 'ode_model_pydeseq2_NAFL_vs_Normal_v10_mean.py')
OUT_DIR = os.path.join(HERE, 'outputs')
os.makedirs(OUT_DIR, exist_ok=True)

# Simulate ONCE over a long horizon (0.1h resolution, same density as the
# manuscript's 300h/3001-point run); AUC windows are then computed as
# sub-slices of this single long trajectory -- no need to re-simulate per
# window.
T_MAX_LONG = 1000.0
N_POINTS_LONG = 10001
T_LONG = np.linspace(0, T_MAX_LONG, N_POINTS_LONG)

CORE_OUTPUTS = ['P_Cell_death', 'P_Hepatocyte_injury', 'P_Inflammation']
KI = 0.3
SILYMARIN_TARGETS = ['CASP3', 'CASP7', 'CASP8', 'CYP2E1', 'IL_8', 'TNFa', 'NF_kB', 'TGF_b1']


def load_model(path):
    spec = importlib.util.spec_from_file_location(
        os.path.splitext(os.path.basename(path))[0] + '_aucwin', path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def make_y0(model):
    y0 = np.zeros(len(model.STATE_VARS))
    for i, sv in enumerate(model.STATE_VARS):
        if sv.endswith('_inactive'):
            node = sv[:-len('_inactive')]
            if node in model.GENE_METADATA:
                y0[i] = model.GENE_METADATA[node]['initial_ratio']
            elif node in model.PATHWAY_METADATA:
                y0[i] = model.PATHWAY_METADATA[node]['initial_level']
            else:
                raise KeyError(node)
    return y0


def build_ki_source(base_src, model, target_nodes, ki):
    idx_pairs = [(model.STATE_VARS.index(n + '_active'),
                  model.STATE_VARS.index(n + '_inactive')) for n in target_nodes]
    lines = base_src.split('\n')
    for ai, ii in idx_pairs:
        for idx in (ii, ai):
            for li, line in enumerate(lines):
                if line.lstrip().startswith(f'dydt[{idx}]'):
                    eq = line.find('=')
                    rhs = line[eq + 1:].strip()
                    if rhs.startswith('-('):
                        sign, bs = '-', 1
                    elif rhs.startswith('('):
                        sign, bs = '', 0
                    else:
                        continue
                    depth, ep = 0, None
                    for p, ch in enumerate(rhs[bs:], start=bs):
                        if ch == '(':
                            depth += 1
                        elif ch == ')':
                            depth -= 1
                            if depth == 0:
                                ep = p
                                break
                    if ep is None:
                        continue
                    act = rhs[bs:ep + 1]
                    rem = rhs[ep + 1:]
                    ind = line[:len(line) - len(line.lstrip())]
                    lines[li] = f"{ind}dydt[{idx}] = {sign}(({act})*{ki:.4f}){rem}"
                    break
    return '\n'.join(lines)


def make_mod_from_source(src, name):
    mod = types.ModuleType(name)
    mod.__file__ = name
    exec(compile(src, name, 'exec'), mod.__dict__)
    return mod


def run_sim_long(model, y0):
    """Single long-horizon solve (0-1000h); reused for every AUC window."""
    return odeint(model.ode_system, y0, T_LONG, args=(None,), mxstep=20000)


def auc_window(traj_col, t_window):
    """Trapezoidal AUC of traj_col restricted to [0, t_window] hours."""
    mask = T_LONG <= t_window
    f = getattr(np, 'trapezoid', None) or np.trapz
    return float(f(traj_col[mask], T_LONG[mask]))


def pct_reduction_window(traj_drug_col, traj_nodrug_col, t_window):
    a0 = auc_window(traj_nodrug_col, t_window)
    a1 = auc_window(traj_drug_col, t_window)
    return (a0 - a1) / a0 * 100 if a0 > 0 else np.nan


def output_active_index(model, pw):
    return model.STATE_VARS.index(pw + '_active')
