#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Shared utilities for alt_initial_conditions (R1-4)."""
import os, types, importlib.util
import numpy as np
from scipy.integrate import odeint

HERE = os.path.dirname(os.path.abspath(__file__))
MODEL_NASH_PATH = os.path.join(HERE, 'models', 'ode_model_pydeseq2_NASH_vs_Normal_v10_mean.py')
MODEL_NAFL_PATH = os.path.join(HERE, 'models', 'ode_model_pydeseq2_NAFL_vs_Normal_v10_mean.py')
OUT_DIR = os.path.join(HERE, 'outputs')
os.makedirs(OUT_DIR, exist_ok=True)

T_MAX, N_POINTS = 300.0, 3001
T = np.linspace(0, T_MAX, N_POINTS)
CORE_OUTPUTS = ['P_Cell_death', 'P_Hepatocyte_injury', 'P_Inflammation']
KI = 0.3
SILYMARIN_TARGETS = ['CASP3', 'CASP7', 'CASP8', 'CYP2E1', 'IL_8', 'TNFa', 'NF_kB', 'TGF_b1']


def _trapz(y):
    f = getattr(np, 'trapezoid', None) or np.trapz
    return float(f(y, T))


def load_model(path):
    spec = importlib.util.spec_from_file_location(
        os.path.splitext(os.path.basename(path))[0] + '_altic', path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def make_y0(model, active_fraction=0.0, p_output_inactive_level=100.0):
    """active_fraction: initialize each molecular node's ACTIVE form as this
    fraction of its INACTIVE form's initial ratio (manuscript default: 0.0,
    i.e. the whole network starts fully inactive). Addresses R1-4's question
    "why should all signaling nodes have no initial activity in established
    NAFL/NASH tissue?".

    p_output_inactive_level: the fixed inactive-form starting value for the
    12 P_* pathway-output nodes (manuscript default: 100.0). Addresses
    R1-4's question about this constant's justification.
    """
    y0 = np.zeros(len(model.STATE_VARS))
    for i, sv in enumerate(model.STATE_VARS):
        if sv.endswith('_inactive'):
            node = sv[:-len('_inactive')]
            if node in model.GENE_METADATA:
                y0[i] = model.GENE_METADATA[node]['initial_ratio']
            elif node in model.PATHWAY_METADATA:
                y0[i] = p_output_inactive_level
            else:
                raise KeyError(node)
        elif sv.endswith('_active'):
            node = sv[:-len('_active')]
            if node in model.GENE_METADATA:
                inactive_val = model.GENE_METADATA[node]['initial_ratio']
                y0[i] = active_fraction * inactive_val
            # P_* active forms always start at 0 regardless of active_fraction
            # (they are the READOUT, not a driver -- there is no biological
            # basis for a nonzero "baseline pathology level" independent of
            # the network's own upstream dynamics).
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


def run_sim(model, y0):
    return odeint(model.ode_system, y0, T, args=(None,), mxstep=10000)


def pct_reduction(traj_drug_col, traj_nodrug_col):
    a0 = _trapz(traj_nodrug_col)
    a1 = _trapz(traj_drug_col)
    return (a0 - a1) / a0 * 100 if a0 > 0 else np.nan


def output_active_index(model, pw):
    return model.STATE_VARS.index(pw + '_active')
