#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Shared utilities for joint_n_ksp_grid (R3 minor-3)."""
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
        os.path.splitext(os.path.basename(path))[0] + '_nksp', path)
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


def parameterize_nksp(src, n, ksp):
    """Rewrite all Hill-kinetics constants in the model source to use Hill
    coefficient n and half-saturation constant ksp, keeping kcat/Vmax fixed
    at 2.0 (matching the manuscript's own Supplementary Table S4 approach,
    which likewise varies only n and ksp).

    Every Hill term in the source has the form
        kcat * Regulator * Substrate**2.0 / (2.0**2.0 + Substrate**2.0)
    i.e. TWO occurrences of the literal exponent 2.0 (n) plus one base 2.0
    (ksp, always immediately followed by **2.0 in the denominator) plus one
    unrelated multiplicative 2.0 (kcat), which must be left untouched.

    Verified (see joint_n_ksp_grid development notes / README) to exactly
    reproduce the manuscript's already-published Supplementary Table S4
    values for (n=1,ksp=2) and (n=2,ksp=1), and Table 4's values for the
    primary (n=2,ksp=2) case.
    """
    # Step 1: ksp**n terms (denominator half-saturation) -- must run BEFORE
    # step 2, since this pattern's base "2.0" would otherwise be mistaken
    # for a bare exponent by step 2's blanket "**2.0"->"**n" replacement.
    src = src.replace('2.0**2.0', f'{ksp}**{n}')
    # Step 2: all remaining bare substrate/product exponents (the Hill
    # coefficient n applied directly to a state variable, e.g. "Bax_inactive**2.0").
    src = src.replace('**2.0', f'**{n}')
    return src


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
