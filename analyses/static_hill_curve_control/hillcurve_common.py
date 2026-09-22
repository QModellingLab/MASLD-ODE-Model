#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Shared utilities for the random_target_control analysis (R1-7 / R2-4 / R3-5).

Question being addressed:
  Is the predicted stage-dependent silymarin effect specific to the 8
  literature-based targets (CASP3, CASP7, CASP8, CYP2E1, IL_8/CXCL8, TNFa,
  NF_kB, TGF_b1), or would inhibiting ANY random set of 8 nodes at the same
  ki=0.3 produce a similar (or stronger) NAFL>NASH advantage -- simply
  because several targets are output-proximal / highly-connected nodes?

Two null models are built here:
  1. "random_target"    -- draw 8 nodes uniformly at random from all 58
                            molecular nodes (same pool the real silymarin
                            targets were drawn from). Addresses R1-7 / R2-4.
  2. "matched_control"   -- draw 8 nodes from the pool EXCLUDING both the
                            real silymarin targets AND the top sensitivity
                            drivers reported in Fig. 6 (CASP7, CASP3, Cytc,
                            Bax, TNFa, FasL, IL_8, TGF_b1). Addresses R3-5's
                            explicit request for a control that excludes
                            "targets or top drivers".

Topology, kinetics (Vmax=ksp=n=kcat=2.0), and initial conditions are the
REAL hsa04932 network and the REAL GSE126848-derived node ratios (baked
into the model files' GENE_METADATA) -- nothing about the network is
scrambled here (that is what edge_scramble_control does). Only the
IDENTITY of the 8 inhibited nodes changes between replicates.

Because the untreated (no-drug) baseline trajectory does not depend on
which nodes are chosen for inhibition, it is computed ONCE and reused
across all replicates -- each replicate only requires 2 new ODE solves
(NASH-drug, NAFL-drug), which is why this analysis is much faster than
the patient-bootstrap or edge-scramble analyses.
"""
import os, sys, json, types, importlib.util
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
N_TARGETS = 8

# The 8 literature-based silymarin targets (KEGG display names, matching
# GENE_METADATA / STATE_VARS naming in the model files).
SILYMARIN_TARGETS = ['CASP3', 'CASP7', 'CASP8', 'CYP2E1', 'IL_8', 'TNFa', 'NF_kB', 'TGF_b1']

# Top local-sensitivity drivers reported in Fig. 6 / Results & Discussion
# (dominant driver per core output; display names):
#   P_Hepatocyte_injury : CASP7 (0.137), CASP3 (0.072), Cytc (0.053), Bax (0.022)
#   P_Cell_death         : TNFa (0.161), FasL (0.110)
#   P_Inflammation       : IL_8 (0.021); TGF_b1 noted as minor contributor
TOP_DRIVERS = ['CASP7', 'CASP3', 'Cytc', 'Bax', 'TNFa', 'FasL', 'IL_8', 'TGF_b1']

EXCLUDE_FOR_MATCHED_CONTROL = sorted(set(SILYMARIN_TARGETS) | set(TOP_DRIVERS))


def _trapz(y):
    f = getattr(np, 'trapezoid', None) or np.trapz
    return float(f(y, T))


def load_model(path):
    spec = importlib.util.spec_from_file_location(
        os.path.splitext(os.path.basename(path))[0], path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def eligible_pool(model):
    """All 58 molecular node names (GENE_METADATA keys), i.e. the same
    universe the real silymarin targets and top drivers were drawn from."""
    return sorted(model.GENE_METADATA.keys())


def make_y0(model):
    """Y0 straight from the model's own baked-in GENE_METADATA ratios and
    PATHWAY_METADATA initial levels -- these are the fixed GSE126848 main-
    cohort values used for the manuscript's Table 4 point estimates. No
    resampling here (that's what patient_bootstrap_resampling does)."""
    y0 = np.zeros(len(model.STATE_VARS))
    for i, sv in enumerate(model.STATE_VARS):
        if sv.endswith('_inactive'):
            node = sv[:-len('_inactive')]
            if node in model.GENE_METADATA:
                y0[i] = model.GENE_METADATA[node]['initial_ratio']
            elif node in model.PATHWAY_METADATA:
                y0[i] = model.PATHWAY_METADATA[node]['initial_level']
            else:
                raise KeyError(f'Unrecognized node {node!r} for {sv!r}')
        # '_active' forms stay at 0 (already zero-initialized)
    return y0


def build_ki_source(base_src, model, target_nodes, ki):
    """Return a modified copy of the model's source with a uniform
    inhibitory multiplier ki applied to the activation-rate term of each
    target node's active/inactive ODEs (same mechanism as the manuscript's
    Eq. 4, generalized here to an arbitrary node list instead of the fixed
    8 silymarin targets)."""
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


def already_done(jsonl_path):
    done = set()
    if os.path.exists(jsonl_path):
        with open(jsonl_path) as fh:
            for line in fh:
                try:
                    done.add(json.loads(line)['replicate'])
                except Exception:
                    pass
    return done
