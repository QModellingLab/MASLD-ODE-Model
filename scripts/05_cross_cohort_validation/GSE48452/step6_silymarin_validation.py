#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
step6_silymarin_validation.py  (GSE48452)  [Fig 6 equivalent]
================================================================
Independent-cohort validation of the silymarin intervention (Fig 6 analogue).

Feeds GSE48452-derived node initial ratios into the EXISTING v10 model, then
applies the SAME silymarin multi-target inhibition used in the main manuscript
(Fig 6) and checks whether the key finding is reproduced:

    silymarin delays pathway saturation MORE at the early (NAFL) stage than at
    the advanced (NASH) stage  ->  "early-stage intervention is more effective".

Silymarin is modeled exactly as in the manuscript: an 8-target inhibitor that
multiplies the activation terms of CASP3, CASP7, CASP8, CYP2E1, IL_8, TNFa,
NF_kB, TGF_b1 by a factor ki (1.0 no drug; 0.7 low; 0.5 medium; 0.3 high).
The model body is UNCHANGED; only the node initial values are replaced by the
independent-cohort ratios (same approach as step5).

Reads ratios from step5 output (GSE48452_node_initial_ratios.xlsx). Run step5
first.

Outputs:
  GSE48452_silymarin_dynamics.xlsx   (t-half per stage x dose, per output)
  Fig_silymarin_validation.png        (Fig 6 style: stages x outputs, dose lines)
================================================================
"""
import os
import importlib.util
import types
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.integrate import odeint

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
MODEL_DIR = SCRIPT_DIR
MODEL_TEMPLATE = os.path.join(MODEL_DIR, 'ode_model_pydeseq2_NASH_vs_Normal_v10_mean.py')
RATIOS_XLSX = os.path.join(SCRIPT_DIR, 'GSE48452_node_initial_ratios.xlsx')

OUT_DYN = os.path.join(SCRIPT_DIR, 'GSE48452_silymarin_dynamics.xlsx')
OUT_FIG = os.path.join(SCRIPT_DIR, 'Fig_silymarin_validation.png')

# Stages to evaluate for the intervention (early vs advanced)
INTERVENTION_STAGES = ['NAFL', 'NASH']     # sheets in the ratios xlsx
# Silymarin doses (ki multiplier on target activation terms)
DOSES = [
    {'ki': 1.0, 'label': 'Disease (no drug)', 'color': '#c1121f', 'ls': '-',  'lw': 1.8},
    {'ki': 0.7, 'label': 'Silymarin low',     'color': '#F39C12', 'ls': '--', 'lw': 1.3},
    {'ki': 0.5, 'label': 'Silymarin medium',  'color': '#16A085', 'ls': '-.', 'lw': 1.3},
    {'ki': 0.3, 'label': 'Silymarin high',    'color': '#2980B9', 'ls': '-',  'lw': 1.6},
]
NORMAL_REF = {'label': 'Normal (reference)', 'color': '#1f3a5f', 'ls': ':', 'lw': 1.3}

# Silymarin 8 targets: (name, inactive_idx, active_idx) -- indices for v10 model
SILYMARIN_TARGETS = [
    ('CASP3', 18, 19), ('CASP7', 20, 21), ('CASP8', 22, 23),
    ('CYP2E1', 26, 27), ('IL_8', 58, 59), ('TNFa', 104, 105),
    ('NF_kB', 80, 81), ('TGF_b1', 100, 101),
]

# Outputs to display (silymarin acts on cell-death / inflammation / injury axes)
DISPLAY_OUTPUTS = ['P_Hepatocyte_injury', 'P_Cell_death', 'P_Inflammation']

T_MAX, N_POINTS = 300.0, 3001
HALF_TARGET = 50.0


def load_model(path):
    spec = importlib.util.spec_from_file_location('base', path)
    mod = importlib.util.module_from_spec(spec); spec.loader.exec_module(mod)
    return mod


def build_silymarin_source(base_source, ki):
    """Multiply the activation block of each silymarin target by ki."""
    lines = base_source.split('\n')
    for _, i_idx, a_idx in SILYMARIN_TARGETS:
        for idx in (i_idx, a_idx):
            for li, line in enumerate(lines):
                if line.lstrip().startswith(f'dydt[{idx}]'):
                    eq = line.find('='); rhs = line[eq + 1:].strip()
                    if rhs.startswith('-('):
                        sign, bs = '-', 1
                    elif rhs.startswith('('):
                        sign, bs = '', 0
                    else:
                        continue
                    depth = 0; ep = None
                    for p, ch in enumerate(rhs[bs:], start=bs):
                        if ch == '(':
                            depth += 1
                        elif ch == ')':
                            depth -= 1
                            if depth == 0:
                                ep = p; break
                    if ep is None:
                        continue
                    act = rhs[bs:ep + 1]; rem = rhs[ep + 1:]
                    ind = line[: len(line) - len(line.lstrip())]
                    lines[li] = f"{ind}dydt[{idx}] = {sign}(({act}) * {ki:.4f}){rem}"
                    break
    return '\n'.join(lines)


def load_module_from_source(source, name):
    mod = types.ModuleType(name); mod.__file__ = name
    exec(compile(source, name, 'exec'), mod.__dict__)
    return mod


def make_y0(model, ratio_map):
    y0 = np.zeros(len(model.STATE_VARS))
    for i, sv in enumerate(model.STATE_VARS):
        if sv.endswith('_inactive'):
            node = sv[:-len('_inactive')]
            if node in ratio_map:
                y0[i] = ratio_map[node]
            elif node in model.PATHWAY_METADATA:
                y0[i] = model.PATHWAY_METADATA[node].get('initial_level', 100.0)
            else:
                y0[i] = 1.0
    return y0


def t_half(t, traj, target=HALF_TARGET):
    above = np.where(traj >= target)[0]
    return float(t[above[0]]) if len(above) > 0 else np.inf


def main():
    print('=' * 64)
    print('GSE48452 silymarin intervention validation (Fig 6 equivalent)')
    print('=' * 64)
    for f in (MODEL_TEMPLATE, RATIOS_XLSX):
        if not os.path.exists(f):
            raise FileNotFoundError(f'Required file not found: {f}\n'
                                    f'(run step5 first to produce the ratios xlsx)')

    model = load_model(MODEL_TEMPLATE)
    with open(MODEL_TEMPLATE, 'r', encoding='utf-8') as fh:
        base_source = fh.read()
    p_idx = {pw: model.STATE_VARS.index(pw + '_active') for pw in model.PATHWAY_METADATA}
    t = np.linspace(0, T_MAX, N_POINTS)

    # Load GSE48452 node ratios per stage (from step5)
    ratios = {}
    xls = pd.ExcelFile(RATIOS_XLSX)
    for stage in INTERVENTION_STAGES:
        sheet = f'{stage}_ratios'
        if sheet not in xls.sheet_names:
            raise ValueError(f'Sheet {sheet} not in {RATIOS_XLSX}')
        df = pd.read_excel(xls, sheet)
        ratios[stage] = dict(zip(df['node'], df['ratio']))
    print(f'Loaded stage ratios: {list(ratios.keys())}')

    # Normal reference (all gene ratios = 1.0), no drug
    ratios_normal = {n: 1.0 for n in model.GENE_METADATA}

    # Pre-build silymarin modules per dose (compile once)
    sily_mods = {}
    for d in DOSES:
        if d['ki'] == 1.0:
            sily_mods[d['ki']] = model       # no drug = base model
        else:
            src = build_silymarin_source(base_source, d['ki'])
            sily_mods[d['ki']] = load_module_from_source(src, f"sily_{d['ki']}")
    print('Silymarin dose modules compiled.')

    # Simulate: Normal ref + each stage x each dose
    print('\nSimulating...')
    traj = {}   # (stage, ki) -> y
    y0_norm = make_y0(model, ratios_normal)
    traj[('Normal', 1.0)] = odeint(model.ode_system, y0_norm, t,
                                   args=(None,), mxstep=10000)
    for stage in INTERVENTION_STAGES:
        y0 = make_y0(model, ratios[stage])
        for d in DOSES:
            mod = sily_mods[d['ki']]
            traj[(stage, d['ki'])] = odeint(mod.ode_system, y0, t,
                                            args=(None,), mxstep=10000)
        print(f'  {stage}: {len(DOSES)} doses done.')

    # t-half table per output, stage, dose; plus drug effect (delay vs no-drug)
    rows = []
    for pw in DISPLAY_OUTPUTS:
        i = p_idx[pw]
        for stage in INTERVENTION_STAGES:
            th_nodrug = t_half(t, traj[(stage, 1.0)][:, i])
            for d in DOSES:
                th = t_half(t, traj[(stage, d['ki'])][:, i])
                delay = (th - th_nodrug) if (np.isfinite(th) and np.isfinite(th_nodrug)) else np.nan
                rows.append({'P_output': pw, 'stage': stage, 'ki': d['ki'],
                             'dose': d['label'], 't_half': th, 'delay_vs_nodrug': delay})
    dyn = pd.DataFrame(rows)
    dyn.to_excel(OUT_DYN, index=False)
    print(f'\nDynamics table saved: {OUT_DYN}')

    # Objective summary: drug effect, NAFL vs NASH.
    # Saturating outputs: t-half delay at ki=0.3 vs no-drug (larger = stronger).
    # Non-saturating outputs: steady-state reduction at ki=0.3 vs no-drug.
    def steady(y, i):
        return float(np.mean(y[-int(N_POINTS * 0.1):, i]))

    print('\n' + '=' * 70)
    print('Silymarin effect at highest dose (ki=0.3) vs no-drug, NAFL vs NASH')
    print('  saturating  -> t-half delay (h, larger = stronger)')
    print('  non-satur.  -> steady-state reduction (a.u., larger = stronger)')
    print('=' * 70)
    print(f'{"P_output":<25}{"metric":>10}{"NAFL":>10}{"NASH":>10}{"NAFL>NASH":>11}')
    for pw in DISPLAY_OUTPUTS:
        i = p_idx[pw]
        base_norm = steady(traj[('Normal', 1.0)], i)
        saturating = base_norm > HALF_TARGET
        eff = {}
        for stage in INTERVENTION_STAGES:
            if saturating:
                th0 = t_half(t, traj[(stage, 1.0)][:, i])
                th1 = t_half(t, traj[(stage, 0.3)][:, i])
                eff[stage] = (th1 - th0) if (np.isfinite(th0) and np.isfinite(th1)) else np.nan
            else:
                s0 = steady(traj[(stage, 1.0)], i)
                s1 = steady(traj[(stage, 0.3)], i)
                eff[stage] = s0 - s1     # reduction
        metric = 't-half' if saturating else 'steady'
        dn, ds = eff.get('NAFL', np.nan), eff.get('NASH', np.nan)
        flag = ('YES' if (np.isfinite(dn) and np.isfinite(ds) and dn > ds) else '-')
        print(f'{pw:<25}{metric:>10}{dn:>10.2f}{ds:>10.2f}{flag:>11}')
    print('\n  (NAFL effect > NASH effect  ->  silymarin more effective early)')

    # Figure: Fig 6 style. Rows = stages (NASH top, NAFL bottom), cols = outputs.
    n_stage = len(INTERVENTION_STAGES)
    n_out = len(DISPLAY_OUTPUTS)
    fig, axes = plt.subplots(n_stage, n_out, figsize=(4.3 * n_out, 3.6 * n_stage),
                             constrained_layout=True, squeeze=False)
    stage_order_fig = ['NASH', 'NAFL']   # advanced top, early bottom (like Fig 6)
    stage_order_fig = [s for s in stage_order_fig if s in INTERVENTION_STAGES]
    for r, stage in enumerate(stage_order_fig):
        for c, pw in enumerate(DISPLAY_OUTPUTS):
            ax = axes[r][c]; i = p_idx[pw]
            # Normal reference
            ax.plot(t, traj[('Normal', 1.0)][:, i], label=NORMAL_REF['label'],
                    color=NORMAL_REF['color'], ls=NORMAL_REF['ls'], lw=NORMAL_REF['lw'])
            for d in DOSES:
                ax.plot(t, traj[(stage, d['ki'])][:, i], label=d['label'],
                        color=d['color'], ls=d['ls'], lw=d['lw'])
            ax.set_title(f'{stage}: {pw.replace("P_", "")}', fontsize=9)
            ax.set_xlabel('Time (h)', fontsize=8)
            ax.set_ylabel('P_active (a.u.)', fontsize=8)
            ax.grid(alpha=0.25)
            if r == 0 and c == 0:
                ax.legend(fontsize=6.5, loc='lower right')
    plt.suptitle('GSE48452-driven silymarin intervention validation',
                 fontsize=12, fontweight='bold')
    plt.savefig(OUT_FIG, dpi=300, bbox_inches='tight')
    plt.close()
    print(f'\nFigure saved: {OUT_FIG}')
    print('\nDone.')


if __name__ == '__main__':
    main()
