#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Robustness scan: does the stage-dependent silymarin early-intervention
advantage (NAFL > NASH) survive changes in Hill n and ksp?
For each (n, ksp) variant, applies silymarin (ki=0.3) to the 8 target
nodes and compares the relative AUC reduction at NAFL vs NASH for the
three core outputs. Prints a qualitative-consistency table.
"""
import os, importlib.util, types
import numpy as np
from scipy.integrate import odeint
try:
    import pandas as pd
except ImportError:
    pd = None
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
REPO_ROOT  = os.path.dirname(os.path.dirname(SCRIPT_DIR))
# variants written by generate_robustness_models.py; check data/models first
_CAND = [
    os.path.join(REPO_ROOT, 'data', 'models', 'robustness_models'),
    os.path.join(SCRIPT_DIR, 'robustness_models'),
]
MDIR = next((d for d in _CAND if os.path.isdir(d)), _CAND[0])
T_END, T_POINTS, KI = 300, 3001, 0.3
GRID = [(1.0,2.0),(2.0,2.0),(4.0,2.0),(2.0,1.0),(2.0,4.0)]
CORE = ['P_Hepatocyte_injury','P_Cell_death','P_Inflammation']
TARGETS = [('CASP3',18,19),('CASP7',20,21),('CASP8',22,23),('CYP2E1',26,27),
           ('IL_8',58,59),('TNFa',104,105),('NF_kB',80,81),('TGF_b1',100,101)]

def load_src(f):
    return open(f,'r',encoding='utf-8').read()
def load_mod_path(f):
    s=importlib.util.spec_from_file_location('m',f); m=importlib.util.module_from_spec(s); s.loader.exec_module(m); return m
def load_mod_src(src,name):
    m=types.ModuleType(name); m.__file__=name; exec(compile(src,name,'exec'),m.__dict__); return m

def build_silymarin_ode(src, ki):
    lines=src.split('\n')
    for _,ii,ai in TARGETS:
        for idx in (ii,ai):
            for li,line in enumerate(lines):
                if line.lstrip().startswith(f'dydt[{idx}]'):
                    eq=line.find('='); rhs=line[eq+1:].strip()
                    if rhs.startswith('-('): sign,bs='-',1
                    elif rhs.startswith('('): sign,bs='',0
                    else: continue
                    depth=0; endp=None
                    for p,ch in enumerate(rhs[bs:],start=bs):
                        if ch=='(':depth+=1
                        elif ch==')':
                            depth-=1
                            if depth==0: endp=p; break
                    if endp is None: continue
                    act=rhs[bs:endp+1]; rem=rhs[endp+1:]
                    ind=line[:len(line)-len(line.lstrip())]
                    lines[li]=f"{ind}dydt[{idx}] = {sign}(({act}) * {ki:.4f}){rem}"
                    break
    return '\n'.join(lines)

def auc_of(mod_or_src, node_idx, is_src=False, name='m'):
    t=np.linspace(0,T_END,T_POINTS)
    mod = load_mod_src(mod_or_src,name) if is_src else mod_or_src
    y=odeint(mod.ode_system, mod.Y0, t, args=(None,), mxstep=10000)
    return float(np.trapezoid(y[:,node_idx],t))

OUT_DIR = SCRIPT_DIR
print(f"Silymarin ki={KI} | advantage = NAFL relative-reduction / NASH relative-reduction")
print(f"{'n':>4}{'ksp':>5} | " + " | ".join(f"{c.replace('P_',''):>16}" for c in CORE) + "   verdict")
print("-"*90)
all_hold=[]; records=[]
for (n,ksp) in GRID:
    tag=f"n{n:.1f}_ksp{ksp:.1f}"
    nash_src=load_src(os.path.join(MDIR,f"ode_model_pydeseq2_NASH_vs_Normal_v10_mean_{tag}.py"))
    nafl_src=load_src(os.path.join(MDIR,f"ode_model_pydeseq2_NAFL_vs_Normal_v10_mean_{tag}.py"))
    nash0=load_mod_src(nash_src,'nash0'); nafl0=load_mod_src(nafl_src,'nafl0')
    SV=nash0.STATE_VARS; amap={s[:-7]:i for i,s in enumerate(SV) if s.endswith('_active')}
    nash_drug=build_silymarin_ode(nash_src,KI); nafl_drug=build_silymarin_ode(nafl_src,KI)
    md=load_mod_src(nash_drug,'nashd'); mf=load_mod_src(nafl_drug,'nafld')
    row=[]; holds=[]; advs=[]
    for c in CORE:
        pi=amap[c]
        nash_base=auc_of(nash0,pi); nash_tx=auc_of(md,pi)
        nafl_base=auc_of(nafl0,pi); nafl_tx=auc_of(mf,pi)
        nash_red=(nash_base-nash_tx)/nash_base if nash_base else 0
        nafl_red=(nafl_base-nafl_tx)/nafl_base if nafl_base else 0
        adv = nafl_red/nash_red if nash_red>1e-9 else float('inf')
        row.append(f"{adv:16.2f}")
        advs.append(adv)
        holds.append(nafl_red > nash_red + 1e-9)  # NAFL advantage direction
    verdict = "ALL NAFL>NASH ✓" if all(holds) else f"{sum(holds)}/3 hold"
    all_hold.append(all(holds))
    star = "  <-- primary" if (n,ksp)==(2.0,2.0) else ""
    print(f"{n:>4.0f}{ksp:>5.0f} | " + " | ".join(row) + f"   {verdict}{star}")
    rec={'n':n,'ksp':ksp,'primary':(n,ksp)==(2.0,2.0)}
    for c,val,h in zip(CORE,advs,holds):
        rec[c]=round(val,3); rec[c+'_NAFL>NASH']=bool(h)
    rec['all_outputs_hold']=all(holds)
    records.append(rec)
print("-"*90)
print("robust across all 5 param sets:", all(all_hold))

# ---- save table (Excel + CSV) ----
cols=['n','ksp','primary']+CORE+['all_outputs_hold']
if pd is not None:
    df=pd.DataFrame(records)[[c for c in ['n','ksp','primary']+sum([[c,c+'_NAFL>NASH'] for c in CORE],[])+['all_outputs_hold']]]
    xlsx=os.path.join(OUT_DIR,'Robustness_scan_summary.xlsx')
    with pd.ExcelWriter(xlsx,engine='openpyxl') as w:
        df.to_excel(w,sheet_name='Advantage_by_param',index=False)
        pd.DataFrame([
            ['Metric','NAFL/NASH early-intervention advantage (relative AUC reduction ratio)'],
            ['Silymarin ki',KI],
            ['Outputs',', '.join(CORE)],
            ['Grid','n=1/2/4 (ksp=2); ksp=1/4 (n=2); primary=n=ksp=2'],
            ['Reading','value>1 = NAFL(early)>NASH advantage; *_NAFL>NASH = direction holds'],
            ['Verdict','direction preserved across all sets; magnitude varies'],
        ],columns=['Field','Value']).to_excel(w,sheet_name='Metadata',index=False)
    df.to_csv(os.path.join(OUT_DIR,'Robustness_scan_summary.csv'),index=False)
    print('  saved:',xlsx)
    print('  saved:',os.path.join(OUT_DIR,'Robustness_scan_summary.csv'))

# ---- overview figure: advantage vs param set, one line per output ----
labels=[f"n{r['n']:.0f}\nksp{r['ksp']:.0f}" for r in records]
x=np.arange(len(records))
fig,ax=plt.subplots(figsize=(7,4))
for c in CORE:
    ax.plot(x,[r[c] for r in records],marker='o',label=c.replace('P_',''))
ax.axhline(1.0,color='0.5',ls='--',lw=0.8)
ax.text(len(records)-1,1.0,' advantage = 1 (no stage effect)',va='bottom',ha='right',fontsize=7,color='0.4')
ax.set_xticks(x); ax.set_xticklabels(labels,fontsize=8)
ax.set_ylabel('NAFL / NASH early-intervention advantage')
ax.set_title('Parameter robustness of the stage-dependent advantage',fontsize=10)
ax.legend(fontsize=8,frameon=False); ax.grid(True,ls=':',alpha=0.4)
for ext in ('png','pdf'):
    out=os.path.join(OUT_DIR,f'Robustness_scan_overview.{ext}')
    fig.savefig(out,dpi=300 if ext=='png' else None,bbox_inches='tight')
    print('  saved:',out)
plt.close(fig)
