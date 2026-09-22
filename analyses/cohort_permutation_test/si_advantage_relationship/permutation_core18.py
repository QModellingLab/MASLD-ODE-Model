import numpy as np
import pandas as pd
from scipy import stats

np.random.seed(42)
N_PERM = 20000

df = pd.read_excel('Table_AllOutputs_SeparationIndex_vs_AUCratio.xlsx', sheet_name='AllOutputs_72combos')
core = df[df['P_output'].isin(['P_Cell_death','P_Hepatocyte_injury','P_Inflammation'])].copy()
valid = core.dropna(subset=['Separation_Index_%', 'NAFL/NASH_ratio']).copy()
print(f"Core N (3 outputs x 6 cohorts) = {len(valid)}")

x = valid['Separation_Index_%'].values
y = valid['NAFL/NASH_ratio'].values
cohort = valid['Dataset'].values
idx_by_cohort = {c: np.where(cohort == c)[0] for c in np.unique(cohort)}

obs_r, obs_p_asymp = stats.spearmanr(x, y)
print(f"Observed Spearman r = {obs_r:.4f} (asymptotic p = {obs_p_asymp:.4f})  [manuscript reports r=0.47, p=0.048]")

# global permutation
null_g = np.array([stats.spearmanr(x, np.random.permutation(y))[0] for _ in range(N_PERM)])
print(f"  global permutation one-sided P = {np.mean(null_g >= obs_r):.5f}")

# cohort-block permutation (within-cohort shuffle of the 3 outputs' ratios)
null_b = np.zeros(N_PERM)
for i in range(N_PERM):
    y_perm = y.copy()
    for c, idxs in idx_by_cohort.items():
        y_perm[idxs] = np.random.permutation(y_perm[idxs])
    null_b[i] = stats.spearmanr(x, y_perm)[0]
print(f"  cohort-block permutation one-sided P = {np.mean(null_b >= obs_r):.5f}")

# Mann-Whitney SI>10%
valid['separated'] = valid['Separation_Index_%'] > 10.0
g1 = valid.loc[valid['separated'], 'NAFL/NASH_ratio'].values
g2 = valid.loc[~valid['separated'], 'NAFL/NASH_ratio'].values
mw_stat, mw_p = stats.mannwhitneyu(g1, g2, alternative='greater')
print(f"\nObserved Mann-Whitney (asymptotic) one-sided P = {mw_p:.5f}  [manuscript reports P=0.0059]")
print(f"  N separated={len(g1)}, N not={len(g2)}")

sep_mask = valid['separated'].values
obs_stat = np.median(y[sep_mask]) - np.median(y[~sep_mask])
null_mw = np.zeros(N_PERM)
for i in range(N_PERM):
    y_perm = y.copy()
    for c, idxs in idx_by_cohort.items():
        y_perm[idxs] = np.random.permutation(y_perm[idxs])
    null_mw[i] = np.median(y_perm[sep_mask]) - np.median(y_perm[~sep_mask])
p_mw_block = np.mean(null_mw >= obs_stat)
print(f"  cohort-block permutation (median diff) one-sided P = {p_mw_block:.5f}")

