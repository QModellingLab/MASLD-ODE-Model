import numpy as np
import pandas as pd
from scipy import stats

np.random.seed(42)
N_PERM = 20000

df = pd.read_excel('Table_AllOutputs_SeparationIndex_vs_AUCratio.xlsx', sheet_name='AllOutputs_72combos')
main = df[~df['P_output'].isin(['P_Improvement_of_NAFLD'])].copy()
valid = main.dropna(subset=['Separation_Index_%', 'NAFL/NASH_ratio']).copy()
print(f"Valid N (11 disease-positive outputs x 6 cohorts, non-missing) = {len(valid)}")

x = valid['Separation_Index_%'].values
y = valid['NAFL/NASH_ratio'].values
cohort = valid['Dataset'].values

obs_r, obs_p_asymp = stats.spearmanr(x, y)
print(f"\nObserved Spearman r = {obs_r:.4f} (asymptotic p = {obs_p_asymp:.4f})")

# --- Test 1: global label permutation (standard nonparametric permutation test for Spearman r) ---
null_global = np.zeros(N_PERM)
for i in range(N_PERM):
    y_perm = np.random.permutation(y)
    null_global[i] = stats.spearmanr(x, y_perm)[0]
p_global = np.mean(null_global >= obs_r)
print(f"\n[Test 1] Global permutation (unconstrained shuffle of ratio labels, N={N_PERM}):")
print(f"  one-sided permutation P = {p_global:.5f}")

# --- Test 2: cohort-block permutation (shuffle ratio only WITHIN each cohort) ---
# Addresses non-independence: outputs sharing a cohort share the same simulated
# trajectory set / transcriptome, so exchangeability should be assessed within cohort.
null_block = np.zeros(N_PERM)
idx_by_cohort = {c: np.where(cohort == c)[0] for c in np.unique(cohort)}
for i in range(N_PERM):
    y_perm = y.copy()
    for c, idxs in idx_by_cohort.items():
        y_perm[idxs] = np.random.permutation(y_perm[idxs])
    null_block[i] = stats.spearmanr(x, y_perm)[0]
p_block = np.mean(null_block >= obs_r)
print(f"\n[Test 2] Cohort-block permutation (shuffle within cohort only, N={N_PERM}):")
print(f"  one-sided permutation P = {p_block:.5f}")

# --- Test 3: Mann-Whitney permutation-based p-value (SI>10% vs <=10%), same block structure ---
valid['separated'] = valid['Separation_Index_%'] > 10.0
obs_mw_stat, obs_mw_p = stats.mannwhitneyu(
    valid.loc[valid['separated'], 'NAFL/NASH_ratio'],
    valid.loc[~valid['separated'], 'NAFL/NASH_ratio'], alternative='greater')
print(f"\n[Reference] Mann-Whitney U (parametric asymptotic), one-sided P = {obs_mw_p:.5f}")

null_mw = np.zeros(N_PERM)
sep_mask = valid['separated'].values
for i in range(N_PERM):
    y_perm = y.copy()
    for c, idxs in idx_by_cohort.items():
        y_perm[idxs] = np.random.permutation(y_perm[idxs])
    g1 = y_perm[sep_mask]
    g2 = y_perm[~sep_mask]
    null_mw[i] = np.median(g1) - np.median(g2)
obs_stat = np.median(y[sep_mask]) - np.median(y[~sep_mask])
p_mw_perm = np.mean(null_mw >= obs_stat)
print(f"\n[Test 3] Cohort-block permutation of median-difference (separated vs not), N={N_PERM}:")
print(f"  observed median diff = {obs_stat:.3f}")
print(f"  one-sided permutation P = {p_mw_perm:.5f}")

