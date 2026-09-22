import numpy as np
import itertools
import pandas as pd

# Table 4 data: cohort -> output -> (NAFL%, NASH%)
data = {
    "GSE126848": {"P_Cell_death": (13.85, 5.67), "P_Hepatocyte_injury": (34.31, 16.30), "P_Inflammation": (9.01, 2.97)},
    "GSE48452":  {"P_Cell_death": (29.98, 22.44), "P_Hepatocyte_injury": (40.39, 46.63), "P_Inflammation": (74.27, 69.85)},
    "GSE89632":  {"P_Cell_death": (11.86, 6.08),  "P_Hepatocyte_injury": (49.05, 55.06), "P_Inflammation": (89.55, 90.04)},
    "GSE130970": {"P_Cell_death": (6.84, 4.12),   "P_Hepatocyte_injury": (93.06, 81.17), "P_Inflammation": (78.56, 7.11)},
    "GSE162694": {"P_Cell_death": (28.21, 21.39), "P_Hepatocyte_injury": (41.32, 40.80), "P_Inflammation": (64.09, 9.77)},
    "GSE213621": {"P_Cell_death": (42.61, 46.36), "P_Hepatocyte_injury": (62.18, 63.49), "P_Inflammation": (28.40, 7.40)},
}

cohorts = list(data.keys())
outputs = ["P_Cell_death", "P_Hepatocyte_injury", "P_Inflammation"]

# log-ratio matrix: rows=cohort, cols=output
logR = np.array([[np.log(data[c][o][0]/data[c][o][1]) for o in outputs] for c in cohorts])
print("log(NAFL/NASH ratio) matrix (rows=cohort, cols=output):")
print(pd.DataFrame(logR, index=cohorts, columns=outputs).round(3))

obs_mean = logR.mean()
obs_median = np.median(logR)
print(f"\nObserved mean(log ratio) = {obs_mean:.4f}  (geometric mean ratio = {np.exp(obs_mean):.3f})")
print(f"Observed median(log ratio) = {obs_median:.4f}")

# ---- Cohort-level sign-flip permutation test (exact, 2^6 = 64 patterns) ----
# Rationale: outputs within a cohort share the same transcriptome/model realization,
# so the exchangeable unit under H0 (no true NAFL>NASH direction) is the whole cohort,
# not each cohort-output cell. Flipping = swapping NAFL<->NASH labels for all 3 outputs at once.
n = len(cohorts)
null_means = []
null_medians = []
for signs in itertools.product([1, -1], repeat=n):
    flipped = logR * np.array(signs)[:, None]
    null_means.append(flipped.mean())
    null_medians.append(np.median(flipped))
null_means = np.array(null_means)
null_medians = np.array(null_medians)

p_mean = np.mean(null_means >= obs_mean)      # one-sided, matches pre-specified direction
p_median = np.mean(null_medians >= obs_median)

print(f"\nExact cohort-level sign-flip permutation test (64 patterns):")
print(f"  one-sided P (mean statistic)   = {p_mean:.4f}")
print(f"  one-sided P (median statistic) = {p_median:.4f}")


print("\n" + "="*60)
print("補充分析")
print("="*60)

# Cohort-level summary (one number per cohort = mean across 3 outputs)
cohort_means = logR.mean(axis=1)
print("\nCohort-level mean log(ratio) (collapsing 3 correlated outputs to 1/cohort):")
for c, v in zip(cohorts, cohort_means):
    print(f"  {c}: {v:.3f}  (ratio={np.exp(v):.2f})")

from scipy.stats import wilcoxon
# one-sided Wilcoxon signed-rank on the 6 cohort-level means vs 0
stat, p_wilcoxon = wilcoxon(cohort_means, alternative='greater')
print(f"\nWilcoxon signed-rank (n=6 cohorts, H1: median log-ratio > 0): P = {p_wilcoxon:.4f}")

# Mixed-effects model: log(ratio) ~ 1 + (1|cohort), output as repeated measure
import statsmodels.formula.api as smf
df = pd.DataFrame({
    "logratio": logR.flatten(),
    "cohort": np.repeat(cohorts, len(outputs)),
    "output": outputs * len(cohorts),
})
md = smf.mixedlm("logratio ~ 1", df, groups=df["cohort"])
mdf = md.fit(reml=True)
print("\nMixed-effects model (random intercept per cohort), REML:")
print(mdf.summary().tables[1])

