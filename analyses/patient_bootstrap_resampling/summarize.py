#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Summarize bootstrap_results.jsonl into median / 95% CI / %(ratio>1) for the
3 core pathway outputs. Run after run_bootstrap.py has finished (or partially
finished -- it will summarize whatever is in the file so far).

Usage:
    python summarize.py
"""
import json, os
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
IN_JSONL = os.path.join(HERE, 'bootstrap_results.jsonl')

CORE_OUTPUTS = ['P_Cell_death', 'P_Hepatocyte_injury', 'P_Inflammation']

# original manuscript point estimates (Table 4, GSE126848, ki=0.3)
ORIG = {'P_Cell_death': 2.44, 'P_Hepatocyte_injury': 2.10, 'P_Inflammation': 3.03}


def main():
    records = []
    with open(IN_JSONL) as f:
        for line in f:
            rec = json.loads(line)
            if 'error' not in rec:
                records.append(rec)

    print(f'Valid records: {len(records)}')
    print()

    for out in CORE_OUTPUTS:
        ratios = np.array([r[f'{out}_ratio'] for r in records
                            if r.get(f'{out}_ratio') is not None])
        if len(ratios) == 0:
            print(f'{out}: no valid ratios')
            continue
        median = np.median(ratios)
        ci_lo, ci_hi = np.percentile(ratios, [2.5, 97.5])
        pct_gt1 = np.mean(ratios > 1) * 100
        mean, sd = np.mean(ratios), np.std(ratios)
        print(f'{out}  (n={len(ratios)}):')
        print(f'  median NAFL/NASH ratio = {median:.3f}')
        print(f'  95% CI (percentile)    = [{ci_lo:.3f}, {ci_hi:.3f}]')
        print(f'  mean +/- SD            = {mean:.3f} +/- {sd:.3f}')
        print(f'  % iterations ratio>1   = {pct_gt1:.1f}%')
        print(f'  original point estimate= {ORIG[out]}')
        print()


if __name__ == '__main__':
    main()
