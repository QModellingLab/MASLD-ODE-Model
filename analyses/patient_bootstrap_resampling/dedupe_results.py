#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Deduplicate bootstrap_results.jsonl by 'iter' (keep first occurrence of each
iter index), report how many duplicates were removed, and overwrite the file
with the deduplicated version (a .bak backup of the original is kept).

Usage:
    python dedupe_results.py
"""
import json, os, shutil

HERE = os.path.dirname(os.path.abspath(__file__))
IN_JSONL = os.path.join(HERE, 'bootstrap_results.jsonl')
BAK = IN_JSONL + '.bak'

seen = {}
total_lines = 0
with open(IN_JSONL) as f:
    for line in f:
        line = line.strip()
        if not line:
            continue
        total_lines += 1
        rec = json.loads(line)
        it = rec['iter']
        if it not in seen:
            seen[it] = line

print(f'Total lines read      : {total_lines}')
print(f'Unique iter indices   : {len(seen)}')
print(f'Duplicate lines removed: {total_lines - len(seen)}')

missing = sorted(set(range(200)) - set(seen.keys()))
if missing:
    print(f'WARNING: missing iter indices (not yet run): {missing}')
else:
    print('All 200 iter indices (0-199) present.')

shutil.copy(IN_JSONL, BAK)
print(f'Backed up original to: {BAK}')

with open(IN_JSONL, 'w') as f:
    for it in sorted(seen.keys()):
        f.write(seen[it] + '\n')
print(f'Wrote deduplicated file back to: {IN_JSONL}')
