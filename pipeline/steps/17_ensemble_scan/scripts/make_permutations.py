#!/usr/bin/env python
"""Write division-label files for the permutation null: perm_00 = real labels, perm_01..perm_N =
labels shuffled among the 632 individuals (division sizes preserved). Usage: make_permutations.py OUTDIR N"""
import sys, os, numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); out, n = sys.argv[1], int(sys.argv[2])
rows = [l.split() for l in open(f"{HERE}/../config/divisions_632.tsv") if not l.startswith("#")]
iids = [r[0] for r in rows]; labs = np.array([r[1] for r in rows])
os.makedirs(out, exist_ok=True)
for k in range(n + 1):
    lab = labs if k == 0 else np.random.default_rng(20251019 + k).permutation(labs)
    with open(f"{out}/perm_{k:02d}.tsv", "w") as fh:
        fh.write("#IID\tPOP\n"); fh.writelines(f"{i}\t{l}\n" for i, l in zip(iids, lab))
print(f"wrote {n + 1} label files to {out}")
