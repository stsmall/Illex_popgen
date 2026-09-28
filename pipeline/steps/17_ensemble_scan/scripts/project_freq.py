#!/usr/bin/env python3
"""Project a SweepFinder2 FreqFile to a constant sample size n0 via hypergeometric
subsampling (unbiased SFS). SF2's spectrum step is ~O(n^3) and sums over every
distinct sample size, so variable n (missing data) makes it crawl; a constant n
makes it fast. Sites with n<n0 are dropped; sites that become fixed-ancestral
(x=0) after projection are dropped; fixed-derived (x=n0) are kept (substitutions
improve CLR power).
Usage: project_freq.py <in.freq> <out.freq> <n0> <seed>
"""
import sys
import numpy as np

inp, out, n0, seed = sys.argv[1], sys.argv[2], int(sys.argv[3]), int(sys.argv[4])
rng = np.random.default_rng(seed)
kept = dropped_lown = dropped_mono = 0
with open(inp) as fh, open(out, "w") as o:
    o.write("position\tx\tn\tfolded\n")
    header = fh.readline()
    for line in fh:
        p = line.split()
        pos, x, n = p[0], int(p[1]), int(p[2])
        if n < n0:
            dropped_lown += 1
            continue
        xp = int(rng.hypergeometric(x, n - x, n0))   # derived count in n0 drawn from n
        if xp == 0:
            dropped_mono += 1
            continue
        o.write(f"{pos}\t{xp}\t{n0}\t0\n")
        kept += 1
sys.stderr.write(f"PROJ kept={kept} dropped_lown={dropped_lown} dropped_mono={dropped_mono} n0={n0}\n")
