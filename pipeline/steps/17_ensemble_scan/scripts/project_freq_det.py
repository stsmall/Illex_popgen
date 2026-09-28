#!/usr/bin/env python3
"""Deterministic projection of a SweepFinder2 FreqFile to constant n0.
Unlike the hypergeometric (random) projection, this rounds each site's derived
count to the nearest integer at resolution 1/n0: xp = round(x*n0/n). This is the
MLE derived-allele frequency at the target sample size, adds NO per-site sampling
noise, and is fully reproducible (no seed) -> a stable, repeatable CLR profile.
Sites with n<n0 are dropped; sites that round to x=0 (mono-ancestral) are dropped;
x=n0 (fixed-derived substitutions) are kept.
Usage: project_freq_det.py <in.freq> <out.freq> <n0>
"""
import sys

inp, out, n0 = sys.argv[1], sys.argv[2], int(sys.argv[3])
kept = dropped_lown = dropped_mono = 0
with open(inp) as fh, open(out, "w") as o:
    o.write("position\tx\tn\tfolded\n")
    fh.readline()
    for line in fh:
        p = line.split()
        pos, x, n = p[0], int(p[1]), int(p[2])
        if n < n0:
            dropped_lown += 1
            continue
        xp = int(round(x * n0 / n))
        if xp <= 0:
            dropped_mono += 1
            continue
        if xp > n0:
            xp = n0
        o.write(f"{pos}\t{xp}\t{n0}\t0\n")
        kept += 1
sys.stderr.write(f"PROJdet kept={kept} dropped_lown={dropped_lown} dropped_mono={dropped_mono} n0={n0}\n")
