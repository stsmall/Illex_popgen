#!/usr/bin/env python3
"""Thin a SweepFinder2 FreqFile. SF2's -l scan is O(grid x SNPs); the CLR is
robust to uniform thinning, so reducing SNP density cuts wall time ~linearly
while preserving the SFS shape and the sweep signal.

Two modes:
  bp   <dist>  keep the first SNP in each <dist>-bp bin (uniform spatial density)
  step <k>     keep every k-th SNP (exact uniform subsample of the SFS)

Usage: thin_freq.py <in.freq> <out.freq> bp   <dist>
       thin_freq.py <in.freq> <out.freq> step <k>
"""
import sys

inp, out, mode, param = sys.argv[1], sys.argv[2], sys.argv[3], int(sys.argv[4])
kept = total = 0
with open(inp) as fh, open(out, "w") as o:
    o.write(fh.readline())               # copy header
    if mode == "bp":
        next_ok = -1
        for line in fh:
            total += 1
            pos = int(line[:line.find("\t")])
            if pos >= next_ok:
                o.write(line)
                kept += 1
                next_ok = pos + param    # advance the window past this SNP
    elif mode == "step":
        for i, line in enumerate(fh):
            total += 1
            if i % param == 0:
                o.write(line)
                kept += 1
    else:
        sys.exit(f"bad mode {mode}")
sys.stderr.write(f"THIN mode={mode} param={param} kept={kept}/{total}\n")
