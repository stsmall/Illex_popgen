"""Throwaway ms.gz maker for the Task-6 dual-fvec toolchain proof.

Uses the committed harness modules (neutral_msprime, ms_export, ascertain) to emit
small neutral ``ms.gz`` files in the exact generator format both fvec toolchains
(diploSHIC + partialSHIC) must ingest. Not part of the training generator -- it only
exercises the toolchains on the generator's ms format.

For each of ``--nseeds`` seeds: simulate one 2.2 Mb neutral locus (n_dip=350 ->
700 haplotypes), slice the focal 1.1 Mb window (subwindow 5 centred, offset 550 kb),
and emit ONE ms replicate into each of two files (the sim is the expensive step, so
both come from the same sims):

  <out>              ascertained (step-13 depth/dropout/missing + MAF>0.01); carries
                     '?' missing calls -- this is the exact generator format.
  <out>.clean.msOut.gz   clean 0/1, MAF>0.01 filtered, NO missing '?'.

Why two files: the fvec toolchains do NOT treat '?' as missing -- diploSHIC maps it
to a spurious allele (ord('?')-ord('0')=15) and partialSHIC's allel.HaplotypeArray
crashes on it. So the dual-fvec proof runs on the clean file; the ascertained file
documents the '?' incompatibility (a Phase-B item: model missingness via mask files,
per diploSHIC convention, not inline '?').

Run with the slim_sim env python (has msprime):
  /home/ssmall/miniforge3/envs/slim_sim/bin/python make_test_ms.py OUT.msOut.gz [--nseeds 3]
"""
from __future__ import annotations

import argparse
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import neutral_msprime
import ms_export
import ascertain as ascertain_mod

MIN_MAF = 0.01


def maf_filter(geno, positions):
    """Keep polymorphic sites with minor-allele freq > MIN_MAF (clean 0/1, no missing)."""
    geno = np.asarray(geno, dtype=np.int8)
    n_hap = geno.shape[0]
    ac = geno.sum(0)
    af = ac / float(n_hap)
    maf = np.minimum(af, 1 - af)
    keep = (maf > MIN_MAF) & (ac > 0) & (ac < n_hap)
    return geno[:, keep], np.asarray(positions)[keep]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("out")
    ap.add_argument("--nseeds", type=int, default=3)
    ap.add_argument("--win", type=int, default=1_100_000)
    ap.add_argument("--sub", type=int, default=100_000)
    ap.add_argument("--n-sub", type=int, default=11)
    ap.add_argument("--L", type=float, default=2_200_000)
    ap.add_argument("--n-dip", type=int, default=350)
    args = ap.parse_args()

    offsets = ms_export.window_offsets(win=args.win, sub=args.sub, n_sub=args.n_sub)
    focal_idx = args.n_sub // 2                       # 5 for n_sub=11
    focal_offset = dict(offsets)[focal_idx]

    model = ascertain_mod.load_model()
    n_hap = 2 * args.n_dip
    clean_out = args.out.replace(".msOut.gz", "") + ".clean.msOut.gz"

    clean_recs, asc_recs = [], []
    for seed in range(1, args.nseeds + 1):
        ts = neutral_msprime.simulate_neutral(seed=seed, L=int(args.L), n_dip=args.n_dip)
        G, pos = ms_export.genotype_and_positions(ts)
        geno, positions_rel = ms_export.window_from(G, pos, focal_offset, args.win)

        cgeno, cpos = maf_filter(geno, positions_rel)
        clean_recs.append((cgeno, cpos))

        rng = np.random.default_rng(seed)
        ageno, apos = ascertain_mod.ascertain(geno, positions_rel, rng, model)
        asc_recs.append((ageno, apos))

        sys.stderr.write(
            "seed %d: raw %d snps -> clean %d snps, ascertained %d snps\n"
            % (seed, geno.shape[1], cgeno.shape[1], ageno.shape[1]))

    ms_export.write_ms_gz(clean_recs, clean_out, n_hap=n_hap, seed=0)
    ms_export.write_ms_gz(asc_recs, args.out, n_hap=n_hap, seed=0)
    sys.stderr.write("wrote %d reps each -> %s (clean) and %s (ascertained)\n"
                     % (len(clean_recs), clean_out, args.out))


if __name__ == "__main__":
    main()
