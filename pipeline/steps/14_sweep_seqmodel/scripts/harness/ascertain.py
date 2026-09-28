"""Ascertainment-matching for sim windows.

Applies the step-13 empirical low-depth model (per-genotype depth from the observed
depth histogram -> depth-scaled het->hom miscall on CALLED genotypes -> MAF/polymorphism
filter) so simulated windows carry the SAME calling error as the real 633-sample callset.
Emits COMPLETE 0/1 -- MISSINGNESS is NOT baked here: diploSHIC applies the empirical
per-genotype missing pattern (miss_rate ~0.294) to the complete sim via
``--vcfForMaskFileName`` (genoMask), and the empirical scan carries the real missingness
natively, so sims and observed data are masked consistently. (The step-13 ``ascertain_ms``
baked missing as ``?``, which BOTH fvec toolchains mishandle -> corrupt; removed here.)

Diploid folding: consecutive haplotype pairs (2i, 2i+1) are one individual, as
produced by the n=700 subsample of individuals in the SLiM reconstruction and by
``neutral_msprime`` (individual-consecutive samples).
"""
from __future__ import annotations

import json
import os

import numpy as np

# Portable: Talapas/HPC sets ASCERTAINMENT_JSON to the copied model; box default kept.
DEFAULT_MODEL = os.environ.get(
    "ASCERTAINMENT_JSON",
    "/sietch_colab/data_share/illex/popgen_data/analysis/"
    "steps/13_diploshic/ascertainment.json")


def load_model(path: str = DEFAULT_MODEL):
    """Load the step-13 ascertainment model into lookup arrays."""
    A = json.load(open(path))
    D = np.array([d for d, _ in A["depth_hist"]], dtype=int)
    P = np.array([p for _, p in A["depth_hist"]], dtype=float)
    P = P / P.sum()                                   # normalise (guards rounding)
    HD = dict((int(d), float(q)) for d, q in A["het_dropout_by_depth"])
    drop_arr = np.array([HD.get(d, 0.0) for d in range(int(D.max()) + 1)], dtype=float)
    min_maf = float(A.get("filters", {}).get("min_maf", 0.01))
    return {"D": D, "P": P, "drop_arr": drop_arr, "min_maf": min_maf}


def ascertain(geno, positions, rng, model=None):
    """Apply depth/dropout/missingness + MAF filter to one window.

    Parameters
    ----------
    geno : (n_hap, n_snp) int array, 0/1 haplotypes; n_hap even, pairs = diploids.
    positions : (n_snp,) array-like.
    rng : numpy Generator.
    model : dict from :func:`load_model` (loaded if None).

    Returns
    -------
    (geno2 [n_hap, n_snp2] int8, COMPLETE 0/1, positions2 [n_snp2]).
    Sites failing MAF>min_maf or monomorphic are dropped. Missingness is NOT baked in;
    diploSHIC applies the empirical missing pattern downstream via its genoMask.
    """
    if model is None:
        model = load_model()
    D, P, drop_arr, min_maf = model["D"], model["P"], model["drop_arr"], model["min_maf"]

    geno = np.asarray(geno, dtype=np.int64)
    n_hap, s = geno.shape
    positions = np.asarray(positions)
    if s == 0:
        return np.zeros((n_hap, 0), dtype=np.int8), positions[:0]

    ndip = n_hap // 2
    g = geno[:ndip * 2].reshape(ndip, 2, s).sum(1)            # (ndip, s) alt count 0/1/2

    depth = D[rng.choice(len(D), size=g.shape, p=P)]          # per-genotype depth
    het = g == 1
    # depth<2 genotypes are MISSING in the real callset, but we do NOT bake that here:
    # diploSHIC applies the empirical per-genotype missing pattern to the COMPLETE sim via
    # --vcfForMaskFileName (genoMask). So we model only the low-depth het->hom miscall ERROR on
    # CALLED genotypes (depth>=2) and emit complete 0/1; the depth<2 genotypes keep their true
    # value (they will be masked by diploSHIC's genoMask, matching the real missing pattern).
    called = depth >= 2
    drop = drop_arr[np.clip(depth.astype(int), 0, len(drop_arr) - 1)]
    het_to_hom = het & called & (rng.random(g.shape) < drop)      # low-depth het->hom on called genos
    g = np.where(het_to_hom, 2 * rng.integers(0, 2, g.shape), g)   # het -> 0 or 2

    # site filters: MAF > min_maf over non-missing calls, and polymorphic
    ac = np.where(g < 0, 0, g).sum(0)
    an = 2 * (g >= 0).sum(0)
    maf = np.minimum(ac, an - ac) / np.maximum(an, 1)
    keep = (maf > min_maf) & (ac > 0) & (ac < an)
    g = g[:, keep]
    positions2 = positions[keep]

    # re-expand diploid genotype -> two haplotypes; missing -> -1 on both
    a = np.where(g < 0, -1, np.where(g == 2, 1, 0)).astype(np.int8)   # hap A: 1 iff hom-alt
    b = np.where(g < 0, -1, np.where(g >= 1, 1, 0)).astype(np.int8)   # hap B: 1 iff carries alt
    geno2 = np.empty((ndip * 2, g.shape[1]), dtype=np.int8)
    geno2[0::2] = a
    geno2[1::2] = b
    return geno2, positions2
