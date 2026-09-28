#!/usr/bin/env python
# ============================================================================
# generate.py -- single-sim end-to-end generator for the Illex diploSHIC /
#                partialSHIC SLiM 9-class training-data harness (Task 5).
#
# One invocation == one simulation of a 2.2 Mb locus with the sweep at the
# centre (XPOS = 1.1 Mb). It composes the already-committed harness modules:
#
#     neutral        -> neutral_msprime.simulate_neutral()          (pure msprime)
#     hard / soft     \  slim.run_slim_sweep.run_slim() (SLiM recipe)
#     partialHard /   -> + slim.run_slim_sweep.reconstruct()
#     partialSoft    /    (subsample n=700 -> recapitate -> overlay; partial
#                          mode strips the segregating beneficial)
#
# then, from the reconstructed n=700 tree sequence, it slices the 11 overlapping
# 1.1 Mb analysis windows (ms_export.window_offsets: the sweep lands in each of
# the 11 subwindows in turn), ascertainment-matches each window to the real
# step-13 callset (ascertain.ascertain), and writes gzipped-ms training windows
# routed by class + a per-window manifest row:
#
#     neutral      : all 11 windows            -> neutral.msOut.gz
#     sweep klass C: focal window (subwin 5)   -> <C>.msOut.gz
#                    the other 10 (subwin != 5) -> linked<C>.msOut.gz
#
# The genotype matrix is computed ONCE and windowed 11 times. The beneficial
# site (+-1 bp of XPOS) is excluded before windowing (in a completed sweep it is
# fixed->monomorphic and MAF-dropped anyway, but conditioning is freq>=0.95 so in
# a subsample it can segregate at high freq; in a partial sweep it is explicitly
# segregating -- either way it is NOT neutral variation and must not leak in).
#
# DETERMINISM (all derived from --seed):
#   ss = numpy.SeedSequence(seed); s_child, asc_child = ss.spawn(2)
#     * s_child   -> the log-uniform s-draw stream  (used only for a sweep when
#                    --s is not supplied)
#     * asc_child -> the ascertainment RNG stream    (used for every class)
#   The SLiM master seed and the reconstruction subsample RNG both use the raw
#   --seed (as in run_slim_sweep), so all four streams are independent and the
#   whole pipeline is reproducible from --seed alone.
#
# Thread-capped (OMP/BLAS = 1) in-process BEFORE numpy import so a campaign can
# fan many of these out under the shared-machine core budget.
#
# CLI:
#   generate.py --klass {neutral,hard,soft,partialHard,partialSoft} --seed S
#               [--s SEL] [--F0 F] --outdir DIR [--L 2200000] [--n-dip 350]
# ============================================================================
from __future__ import annotations

import os

# --- thread-cap BEFORE numpy/BLAS import (shared machine; campaigns fan out) --
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
           "NUMEXPR_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    os.environ.setdefault(_v, "1")

import argparse            # noqa: E402
import csv                 # noqa: E402
import sys                 # noqa: E402

import numpy as np         # noqa: E402

# --- import the harness modules by ABSOLUTE path (CWD is not the harness dir) -
_HERE = os.path.dirname(os.path.abspath(__file__))
_SLIM_DIR = os.path.join(_HERE, "slim")
for _p in (_HERE, _SLIM_DIR):
    if _p not in sys.path:
        sys.path.insert(0, _p)

from ms_export import (window_offsets, genotype_and_positions,  # noqa: E402
                       window_from, write_ms_gz)
from neutral_msprime import simulate_neutral                    # noqa: E402
from ascertain import load_model, ascertain                     # noqa: E402
import run_slim_sweep as rss                                     # noqa: E402

# ---------------------------------------------------------------------------
# constants (PLAN_03 Global Constraints)
# ---------------------------------------------------------------------------
SLIM = os.environ.get("SLIM_BIN", "/home/ssmall/miniforge3/envs/slim_sim/bin/slim")
WIN = 1_100_000            # analysis window (bp)
SUB = 100_000             # subwindow (bp)
N_SUB = 11                # number of subwindows / analysis windows
FOCAL_SUBWIN = 5          # central subwindow == the focal class
S_LO, S_HI = 1e-4, 1e-2   # selection prior: s ~ log-U(1e-4, 1e-2)

# demography / rescaling (verbatim from run_slim_sweep defaults / the recipes)
Q = 100.0
N0 = 6_808_096.0
NREF = 547_928.0
TGROW = 769_519.0
MU = 3e-9
R = 2.1e-9
TAU_GEN = 20_000.0
FIX_THRESH = 0.95
F_LO, F_HI, F_TARGET = 0.20, 0.99, 0.60
SOFT_F0_DEFAULT = 0.05

# klass -> (sweep_type, MODE); neutral has no SLiM source
KLASS_MAP = {
    "neutral":     (None, None),
    "hard":        ("hard", "completed"),
    "soft":        ("soft", "completed"),
    "partialHard": ("hard", "partial"),
    "partialSoft": ("soft", "partial"),
}

MANIFEST_HEADER = ["klass", "file", "subwin", "is_focal", "seed", "s", "F0",
                   "f0_realized", "mode", "present_freq", "onset_tick",
                   "nrestart", "n_snp"]


def _fmt(x):
    """Manifest cell: empty for None, full-precision repr for floats/ints, else str."""
    if x is None or x == "":
        return ""
    if isinstance(x, float):
        return repr(x)
    return str(x)


def tau_gen_for_s(s, q=Q, n0=N0, tgrow=TGROW, tau_floor_gen=TAU_GEN, margin=2.0):
    """Sweep onset (generations before present) scaled so a sweep of coefficient s
    can FIX before the present -- i.e. weaker sweeps are simply OLDER.

    A hard sweep's logistic fixation time is ~ (2 / (s*q)) * ln(2 * n0/q) rescaled
    ticks; the onset is set to `margin`x that, floored at ``tau_floor_gen`` (strong
    sweeps stay recent) and capped at 70% of the growth epoch so it still fits. This
    keeps the full wide weak->strong s prior while GUARANTEEING completion: a fixed
    TAU_GEN=20000 (200 ticks) fails for ~half of log-U(1e-4,1e-2) -- any s below
    ~0.0012 needs >200 ticks to fix and restart-exhausts (SLiM exit 1). Soft sweeps
    start from f0>1/2N so fix faster; this hard-based onset over-margins them safely.
    """
    n0s = n0 / q
    sel_s = s * q
    dur_ticks = (2.0 / sel_s) * float(np.log(2.0 * n0s))
    tau_ticks = max(tau_floor_gen / q, margin * dur_ticks)
    tau_ticks = min(tau_ticks, 0.70 * (tgrow / q))     # leave room in the growth epoch
    return tau_ticks * q


def _build_slim_args(seed, sweep_type, mode, s, F0, L, n_dip):
    """A run_slim_sweep-compatible argparse.Namespace for run_slim + reconstruct.

    TAU_GEN is scaled to s (:func:`tau_gen_for_s`) so the sweep completes.
    """
    return argparse.Namespace(
        seed=seed, Q=Q, N0=N0, NREF=NREF, TGROW=TGROW, MU=MU, R=R,
        L=L, SEL=s, TAU_GEN=tau_gen_for_s(s), sweep_type=sweep_type, F0=F0,
        FIX_THRESH=FIX_THRESH, MODE=mode, F_LO=F_LO, F_HI=F_HI,
        F_TARGET=F_TARGET, n_dip=n_dip, slim=SLIM,
        slim_recipe=os.path.join(_SLIM_DIR, f"{sweep_type}_sweep.slim"),
    )


def write_windows(klass, seed, windows, outdir, n_dip, *, mode="", s=None,
                  F0=None, sweep_info=None):
    """Route the 11 ascertained windows to per-class ms files + append manifest rows.

    ``windows``: list of ``(subwin, geno2, pos_rel2)`` in subwindow order.
    Routing: neutral -> all -> ``neutral.msOut.gz``; a sweep klass C -> the focal
    window (subwin 5) -> ``<C>.msOut.gz`` and the other 10 -> ``linked<C>.msOut.gz``.
    Each file is created (header) if absent, else appended. One manifest row per
    window (schema :data:`MANIFEST_HEADER`). Pure glue over the sim -- unit-tested
    with synthetic windows. Returns ``{"target_of","files","rows","manifest"}``.
    """
    os.makedirs(outdir, exist_ok=True)
    n_hap = 2 * n_dip
    sweep_info = sweep_info or {}

    if klass == "neutral":
        target_of = {sw: "neutral.msOut.gz" for sw, _, _ in windows}
    else:
        # Per-subwindow-position files so diploSHIC makeTrainingSets can read
        # <klass>_<i>.msOut.gz: window i has the sweep in subwindow i, so _5 is the
        # focal sweep class and _0.._4,_6.._10 are the linked class (the sweep sits in
        # a flanking subwindow). Lumping them (focal vs linked) loses the per-position
        # structure makeTrainingSets requires.
        target_of = {sw: f"{klass}_{sw}.msOut.gz" for sw, _, _ in windows}

    by_file = {}                              # fname -> [(geno2, pos_rel2), ...] in subwin order
    for subwin, geno2, pos_rel2 in windows:
        by_file.setdefault(target_of[subwin], []).append((geno2, pos_rel2))
    for fname, recs in by_file.items():
        path = os.path.join(outdir, fname)
        write_ms_gz(recs, path, n_hap=n_hap, seed=seed, append=os.path.exists(path))

    f0_realized = sweep_info.get("f0_realized", "")
    present_freq = sweep_info.get("present_freq", sweep_info.get("freq", ""))
    onset_tick = sweep_info.get("onset_tick", "")
    nrestart = sweep_info.get("nrestart", "")

    manifest_path = os.path.join(outdir, "manifest.csv")
    write_header = not os.path.exists(manifest_path)
    rows = []
    with open(manifest_path, "a", newline="") as fh:
        w = csv.writer(fh)
        if write_header:
            w.writerow(MANIFEST_HEADER)
        for subwin, geno2, _pos in windows:
            row = [klass, target_of[subwin], subwin, subwin == FOCAL_SUBWIN, seed,
                   _fmt(s), _fmt(F0), _fmt(f0_realized), mode or "",
                   _fmt(present_freq), _fmt(onset_tick), _fmt(nrestart),
                   int(geno2.shape[1])]
            w.writerow(row)
            rows.append(row)
    return {"target_of": target_of, "files": sorted(set(target_of.values())),
            "rows": rows, "manifest": manifest_path}


def _run_source(klass, seed, s, F0, L, n_dip, outdir, keep_slim_trees):
    """Produce the reconstructed n=700 tree sequence + SWEEP_RESULT provenance.

    neutral -> pure msprime; sweeps -> SLiM recipe (run_slim) then reconstruct
    (subsample -> recapitate -> overlay). Raw SLiM .trees are cleaned up (unless
    keep_slim_trees) so a campaign does not fill the disk. Returns (ts, sweep_info).
    """
    sweep_type, mode = KLASS_MAP[klass]
    if klass == "neutral":
        return simulate_neutral(seed, L=L, n_dip=n_dip), {}

    args = _build_slim_args(seed, sweep_type, mode, s, F0, L, n_dip)
    work = os.path.join(outdir, "_slim_work")
    os.makedirs(work, exist_ok=True)
    tag = f"{klass}_seed{seed}"
    slim_trees = os.path.join(work, f"{tag}_raw.trees")
    save_trees = os.path.join(work, f"{tag}_restart.trees")
    est_trees = os.path.join(work, f"{tag}_est.trees")   # soft only
    try:
        sweep_info = rss.run_slim(args, slim_trees, save_trees, est_trees=est_trees)
        ts = rss.reconstruct(args, slim_trees)
    finally:
        if not keep_slim_trees:
            for f in (slim_trees, save_trees, est_trees):
                if os.path.exists(f):
                    os.remove(f)
            try:
                os.rmdir(work)          # only succeeds if now empty
            except OSError:
                pass
    return ts, sweep_info


def run_generate(klass, seed, outdir, s=None, F0=None, L=2_200_000,
                 n_dip=350, model=None, keep_slim_trees=True, min_maf=None):
    """Run one simulation end-to-end; write the routed ms files + manifest rows.

    Returns a summary dict (per-window n_snp, drawn s, mode, files) for callers /
    tests. See the module docstring for the class->source mapping, routing, the
    +-1 bp beneficial exclusion, and the seed-derivation of the s-draw / RNG.
    """
    if klass not in KLASS_MAP:
        raise ValueError(f"unknown klass {klass!r}; choose from {list(KLASS_MAP)}")
    sweep_type, mode = KLASS_MAP[klass]
    os.makedirs(outdir, exist_ok=True)

    # ---- deterministic sub-streams from --seed ----------------------------
    s_ss, asc_ss = np.random.SeedSequence(seed).spawn(2)
    s_rng = np.random.default_rng(s_ss)
    asc_rng = np.random.default_rng(asc_ss)

    # ---- resolve the selection / standing-frequency knobs -----------------
    if sweep_type is None:                    # neutral: no selection knobs
        s, F0 = None, None
    else:
        if s is None:                         # log-uniform s ~ U(1e-4, 1e-2)
            s = float(10.0 ** s_rng.uniform(np.log10(S_LO), np.log10(S_HI)))
        if sweep_type == "soft":
            if F0 is None:
                F0 = SOFT_F0_DEFAULT
        else:                                 # hard / partialHard never take F0
            F0 = None

    # ---- 1. simulate + reconstruct (or msprime neutral) -------------------
    ts, sweep_info = _run_source(klass, seed, s, F0, L, n_dip, outdir,
                                 keep_slim_trees)
    n_hap = 2 * n_dip
    if ts.num_samples != n_hap:
        raise RuntimeError(f"reconstructed {ts.num_samples} haplotypes, expected {n_hap}")
    if keep_slim_trees:                       # persist the EXACT sample tree seq for later re-ascertainment
        ts.dump(os.path.join(outdir, f"{klass}_seed{seed}_recon.trees"))

    # ---- 2. genotype matrix ONCE ------------------------------------------
    G, pos = genotype_and_positions(ts)       # (n_site, n_hap), (n_site,)

    # ---- 3. exclude the beneficial site (+-1 bp of XPOS) before windowing --
    xpos = L // 2
    keep_site = np.abs(pos - xpos) > 1.0
    Gf, posf = G[keep_site, :], pos[keep_site]

    # ---- 4. slice the 11 windows + ascertain each -------------------------
    if model is None:
        model = load_model()
    if min_maf is not None:                   # e.g. --min-maf 0 for a full-SFS training set
        model = dict(model); model["min_maf"] = float(min_maf)
    offsets = window_offsets(win=WIN, sub=SUB, n_sub=N_SUB, sweep_pos=xpos)
    windows = []                              # (subwin, geno2, pos_rel2)
    for subwin, offset in offsets:
        geno, pos_rel = window_from(Gf, posf, offset, WIN)
        geno2, pos_rel2 = ascertain(geno, pos_rel, asc_rng, model)
        windows.append((subwin, geno2, pos_rel2))

    # ---- 5 + 6. route + write per-class ms files + append manifest rows ----
    written = write_windows(klass, seed, windows, outdir, n_dip, mode=mode,
                            s=s, F0=F0, sweep_info=sweep_info)

    summary = {
        "klass": klass, "seed": seed, "mode": mode or "", "s": s, "F0": F0,
        "f0_realized": sweep_info.get("f0_realized", ""),
        "present_freq": sweep_info.get("present_freq", sweep_info.get("freq", "")),
        "onset_tick": sweep_info.get("onset_tick", ""),
        "nrestart": sweep_info.get("nrestart", ""),
        "n_hap": n_hap, "n_sites_total": int(ts.num_sites),
        "files": written["files"],
        "n_snp_by_subwin": {int(sw): int(g.shape[1]) for sw, g, _ in windows},
        "manifest": written["manifest"], "outdir": outdir,
    }
    print("[generate] klass=%s seed=%s mode=%s s=%s F0=%s -> %s"
          % (klass, seed, mode or "-", s, F0, summary["files"]), flush=True)
    print("[generate] n_snp by subwin: "
          + ", ".join(f"{sw}:{summary['n_snp_by_subwin'][sw]}" for sw in range(N_SUB)),
          flush=True)
    return summary


def parse_args(argv=None):
    p = argparse.ArgumentParser(
        description="Single-sim end-to-end diploSHIC/partialSHIC training-window generator.")
    p.add_argument("--klass", required=True, choices=list(KLASS_MAP),
                   help="training class ('class' is a Python keyword; use --klass).")
    p.add_argument("--seed", type=int, required=True,
                   help="master seed; everything (SLiM, s-draw, ascertainment) is derived from it.")
    p.add_argument("--s", type=float, default=None,
                   help="selection coefficient (unscaled). If omitted for a sweep, drawn "
                        "log-uniform in [1e-4, 1e-2] deterministically from --seed.")
    p.add_argument("--F0", type=float, default=None,
                   help="soft-sweep target standing frequency (default 0.05 for soft/partialSoft; "
                        "ignored for hard/partialHard/neutral).")
    p.add_argument("--outdir", required=True, help="output directory (absolute).")
    p.add_argument("--L", type=int, default=2_200_000, help="simulated locus length (bp).")
    p.add_argument("--n-dip", type=int, default=350, dest="n_dip",
                   help="diploids sampled at present (=> 2x haplotypes; default 350 -> 700 hap).")
    p.add_argument("--keep-slim-trees", action=argparse.BooleanOptionalAction, default=True,
                   help="keep the raw SLiM .trees AND dump the reconstructed n=700 tree sequence "
                        "(<klass>_seed<seed>_recon.trees -- the one you re-ascertain from). "
                        "Default ON; pass --no-keep-slim-trees to clean them up.")
    p.add_argument("--min-maf", type=float, default=None, dest="min_maf",
                   help="override the ascertainment MAF cutoff (0 -> full-SFS training set; "
                        "default: the config's min_maf, currently 0.01).")
    p.add_argument("--model", default=None,
                   help="ascertainment model json (default: step-13 ascertainment.json).")
    return p.parse_args(argv)


def main(argv=None):
    a = parse_args(argv)
    model = load_model(a.model) if a.model else None
    run_generate(a.klass, a.seed, a.outdir, s=a.s, F0=a.F0, L=a.L,
                 n_dip=a.n_dip, model=model, keep_slim_trees=a.keep_slim_trees,
                 min_maf=a.min_maf)


if __name__ == "__main__":
    main()
