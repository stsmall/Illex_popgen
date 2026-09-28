#!/usr/bin/env python3
"""calib_run_sim.py -- one CONTIGUOUS calibration segment (neutral|hard|soft) under a VARIABLE
ReLERNN recombination map, for the HMM calibration harness (PLAN_02 Phase A3).

Reuses the training machinery's approach (run SLiM recipe -> subsample -> recapitate -> overlay
neutral mutations) but with TWO changes: (1) the segment recipe (calib_*_segment.slim / neutral_
segment.slim) reads a per-segment SLiM map; (2) recapitation uses an msprime.RateMap (variable
map x Q) instead of a constant rate, so the deep ancestry matches the forward-sim map.

Output: the recapitated+overlaid tree sequence (.trees) + a VCF (for the empirical-matched fvecVcf
scan) OR left as ts for B's fvecSim tiling. Deterministic from --seed.
"""
import argparse, os, subprocess, sys, warnings
import numpy as np
import msprime, tskit, pyslim

HERE = os.path.dirname(os.path.abspath(__file__))
SLIMDIR = os.path.join(HERE, "slim")
SLIM = os.environ.get("SLIM_BIN", "/home/ssmall/miniforge3/envs/slim_sim/bin/slim")
# demography / rescaling (verbatim from the recipes)
Q, N0, NREF, TGROW, MU = 100.0, 6808096.0, 547928.0, 769519.0, 3e-9

RECIPE = {"neutral": "neutral_segment.slim",
          "hard":    "calib_hard_segment.slim",
          "soft":    "calib_soft_segment.slim"}


def build_map(chrom, start, end, out):
    """Call the validated A1 converter -> SLiM map file (end<TAB>rate) for the segment."""
    R = "/sietch_colab/data_share/illex/popgen_data/analysis/steps/11_relernn"
    subprocess.run([sys.executable, os.path.join(SLIMDIR, "relernn_to_slim_map.py"),
                    "--male", f"{R}/run_male_auto/proj/male.kept.PREDICT.BSCORRECTED.txt",
                    "--female", f"{R}/run_female_auto/proj/female.kept.PREDICT.BSCORRECTED.txt",
                    "--chrom", str(chrom), "--start", str(start), "--end", str(end),
                    "--out", out], check=True)


def ratemap_from_file(mapfile):
    """MAPFILE rows 'end<TAB>rate' -> msprime.RateMap over [0,L] with rate*Q (recap in rescaled units)."""
    ends, rates = [], []
    with open(mapfile) as fh:
        for line in fh:
            if not line.strip():
                continue
            e, r = line.split()
            ends.append(int(e)); rates.append(float(r) * Q)
    L = ends[-1] + 1
    positions = [0] + [e + 1 for e in ends]          # interval boundaries, len = len(rates)+1
    return msprime.RateMap(position=positions, rate=rates), L


def run_recipe(klass, seed, L, mapfile, outtrees, sel=None, tau=None, f0=None):
    d = {"Q": Q, "N0": N0, "NREF": NREF, "TGROW": TGROW, "MU": MU, "L": L,
         "MAPFILE": f"'{mapfile}'", "OUTPATH": f"'{outtrees}'"}
    if klass != "neutral":
        d["SEL"] = sel; d["XPOS"] = L // 2; d["SAVEPATH"] = f"'{outtrees}.restart'"
        if tau is not None: d["TAU_GEN"] = tau
        if f0 is not None:  d["F0"] = f0
    cmd = [SLIM, "-s", str(seed)]
    for k, v in d.items():
        cmd += ["-d", f"{k}={v}"]
    cmd.append(os.path.join(SLIMDIR, RECIPE[klass]))
    r = subprocess.run(cmd, capture_output=True, text=True)
    if r.returncode != 0 or not os.path.exists(outtrees):
        sys.stderr.write(r.stdout[-2000:] + "\n" + r.stderr[-2000:] + "\n")
        raise RuntimeError(f"SLiM {klass} failed (rc={r.returncode})")
    return r.stdout


def reconstruct(trees, ratemap, n_dip, seed):
    """subsample n_dip diploids at present -> simplify -> recapitate (VARIABLE RateMap) -> overlay."""
    rng = np.random.default_rng(seed)
    ts = tskit.load(trees)
    alive = pyslim.individuals_alive_at(ts, 0)
    if len(alive) < n_dip:
        raise RuntimeError(f"only {len(alive)} alive; need {n_dip}")
    chosen = rng.choice(alive, size=n_dip, replace=False)
    nodes = np.concatenate([ts.individual(i).nodes for i in chosen])
    ts = ts.simplify(samples=nodes, keep_input_roots=True)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", msprime.TimeUnitsMismatchWarning)
        ts = pyslim.recapitate(ts, ancestral_Ne=NREF / Q, recombination_rate=ratemap,
                               random_seed=seed + 1)                     # <-- variable map
        ts = msprime.sim_mutations(ts, rate=MU * Q, random_seed=seed + 2)
    return ts


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--klass", required=True, choices=list(RECIPE))
    ap.add_argument("--chrom", required=True)
    ap.add_argument("--start", type=int, required=True)
    ap.add_argument("--L", type=int, required=True)
    ap.add_argument("--seed", type=int, required=True)
    ap.add_argument("--n-dip", type=int, default=350)
    ap.add_argument("--sel", type=float, default=0.01)
    ap.add_argument("--outdir", required=True)
    a = ap.parse_args()
    os.makedirs(a.outdir, exist_ok=True)
    tag = f"{a.klass}_{a.chrom}_{a.start}_s{a.seed}"
    mapfile = os.path.join(a.outdir, tag + ".slimmap")
    trees = os.path.join(a.outdir, tag + ".trees")
    build_map(a.chrom, a.start, a.start + a.L, mapfile)
    ratemap, Lmap = ratemap_from_file(mapfile)
    assert Lmap == a.L, f"map L {Lmap} != {a.L}"
    print(f"[calib] SLiM {a.klass} L={a.L} seg {a.chrom}:{a.start}-{a.start + a.L}", flush=True)
    run_recipe(a.klass, a.seed, a.L, mapfile, trees, sel=a.sel)
    ts = reconstruct(trees, ratemap, a.n_dip, a.seed)
    outts = os.path.join(a.outdir, tag + ".recap.trees")
    ts.dump(outts)
    thetaW = ts.num_sites / sum(1.0 / i for i in range(1, ts.num_samples))
    print(f"[calib] DONE {tag}: {ts.num_samples} haps, {ts.num_sites} sites, "
          f"pi={ts.diversity(mode='site'):.5f}, thetaW/site={thetaW / a.L:.5f} -> {outts}", flush=True)


if __name__ == "__main__":
    main()
