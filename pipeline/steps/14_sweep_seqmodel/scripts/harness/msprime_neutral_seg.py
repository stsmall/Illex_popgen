#!/usr/bin/env python3
"""msprime_neutral_seg.py -- FAST pure-msprime neutral segment under the Illex growth demography,
with the REAL ReLERNN recombination landscape (RateMap) OR a constant rate. For the demographic-
model-adequacy check (user 2026-08-11): do neutral sims reproduce the real genome's fvec summary-stat
DISTRIBUTION? If not even variable-recomb msprime matches, the demographic model is wrong.

No SLiM, no rescaling -- exact coalescent, seconds per segment. Output a mutated tree sequence
(.trees, same interface calib_scan_segment.py consumes) so it goes straight through the empirical
fvecVcf path (ascertain -> unphased VCF -> fvecVcf).

Usage: msprime_neutral_seg.py --chrom 43 --start 60000000 --L 3000000 --seed S --n-dip 350
       [--recomb variable|constant] --outdir DIR
"""
import argparse, os, subprocess, sys
import numpy as np, msprime

HERE = os.path.dirname(os.path.abspath(__file__)); SLIMDIR = os.path.join(HERE, "slim")
N0, NREF, TGROW, MU, RCONST = 6_808_096.0, 547_928.0, 769_519.0, 3e-9, 2.1e-9
R = "/sietch_colab/data_share/illex/popgen_data/analysis/steps/11_relernn"


def growth_demography():
    g = np.log(N0 / NREF) / TGROW
    dem = msprime.Demography()
    dem.add_population(name="pop", initial_size=N0, growth_rate=g)
    dem.add_population_parameters_change(time=TGROW, initial_size=NREF, growth_rate=0)
    return dem


def relernn_ratemap(chrom, start, end, outdir):
    """sex-averaged ReLERNN map for the segment -> msprime.RateMap (UNscaled per-bp rate)."""
    mapfile = os.path.join(outdir, f"seg_{chrom}_{start}.slimmap")
    subprocess.run([sys.executable, os.path.join(SLIMDIR, "relernn_to_slim_map.py"),
                    "--male", f"{R}/run_male_auto/proj/male.kept.PREDICT.BSCORRECTED.txt",
                    "--female", f"{R}/run_female_auto/proj/female.kept.PREDICT.BSCORRECTED.txt",
                    "--chrom", str(chrom), "--start", str(start), "--end", str(end),
                    "--out", mapfile], check=True)
    ends, rates = [], []
    for line in open(mapfile):
        if not line.strip():
            continue
        e, r = line.split(); ends.append(int(e)); rates.append(float(r))
    L = ends[-1] + 1
    positions = [0] + [e + 1 for e in ends]
    return msprime.RateMap(position=positions, rate=rates), L


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--chrom", required=True); ap.add_argument("--start", type=int, required=True)
    ap.add_argument("--L", type=int, required=True); ap.add_argument("--seed", type=int, required=True)
    ap.add_argument("--n-dip", type=int, default=350)
    ap.add_argument("--recomb", choices=["variable", "constant"], default="variable")
    ap.add_argument("--outdir", required=True)
    a = ap.parse_args(); os.makedirs(a.outdir, exist_ok=True)
    if a.recomb == "variable":
        rmap, Lm = relernn_ratemap(a.chrom, a.start, a.start + a.L, a.outdir)
        assert Lm == a.L, f"map L {Lm} != {a.L}"
        rec = rmap
    else:
        rec = RCONST
    ts = msprime.sim_ancestry(samples={"pop": a.n_dip}, demography=growth_demography(),
                              sequence_length=a.L, recombination_rate=rec, random_seed=a.seed)
    ts = msprime.sim_mutations(ts, rate=MU, random_seed=a.seed + 1)
    tag = f"msp_{a.recomb}_{a.chrom}_{a.start}_s{a.seed}"
    out = os.path.join(a.outdir, tag + ".recap.trees")     # .recap.trees name so calib_scan_segment reads it
    ts.dump(out)
    print(f"[msp] {a.recomb} {a.chrom}:{a.start}-{a.start+a.L} {ts.num_samples} haps, {ts.num_sites} sites, "
          f"pi={ts.diversity(mode='site'):.5f} -> {out}", flush=True)


if __name__ == "__main__":
    main()
