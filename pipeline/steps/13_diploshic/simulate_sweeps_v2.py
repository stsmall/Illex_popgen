"""
simulate_sweeps_v2.py — Task 6b rebuild of diploSHIC training sweep sims.

Generates:
  Center (window 5): 2000 hard + 2000 soft
    sweep -x jittered uniformly in [5/11, 6/11]
    files: hard_w5_r####.ms, soft_w5_r####.ms
  Linked (windows 0-4, 6-10): 200 hard + 200 soft each (2000 hard + 2000 soft total)
    sweep -x = (pos+0.5)/11 (fixed per window)
    files: hard_w<pos>_r####.ms, soft_w<pos>_r####.ms

All sims use demography.discoal_demog() growth flags.
Sweep priors: alpha=2*N0*10^U(-4,-2), tau~U(0,0.05), soft freq=10^U(-3,-1).
Parallelism: P8 (shared machine constraint).

Output raw ms: WD/sims/sweep2/raw/
"""

import sys, os, subprocess, numpy as np
from multiprocessing import Pool

WD = "/sietch_colab/data_share/illex/popgen_data/analysis/steps/13_diploshic"
sys.path.insert(0, WD)
import demography as D
import simulate as S  # for rho_draw, NHAP, WINLEN, theta()

DISCOAL = "/home/ssmall/programs/discoal_oobfix/discoal"
OUTDIR = f"{WD}/sims/sweep2/raw"
NHAP = S.NHAP      # 700
WINLEN = S.WINLEN  # 44000
NSUB = S.NSUB      # 11

os.makedirs(OUTDIR, exist_ok=True)

# ─────────────────────────────────────────────────────────────
# Build sim spec list: (outpath, seed, kind, window_pos, jitter)
# seed is deterministic from (kind, window, rep)
# ─────────────────────────────────────────────────────────────

def make_spec(kind, win, rep):
    """Return (outpath, seed, kind, win) for one simulation."""
    outpath = f"{OUTDIR}/{kind}_w{win}_r{rep:04d}.ms"
    # deterministic seed: use hash of tuple to keep distinct from old sims
    seed = hash(("v2", kind, win, rep)) % (2**31)
    return (outpath, seed, kind, win)


specs = []

# Center window 5: 2000 hard + 2000 soft
for rep in range(2000):
    specs.append(make_spec("hard", 5, rep))
    specs.append(make_spec("soft", 5, rep))

# Linked windows 0-4, 6-10: 200 hard + 200 soft each
for win in list(range(5)) + list(range(6, 11)):
    for rep in range(200):
        specs.append(make_spec("hard", win, rep))
        specs.append(make_spec("soft", win, rep))

print(f"Total specs to generate: {len(specs)}")


# ─────────────────────────────────────────────────────────────
# Worker: generate one sim
# ─────────────────────────────────────────────────────────────

def run_one(spec):
    outpath, seed, kind, win = spec
    if os.path.exists(outpath):
        return "skip"

    rng = np.random.default_rng(seed)

    # rho
    rho = S.rho_draw(rng)
    theta = S.theta()

    # sweep position
    if win == 5:
        # jitter uniformly in [5/11, 6/11]
        x = rng.uniform(5.0/NSUB, 6.0/NSUB)
    else:
        x = (win + 0.5) / NSUB

    # sweep strength: alpha = 2*N0 * s, s ~ 10^U(-4,-2)
    alpha = 2.0 * D.N0 * 10 ** rng.uniform(-4, -2)
    tau = rng.uniform(0, 0.05)

    # base discoal args
    growth_flags = D.discoal_demog()
    cmd = [DISCOAL, str(NHAP), "1", str(WINLEN),
           "-t", f"{theta:.3f}",
           "-r", f"{rho:.3f}",
           "-ws", f"{tau:.4f}",
           "-a", f"{alpha:.2f}",
           "-x", f"{x:.4f}"]

    if kind == "soft":
        freq = 10 ** rng.uniform(-3, -1)
        cmd += ["-f", f"{freq:.5f}"]

    cmd += growth_flags

    try:
        result = subprocess.run(cmd, capture_output=True, text=True, check=True)
        with open(outpath, "w") as f:
            f.write(result.stdout)
        return "ok"
    except subprocess.CalledProcessError as e:
        return f"FAIL {outpath}: {e.stderr[:200]}"


# ─────────────────────────────────────────────────────────────
# Main: run at P8
# ─────────────────────────────────────────────────────────────

if __name__ == "__main__":
    print(f"Running {len(specs)} sims at P8...")
    with Pool(processes=8) as pool:
        results = pool.map(run_one, specs)

    ok = sum(1 for r in results if r == "ok")
    skip = sum(1 for r in results if r == "skip")
    fail = [r for r in results if r not in ("ok", "skip")]
    print(f"Done: {ok} new, {skip} skipped, {len(fail)} failures")
    if fail:
        for f in fail[:10]:
            print(f"  {f}")
        sys.exit(1)
    print("All sims complete.")
