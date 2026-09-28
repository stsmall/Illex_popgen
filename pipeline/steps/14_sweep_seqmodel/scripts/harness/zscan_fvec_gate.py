#!/usr/bin/env python3
"""zscan_fvec_gate.py -- decide whether subsampling the autosomal sims to n=330 is an adequate
model for the lower-Ne (3/4 autosomal) chrZ, or whether Z-demography re-sims are needed.

Rationale (validated on empirical windowed_diversity): the Z theta LEVEL is ~2/3 of autosomes
(pi 0.0057 vs 0.0084), but diploSHIC normalizes each stat's 11-subwindow vector WITHIN each
window, so a uniform theta scaling cancels -- the CNN keys on the RELATIVE spatial profile.
And the one non-absorbed feature, Tajima's D, is identical Z-vs-auto (-2.10 vs -2.09).

KEY: an ABSOLUTE in-distribution test is the WRONG gate -- EVERY empirical chromosome sits
partly outside the flat-sim cloud because of the known within-window spatial-heterogeneity gap
(the autosomal over-call story; distVar ~0.33 out even for chr1). We already ACCEPTED that for
the autosomes (outlier-ranking framing). So the gate is COMPARATIVE: it measures the chrZ-vs-simZ
mismatch AND the autosome-vs-simAuto baseline mismatch, and flags only stats where chrZ is
MATERIALLY WORSE than the autosomes -- i.e. Z-SPECIFIC degradation introduced by the 3/4-Ne /
subsampling, over and above the shared baseline we already live with.

Usage: zscan_fvec_gate.py --z-sim-dir D --z-emp chrZ.fvec --auto-sim-dir D --auto-emp chr1.fvec --out V
Exit 0 = PASS (Z no worse than autosomes), 3 = FLAG (Z-specific offset -> re-sim under 3/4 Ne).
"""
import argparse, glob, os, sys
import numpy as np

# LD/haplotype stats are known-degenerate on UNPHASED data (diploshic-zns-train-scan-mismatch,
# illex-unphased-no-haplotype-stats) -> reported but NOT used for the pass/fail decision.
LD_STATS = {"diplo_ZnS", "diplo_Omega", "diplo_H1", "diplo_H12", "diplo_H2/H1", "diplo_H", "nDiplos"}
FLAG_DELTA = 0.15  # per-stat: FLAG if (chrZ out-frac) - (autosome out-frac) exceeds this = Z-specific

def load_fvec(path, strip_meta):
    with open(path) as fh:
        hdr = fh.readline().rstrip("\n").split("\t")
    # skip the 4 non-numeric meta cols (incl bigWinRange "a-b") DURING load, not after
    usecols = range(4, len(hdr)) if strip_meta else None
    data = np.loadtxt(path, delimiter="\t", skiprows=1, ndmin=2, usecols=usecols)
    if strip_meta:
        hdr = hdr[4:]
    return hdr, data

def stat_groups(hdr):
    """map stat-name -> column indices, from '<stat>_win\\d+' headers (stat-major layout)."""
    import re
    g = {}
    for i, name in enumerate(hdr):
        s = re.sub(r"_win\d+$", "", name)
        g.setdefault(s, []).append(i)
    return g

def l1_normalize_per_window(mat, cols):
    """within-window relative profile of one stat: each row's cols / sum|cols| (diploSHIC-style)."""
    sub = mat[:, cols].astype(float)
    denom = np.abs(sub).sum(axis=1, keepdims=True)
    denom[denom == 0] = 1.0
    return sub / denom

def load_sim_cloud(sim_dir):
    files = sorted(glob.glob(os.path.join(sim_dir, "*.fvec")))
    if not files:
        sys.exit(f"no sim fvecs in {sim_dir}")
    hdr = None; parts = []
    for f in files:
        h, d = load_fvec(f, strip_meta=False)
        hdr = hdr or h; parts.append(d)
    return hdr, np.vstack(parts)

def out_frac_per_stat(sim, emp, cols):
    """frac of empirical (win x subwin) normalized values outside the sim [0.5,99.5] envelope."""
    ns = l1_normalize_per_window(sim, cols)
    ne = l1_normalize_per_window(emp, cols)
    lo = np.percentile(ns, 0.5, axis=0); hi = np.percentile(ns, 99.5, axis=0)
    return float(((ne < lo) | (ne > hi)).mean())

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--z-sim-dir", required=True)      # Z330 fvec_nometa/
    ap.add_argument("--z-emp", required=True)           # empirical chrZ.fvec
    ap.add_argument("--auto-sim-dir", required=True)    # autosomal (n=350) fvec_nometa/  [baseline]
    ap.add_argument("--auto-emp", required=True)        # an autosomal empirical fvec [baseline]
    ap.add_argument("--out", required=True)
    a = ap.parse_args()

    zhdr, zsim = load_sim_cloud(a.z_sim_dir)
    ahdr, asim = load_sim_cloud(a.auto_sim_dir)
    zehdr, zemp = load_fvec(a.z_emp, strip_meta=True)
    aehdr, aemp = load_fvec(a.auto_emp, strip_meta=True)
    assert zhdr == zehdr == ahdr == aehdr, "column mismatch across inputs"

    groups = stat_groups(zhdr)
    lines = []
    lines.append(f"# chrZ subsample-adequacy gate (COMPARATIVE)")
    lines.append(f"# Z: sim n_win={zsim.shape[0]} emp n_win={zemp.shape[0]} | AUTO baseline: sim n_win={asim.shape[0]} emp n_win={aemp.shape[0]}")
    lines.append(f"# out-frac = frac of empirical windows outside its OWN sim [0.5,99.5] envelope (within-window normalized)")
    lines.append(f"# FLAG only if (Z_out - AUTO_out) > {FLAG_DELTA} on a NON-LD stat = Z-SPECIFIC degradation beyond the accepted autosomal baseline")
    lines.append("")
    lines.append(f"{'stat':<14}{'Z_lvl(emp/sim)':>15}{'Z_out':>8}{'AUTO_out':>10}{'delta':>8}  {'used?':>6}  flag")
    decisive = []
    for s, cols in groups.items():
        z_ratio = (np.median(zemp[:, cols]) / np.median(zsim[:, cols])) if np.median(zsim[:, cols]) != 0 else float("nan")
        z_out = out_frac_per_stat(zsim, zemp, cols)
        a_out = out_frac_per_stat(asim, aemp, cols)
        delta = z_out - a_out
        used = "no(LD)" if s in LD_STATS else "yes"
        flag = ""
        if used == "yes" and delta > FLAG_DELTA:
            flag = "<<< Z-WORSE"; decisive.append((s, delta))
        lines.append(f"{s:<14}{z_ratio:>15.3f}{z_out:>8.3f}{a_out:>10.3f}{delta:>+8.3f}  {used:>6}  {flag}")
    lines.append("")
    if decisive:
        lines.append("VERDICT: FLAG -- chrZ materially WORSE than the autosomal baseline on: " +
                     ", ".join(f"{s}(+{d:.2f})" for s, d in decisive))
        lines.append("  -> subsampling+3/4-Ne introduces Z-SPECIFIC feature mismatch beyond the accepted autosomal gap.")
        lines.append("  -> regenerate chrZ sims under 3/4-Ne demography + Z recomb map before trusting chrZ calls.")
        verdict = 3
    else:
        lines.append("VERDICT: PASS -- chrZ is NO WORSE than the autosomes we already accepted (no non-LD stat exceeds")
        lines.append(f"  the +{FLAG_DELTA} delta). The theta LEVEL differs (~2/3, see Z_lvl) but is absorbed by")
        lines.append("  within-window normalization; the residual mismatch is the SAME spatial-het gap as the autosomes.")
        lines.append("  The n=330-subsampled model is adequate for chrZ; use the same outlier-ranking framing as autosomes.")
        verdict = 0
    txt = "\n".join(lines) + "\n"
    open(a.out, "w").write(txt)
    sys.stdout.write(txt)
    sys.exit(verdict)

if __name__ == "__main__":
    main()
