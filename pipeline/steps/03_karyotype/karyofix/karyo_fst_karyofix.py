"""Between-karyotype Hudson FST and near-fixed differences with the karyotype_final groups (2a/2a n=223, 2+/2+ n=97).

Reconstructs the recipe of steps/03_karyotype/inversion_content/chr2_karyo_fst_windows.tsv (its streaming awk pass is
not saved; this reproduces that table's num/den/n to <1e-6 from aa_chr2_af.txt / bb_chr2_af.txt, see --check):
  per SNP of variants_filt.vcf.gz (chr2), p = bcftools +fill-tags AF within each group (called genotypes), and
    num = (p1 - p2)^2 - p1(1-p1)/(N1-1) - p2(1-p2)/(N2-1),  den = p1(1-p2) + p2(1-p1),  N = 2 x group size (fixed),
  summed in 100-kb windows (win = floor(pos / 1e5) * 1e5) over SNPs with den > 0; columns win, num, den, n.
  Plot rule (plot_chr2_karyotype_fst.py): keep n >= 10 and den > 0, FST = max(num/den, 0), window mid = win + 50 kb,
  inversion window if 60.54 Mb <= mid <= 79.50 Mb; medians inside / outside.
Near-fixed differences (analyze_inversion_fixed_diffs.py): |p1 - p2| >= 0.9 at SNPs with 60.54 Mb <= pos <= 79.50 Mb.
"""
import sys

import numpy as np
import pandas as pd

INV_S, INV_E = 60_540_000, 79_500_000
WIN = 100_000


def load(fa, fb):
    a = pd.read_csv(fa, sep="\t", header=None, names=["c", "pos", "p1"])
    b = pd.read_csv(fb, sep="\t", header=None, names=["c", "pos", "p2"])
    assert (a.pos.values == b.pos.values).all()
    return a.pos.values, pd.to_numeric(a.p1, errors="coerce").values, pd.to_numeric(b.p2, errors="coerce").values


def run(fa, fb, n1, n2, out=None):
    pos, p1, p2 = load(fa, fb)
    N1, N2 = 2 * n1, 2 * n2
    num = (p1 - p2) ** 2 - p1 * (1 - p1) / (N1 - 1) - p2 * (1 - p2) / (N2 - 1)
    den = p1 * (1 - p2) + p2 * (1 - p1)
    ok = np.isfinite(num) & np.isfinite(den) & (den > 0)
    w = (pd.DataFrame({"win": (pos[ok] // WIN) * WIN, "num": num[ok], "den": den[ok]})
         .groupby("win").agg(num=("num", "sum"), den=("den", "sum"), n=("num", "size")).reset_index())
    if out:
        # the original table holds only windows with n >= 10 (identical to it for the old groups)
        w[w.n >= 10].to_csv(out, sep="\t", header=False, index=False, float_format="%.6f")
    k = w[(w.n >= 10) & (w.den > 0)].copy()
    k["fst"] = (k.num / k.den).clip(lower=0)
    mid = k.win + WIN / 2
    inv = (mid >= INV_S) & (mid <= INV_E)
    res = {"n_snps": int(np.isfinite(p1).sum()), "nan_af": int((~np.isfinite(p1) | ~np.isfinite(p2)).sum())}
    for lab, lo, hi in (("inv_60.54_79.50", INV_S, INV_E), ("body_60.55_79.52", 60_550_000, 79_520_000),
                        ("chr2_60_80", 60_000_000, 80_000_000)):
        s = ok & (pos >= lo) & (pos <= hi)
        res[f"fst_ratio_of_averages_{lab}"] = float(num[s].sum() / den[s].sum())
    res["fst_window_median_inversion"] = float(k.fst[inv].median())
    res["fst_window_median_collinear"] = float(k.fst[~inv].median())
    res["n_windows_inversion"], res["n_windows_collinear"] = int(inv.sum()), int((~inv).sum())
    s = (pos >= INV_S) & (pos <= INV_E)
    d = np.abs(p1 - p2)[s]
    res["inversion_snps"] = int(s.sum())
    for t in (0.8, 0.9, 0.95, 0.99):
        res[f"near_fixed_dAF_ge_{t}"] = int((d >= t).sum())
    return res


if __name__ == "__main__":
    import json

    fa, fb, n1, n2 = sys.argv[1], sys.argv[2], int(sys.argv[3]), int(sys.argv[4])
    out = sys.argv[5] if len(sys.argv) > 5 else None
    print(json.dumps(run(fa, fb, n1, n2, out), indent=1))
