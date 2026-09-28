#!/usr/bin/env python3
"""Fig 4 — single composite recombination figure.

A: genome-wide autosomal landscape, SEX-AVERAGED (male+female pooled -> rolling
   median), all autosomes concatenated, chr2 inserted (from the autosomal-net
   maps applied to chr2), chr2 inversion 60.54-79.50 Mb shaded.
B: chrZ, male vs female kept separate (genuinely sex-biased).

Heterochiasmy is not identifiable on autosomes (both sexes' LD-based maps
estimate the sex-averaged rate), so autosomes are combined; the Z is not.
"""
import glob, os, re, sys
import numpy as np, pandas as pd
sys.path.insert(0, "/sietch_colab/data_share/illex/popgen_data/analysis/manuscript")
from figstyle import apply, C, SEXES, despine
apply()
import matplotlib.pyplot as plt

RD = "/sietch_colab/data_share/illex/popgen_data/analysis/steps/11_relernn"
OUT = "/sietch_colab/data_share/illex/popgen_data/analysis/manuscript/figures/fig6_recombination.png"
CM = 1e8
INV = (60.54, 79.50)  # chr2 inversion (Mb)


def _load(path):
    if not path or not os.path.exists(path):
        return None
    df = pd.read_csv(path, sep="\t")
    df["chrom"] = df["chrom"].astype(str).str.replace(r"^b'|'$", "", regex=True)
    df["cM_Mb"] = df["recombRate"].astype(float) * CM
    df["mid"] = (df["start"] + df["end"]) / 2.0 / 1e6  # Mb
    return df[["chrom", "mid", "cM_Mb"]]


def _bscorr(projdir):
    hits = glob.glob(os.path.join(RD, projdir, "proj", "*.PREDICT.BSCORRECTED.txt"))
    return _load(hits[0]) if hits else None


def rollmed(x, y, k=25):
    o = np.argsort(x); x, y = np.asarray(x)[o], np.asarray(y)[o]
    if len(y) < 3:
        return x, y
    k = max(3, min(k, len(y) // 2 * 2 + 1))
    ys = pd.Series(y).rolling(k, center=True, min_periods=max(3, k // 3)).median().to_numpy()
    return x, ys


# ---- autosomes: pool male+female (sex-average) --------------------------
male = _bscorr("run_male_auto"); female = _bscorr("run_female_auto")
auto = pd.concat([d for d in (male, female) if d is not None], ignore_index=True)
# add chr2 from the autosomal-net predictions (male + female), which the
# autosomal ReLERNN run excluded
for f in ("chr2_male.autonet.PREDICT.txt", "chr2_female.autonet.PREDICT.txt"):
    d = _load(os.path.join(RD, "chr2_autosomal_predict", f))
    if d is not None:
        auto = pd.concat([auto, d], ignore_index=True)

def chrkey(c):
    return int(c) if c.isdigit() else 999
chroms = sorted(auto["chrom"].unique(), key=chrkey)

fig, (axA, axB) = plt.subplots(
    2, 1, figsize=(13, 6.2), gridspec_kw={"height_ratios": [3, 1], "hspace": 0.42})

# panel A: concatenated sex-averaged landscape
offset = 0.0; ticks = []; ticklab = []; GAP = 3.0
for i, c in enumerate(chroms):
    sub = auto[auto["chrom"] == c]
    x, y = rollmed(sub["mid"].to_numpy(), sub["cM_Mb"].to_numpy())
    gx = x + offset
    span = float(np.nanmax(sub["mid"]))
    if i % 2 == 0:
        axA.axvspan(offset, offset + span, color="#f4f4f4", zorder=0)
    if c == "2":  # shade the inversion
        axA.axvspan(offset + INV[0], offset + INV[1], color=C["shade"], zorder=1)
        axA.annotate("chr2 inversion", xy=(offset + (INV[0] + INV[1]) / 2, 0.30),
                     xytext=(offset + 40, 0.355), fontsize=8, color=C["warm"],
                     ha="left", arrowprops=dict(arrowstyle="-", color=C["warm"], lw=0.7))
    axA.plot(gx, y, color=C["ink"], lw=0.7)
    ticks.append(offset + span / 2); ticklab.append(c)
    offset += span + GAP

axA.set_xticks(ticks); axA.set_xticklabels(ticklab, fontsize=6)
axA.tick_params(axis="x", length=0)
axA.set_xlim(-GAP, offset)
axA.set_ylim(0, 0.37)
axA.set_ylabel("recombination rate\n(cM/Mb, sex-averaged)")
axA.set_title("a   Autosomal recombination landscape", loc="left")
axA.set_xlabel("chromosome")
despine(axA)

# panel B: chrZ male vs female
zm, zf = _bscorr("run_Z_male"), _bscorr("run_Z_female")
for d, col, lab in ((zm, SEXES["male"], "male (ZZ)"), (zf, SEXES["female"], "female (ZW)")):
    if d is not None:
        x, y = rollmed(d["mid"].to_numpy(), d["cM_Mb"].to_numpy())
        axB.plot(x, y, color=col, lw=1.6, label=lab)
        axB.text(x[np.isfinite(y)][-1], y[np.isfinite(y)][-1], "  " + lab.split()[0],
                 color=col, fontsize=8, va="center")
axB.set_title("b   chrZ recombination (sex-biased)", loc="left")
axB.set_xlabel("chrZ position (Mb)")
axB.set_ylabel("cM/Mb")
axB.set_ylim(0, 0.38)
despine(axB)

fig.savefig(OUT)
print("wrote", OUT, "| autosomes:", ",".join(chroms), "| chr2 included:", "2" in chroms)
