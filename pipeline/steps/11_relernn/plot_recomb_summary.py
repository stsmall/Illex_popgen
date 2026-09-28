"""Recombination figures — SUMMARY (not a landscape).

Main (fig_recomb_summary.png):
  a. OVERALL autosomal recombination rate as a single MEAN LINE vs physical
     position (Mb), averaged across all autosomes (chr2 + Z excluded), male+female
     combined, with a +/-1 SD variance ribbon.
  b. chrZ male vs female (sex-biased).
Supplement (supp_recomb_boxplots.png):
  per-chromosome box plot of window rates (chr2 colinear vs inverted), + Z male/female.
"""
import glob, os, sys
import numpy as np, pandas as pd
sys.path.insert(0, "/sietch_colab/data_share/illex/popgen_data/analysis/manuscript")
from figstyle import apply, C, SEXES, despine
apply()
import matplotlib.pyplot as plt

RD = "/sietch_colab/data_share/illex/popgen_data/analysis/steps/11_relernn"
FIG = "/sietch_colab/data_share/illex/popgen_data/analysis/manuscript/old_figs/fig_recomb_summary.png"
SUP = "/sietch_colab/data_share/illex/popgen_data/analysis/manuscript/supp_figures/supp_recomb_boxplots.png"


def load(path):
    if not path or not os.path.exists(path):
        return None
    d = pd.read_csv(path, sep="\t")
    d["chrom"] = d["chrom"].astype(str).str.replace(r"^b'|'$", "", regex=True)
    d["cM"] = d["recombRate"].astype(float) * 1e8
    d["mid"] = (d["start"] + d["end"]) / 2.0 / 1e6      # Mb
    return d[["chrom", "mid", "cM"]]


def bscorr(proj):
    h = glob.glob(os.path.join(RD, proj, "proj", "*.PREDICT.BSCORRECTED.txt"))
    return load(h[0]) if h else None


# pooled male+female autosomes (chr2 excluded from these runs already)
auto = pd.concat([d for d in (bscorr("run_male_auto"), bscorr("run_female_auto"))
                  if d is not None], ignore_index=True)
chr2 = pd.concat([load(os.path.join(RD, "chr2_autosomal_predict", f)) for f in
                  ("chr2_male.autonet.PREDICT.txt", "chr2_female.autonet.PREDICT.txt")],
                 ignore_index=True)
zm, zf = bscorr("run_Z_male"), bscorr("run_Z_female")
autos = sorted(auto["chrom"].unique(), key=int)

# ---------------- MAIN: mean line vs position, with variance ribbon --------
BIN = 2.0  # Mb
auto["bin"] = (auto["mid"] // BIN) * BIN + BIN / 2.0
# only bins with enough autosomes contributing (>=5) so the far tail isn't noisy
g = auto.groupby("bin")["cM"]
prof = pd.DataFrame({"mean": g.mean(), "sd": g.std(), "n": g.count()}).reset_index()
prof = prof[prof["n"] >= 5]
gmean = auto["cM"].mean()

fig, (axa, axb) = plt.subplots(1, 2, figsize=(10.5, 4.0),
                               gridspec_kw={"width_ratios": [3.6, 1]})
axa.fill_between(prof["bin"], prof["mean"] - prof["sd"], prof["mean"] + prof["sd"],
                 color=C["accent"], alpha=0.18, lw=0, label="$\\pm$1 SD across autosomes")
axa.plot(prof["bin"], prof["mean"], color=C["accent"], lw=1.8, label="mean autosomal rate")
axa.axhline(gmean, color=C["muted"], lw=0.8, ls="--")
axa.text(prof["bin"].max(), gmean, f"  genome mean {gmean:.3f}", va="center",
         fontsize=8, color=C["muted"])
axa.set_xlim(0, prof["bin"].max())
axa.set_ylim(0, 0.36)
axa.set_xlabel("position along chromosome (Mb)")
axa.set_ylabel("recombination rate (cM/Mb)")
axa.set_title("a   Autosomal recombination rate (chr2 and Z excluded; male + female combined)",
              loc="left", fontsize=10.5)
axa.legend(loc="lower left", fontsize=8)
despine(axa)

# chrZ male vs female as points +/- SD
for i, (lab, d, col) in enumerate([("male", zm, SEXES["male"]), ("female", zf, SEXES["female"])]):
    if d is not None:
        v = d["cM"].to_numpy()
        axb.errorbar([i], [v.mean()], yerr=[v.std()], fmt="o", ms=7, color=col,
                     ecolor=col, elinewidth=1.2, capsize=4)
        axb.text(i, v.mean() + v.std() + 0.012, lab, ha="center", fontsize=8.5, color=col)
axb.set_xlim(-0.6, 1.6); axb.set_xticks([]); axb.set_ylim(0, 0.38)
axb.set_ylabel("cM/Mb")
axb.set_title("b   chrZ", loc="left")
despine(axb)
fig.savefig(FIG, dpi=220, bbox_inches="tight")
print("wrote", FIG, "| autosomal mean %.3f cM/Mb; %d Mb-bins" % (gmean, len(prof)))

# ---------------- SUPPLEMENT: per-chrom boxplots ----------------
per = {c: auto[auto["chrom"] == c]["cM"].to_numpy() for c in autos}
# split chr2 into colinear flanks vs the inverted region (mid in Mb; breakpoints 60.54-79.50)
BP_L2, BP_R2 = 60.54, 79.50
chr2_inv = chr2[(chr2["mid"] >= BP_L2) & (chr2["mid"] <= BP_R2)]["cM"].to_numpy()
chr2_col = chr2[(chr2["mid"] < BP_L2) | (chr2["mid"] > BP_R2)]["cM"].to_numpy()
data = [per[c] for c in autos] + [chr2_col, chr2_inv]
labels = autos + ["2 col", "2 inv"]
figs, axs = plt.subplots(figsize=(9.6, 5.0))
bp = axs.boxplot(data, positions=np.arange(len(data)), widths=0.72, showfliers=False,
                 patch_artist=True, medianprops=dict(color=C["warm"], lw=1.1))
for i, patch in enumerate(bp["boxes"]):
    fc = C["faint"]
    if labels[i] == "2 col":
        fc = "#f6d5c4"          # chr2 colinear flanks (light)
    elif labels[i] == "2 inv":
        fc = "#D55E00"          # chr2 inverted region (suppressed recomb) flagged
    patch.set(facecolor=fc, edgecolor=C["muted"], lw=0.7)
base = len(data)
for j, (lab, d, col) in enumerate([("Z♂", zm, SEXES["male"]), ("Z♀", zf, SEXES["female"])]):
    if d is not None:
        b = axs.boxplot([d["cM"].to_numpy()], positions=[base + j], widths=0.72,
                        showfliers=False, patch_artist=True,
                        medianprops=dict(color="black", lw=1.1))
        b["boxes"][0].set(facecolor=col, alpha=0.5, edgecolor=col, lw=0.8)
        labels.append(lab)
# autosomal expectation line (median over the other autosomes) -- chr2 colinear AND
# inverted both sit well below it (chromosome-wide suppression; Mann-Whitney p<1e-300,
# not a chromosome-size effect: no length-rate relationship across autosomes)
auto_med = float(np.median(np.concatenate([per[c] for c in autos])))
axs.axhline(auto_med, ls="--", lw=1.0, color=C["ink"], zorder=0)
axs.text(len(autos) - 0.5, auto_med + 0.006, f"autosomal median {auto_med:.2f}",
         fontsize=7.5, color=C["ink"], ha="right", va="bottom")
axs.set_xlim(-0.7, len(labels) - 0.3)
axs.set_xticks(np.arange(len(labels)))
axs.set_xticklabels(labels, fontsize=6.5, rotation=90)
axs.tick_params(axis="x", length=0)
axs.set_ylim(0, 0.42)
axs.set_ylabel("recombination rate (cM/Mb)")
axs.set_xlabel("chromosome  (chr2 split: colinear flanks vs inverted region; Z shown male/female)")
axs.set_title("Per-chromosome recombination-rate distribution", loc="left")
despine(axs)
figs.savefig(SUP, dpi=220, bbox_inches="tight")
print("wrote", SUP)
