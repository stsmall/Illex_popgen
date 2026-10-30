"""Supp Fig S12: per-chromosome male ReLERNN recombination-rate distributions, raw and bootstrap-corrected.
Autosomes: the network trained on chromosomes 1, 13 and 45 (run_Achk_male_uf, 6 kb windows) applied to every
autosome except chromosome 2 (run_male_uf_achknet), corrected with that network's 200-replicate bootstrap table
(BSCORRECTED200, analysis/steps/11_relernn/bs200_correct_genomewide.py).
Z: the Z network trained under the autosomal history with Ne x 0.60 (run_Z_male_uf, 30 kb windows), so its rate
is already on the autosomal scale, corrected with its own 200-replicate table.
r_Z/r_A is the ratio of the Z median to the median of all autosomal windows."""
import os, sys
os.environ.setdefault("MPLCONFIGDIR", "/dev/shm/mplcache")
_HERE = os.path.dirname(os.path.abspath(__file__)); _ROOT = os.path.dirname(_HERE) if os.path.basename(_HERE) == "scripts" else _HERE
sys.path.insert(0, _HERE)
from figstyle import apply, C, despine
apply()
import numpy as np, pandas as pd, matplotlib.pyplot as plt
RD = "/sietch_colab/data_share/illex/popgen_data/analysis/steps/11_relernn"
OUT = os.path.join(_ROOT, "supp_figures", "figS12_recomb_boxplots.png")


def load(proj, bn, bsfile):
    p = f"{RD}/{proj}/proj/{bn}.PREDICT"
    d = pd.read_csv(p + ".txt", sep="\t")
    d["chrom"] = d["chrom"].astype(str).str.replace(r"^b'|'$", "", regex=True)
    d["raw"] = d["recombRate"].astype(float) * 1e8
    d["bs"] = pd.read_csv(p + f".{bsfile}.txt", sep="\t")["recombRate"].to_numpy() * 1e8
    return d


A = load("run_male_uf_achknet", "male_fh.kept", "BSCORRECTED200")
Z = load("run_Z_male_uf", "Z_male_uf.kept", "BSCORRECTED")
autos = sorted(A.chrom.unique(), key=int)
labels = autos + ["Z"]

fig, axs = plt.subplots(2, 1, figsize=(13.5, 8.4), sharex=True)
for ax, col, lab, tag in zip(axs, ["raw", "bs"], ["raw", "bootstrap-corrected"], "ab"):
    data = [A.loc[A.chrom == c, col].to_numpy() for c in autos] + [Z[col].to_numpy()]
    bp = ax.boxplot(data, positions=np.arange(len(data)), widths=0.7, showfliers=False, patch_artist=True,
                    medianprops=dict(color=C["warm"], lw=1.1), whiskerprops=dict(lw=0.7), capprops=dict(lw=0.7))
    for patch in bp["boxes"]:
        patch.set(facecolor=C["faint"], edgecolor=C["muted"], lw=0.7)
    bp["boxes"][-1].set(facecolor=C["accent"], alpha=0.5, edgecolor=C["accent"])
    amed, zmed = float(A[col].median()), float(Z[col].median())
    ax.axhline(amed, ls="--", lw=1.0, color=C["ink"], zorder=0)
    ax.text(1.0, 1.01, f"dashed line: autosomal median {amed:.3f} cM/Mb   Z median {zmed:.3f}   "
            f"$r_Z/r_A$ = {zmed / amed:.2f}", transform=ax.transAxes, fontsize=7.5, color=C["muted"],
            va="bottom", ha="right")
    ax.set_ylabel("recombination rate (cM/Mb)"); ax.set_ylim(0, None)
    ax.set_title(f"({tag}) male, {lab}", loc="left", fontweight="bold")
    despine(ax)
    print(f"{lab}: autosomal median {amed:.4f} Z median {zmed:.4f} ratio {zmed / amed:.3f}")
axs[-1].set_xticks(np.arange(len(labels))); axs[-1].set_xticklabels(labels, fontsize=6.5, rotation=90)
axs[-1].tick_params(axis="x", length=0)
axs[-1].set_xlabel("chromosome")
fig.tight_layout()
fig.savefig(OUT, dpi=200, bbox_inches="tight")
print("wrote", OUT)
