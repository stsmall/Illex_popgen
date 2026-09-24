"""Supp Fig S12: per-chromosome ReLERNN recombination-rate distributions on a CONSISTENT footing.
All boxes are the RAW (uncorrected) PREDICT output of the male autosomal network: 43 scanned autosomes,
chromosome 2 restricted to collinear sequence (predicted with the same male network), and
the Z (male and female raw predictions from the Z-specific runs). The previous version pooled
bootstrap-corrected male+female autosomal maps (the degenerate female map inflates the pooled median to 0.28)
against RAW chr2 predictions, which made chr2 look suppressed chromosome-wide. It is not: chr2 median 0.202
vs autosomal 0.200 (60th percentile of per-chromosome medians); inverted 0.2015 vs collinear 0.2023."""
import os, sys, glob
os.environ.setdefault("MPLCONFIGDIR", "/dev/shm/mplcache")
_HERE = os.path.dirname(os.path.abspath(__file__)); _ROOT = os.path.dirname(_HERE) if os.path.basename(_HERE) == "scripts" else _HERE
sys.path.insert(0, _HERE)
from figstyle import apply, C, SEXES, despine
apply()
import numpy as np, pandas as pd, matplotlib.pyplot as plt
RD = "/sietch_colab/data_share/illex/popgen_data/analysis/steps/11_relernn"
OUT = os.path.join(_ROOT, "supp_figures", "figS12_recomb_boxplots.png")
BP_L, BP_R = 60.54, 79.50

def load(path):
    d = pd.read_csv(path, sep="\t")
    d["chrom"] = d["chrom"].astype(str).str.replace(r"^b'|'$", "", regex=True)
    d["cM"] = d["recombRate"].astype(float) * 1e8
    d["mid"] = (d["start"] + d["end"]) / 2 / 1e6
    return d

auto = load(f"{RD}/run_male_auto/proj/male.kept.PREDICT.txt")
chr2 = load(f"{RD}/chr2_autosomal_predict/pred_male/chr2_male.PREDICT.txt")
def zraw(run):
    f = [p for p in glob.glob(f"{RD}/{run}/proj/*.PREDICT.txt") if "BSCORR" not in p]
    return load(f[0]) if f else None
zm, zf = zraw("run_Z_male"), zraw("run_Z_female")

autos = sorted(auto.chrom.unique(), key=int)
per = {c: auto.loc[auto.chrom == c, "cM"].to_numpy() for c in autos}
c2_inv = chr2[(chr2.mid >= BP_L) & (chr2.mid <= BP_R)].cM.to_numpy()
c2_col = chr2[(chr2.mid < BP_L) | (chr2.mid > BP_R)].cM.to_numpy()
data = [per[c] for c in autos] + [c2_col]
labels = autos + ["2*"]

fig, ax = plt.subplots(figsize=(13.5, 5.2))
bp = ax.boxplot(data, positions=np.arange(len(data)), widths=0.7, showfliers=False, patch_artist=True,
                medianprops=dict(color=C["warm"], lw=1.1), whiskerprops=dict(lw=0.7), capprops=dict(lw=0.7))
for i, patch in enumerate(bp["boxes"]):
    fc = "#f6d5c4" if labels[i] == "2*" else C["faint"]
    patch.set(facecolor=fc, edgecolor=C["muted"], lw=0.7)
pos = len(data)
for lab, d, col in (("Z", zm, C["accent"]),):
    if d is not None:
        b = ax.boxplot([d.cM.to_numpy()], positions=[pos], widths=0.7, showfliers=False, patch_artist=True,
                       medianprops=dict(color="black", lw=1.1))
        b["boxes"][0].set(facecolor=col, alpha=0.5, edgecolor=col, lw=0.8)
        labels.append(lab); pos += 1
auto_med = float(np.median(np.concatenate([per[c] for c in autos])))
ax.axhline(auto_med, ls="--", lw=1.0, color=C["ink"], zorder=0)
ax.text(1.0, 1.01, f"dashed line: autosomal median {auto_med:.3f} cM/Mb", transform=ax.transAxes, fontsize=7.5, color=C["muted"], va="bottom", ha="right")
pct = (pd.Series({c: np.median(per[c]) for c in autos}) < np.median(c2_col)).mean() * 100
ax.set_xticks(np.arange(len(labels))); ax.set_xticklabels(labels, fontsize=6.5, rotation=90)
ax.set_ylabel("recombination rate (cM/Mb)"); ax.set_ylim(0, None)
ax.set_xlabel("chromosome  (2* = chromosome 2 collinear sequence only; inversion excluded, see LD in Supplementary Fig. S4)")
ax.set_title("Per-chromosome recombination-rate distribution", loc="left", fontweight="bold")
despine(ax); fig.savefig(OUT, dpi=200, bbox_inches="tight")
print("wrote", OUT, f"| auto median {auto_med:.4f} chr2 collinear {np.median(c2_col):.4f} pct {pct:.0f}")