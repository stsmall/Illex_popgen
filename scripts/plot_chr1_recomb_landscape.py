"""Supp Fig S22: the chromosome-1 male ReLERNN recombination landscape (bootstrap-corrected, 200 replicates) with
95% bootstrap intervals, and what it does not show. Windows of 3-6 kb from the network trained on chromosomes 1,
13 and 45 (run_male_uf_achknet, BSCORRECTED200). Rolling medians over 100 windows (~0.6 Mb) and 333 windows
(~2 Mb). Shaded: the stretches with no windows, 24.0-31.3 Mb (a masked pi-elevated block) and 45.25-46.39 Mb
(inside the 80%-repeat block at 45.0-47.5 Mb, the centromere). The accessible pericentromeric flanks dip only a
few percent, so the map's chromosome-scale dynamic range is compressed."""
import os
os.environ.setdefault("MPLCONFIGDIR", "/dev/shm/mplcache")
import sys
_HERE = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.dirname(_HERE) if os.path.basename(_HERE) == "scripts" else _HERE
sys.path.insert(0, _HERE)
from figstyle import apply, C, despine
apply()
import numpy as np, pandas as pd, matplotlib.pyplot as plt

RD = "/sietch_colab/data_share/illex/popgen_data/analysis/steps/11_relernn"
OUT = os.path.join(_ROOT, "supp_figures", "figS22_chr1_recomb_landscape.png")
CHROM = "1"
CEN = (45.0, 47.5)          # 80%-repeat block (centromere)

d = pd.read_csv(f"{RD}/run_male_uf_achknet/proj/male_fh.kept.PREDICT.BSCORRECTED200.txt", sep="\t")
d["chrom"] = d["chrom"].astype(str).str.replace(r"^b'|'$", "", regex=True)
d = d[d.chrom == CHROM].copy().sort_values("start")
for c in ("recombRate", "CI95LO", "CI95HI"):
    d[c] = d[c].astype(float) * 1e8
d["mid"] = (d.start + d.end) / 2 / 1e6
gap = (d.start.values[1:] - d.end.values[:-1]) / 1e6
holes = [(d.end.values[i] / 1e6, d.start.values[i + 1] / 1e6) for i in np.where(gap > 1.0)[0]]
med = d.recombRate.median()
cs, ce = [h for h in holes if CEN[0] <= h[0] <= CEN[1]][0]
fl = d.recombRate[((d.mid > cs - 1) & (d.mid < cs)) | ((d.mid > ce) & (d.mid < ce + 1))].median() / med - 1

fig, ax = plt.subplots(1, 1, figsize=(13, 3.6))
ax.axvspan(*CEN, color=C["faint"], lw=0, zorder=0, alpha=0.6)
for h in holes:
    ax.axvspan(*h, color=C["faint"], lw=0, zorder=0)
col = C["accent"]
seg = np.cumsum(np.r_[0, (d.start.values[1:] - d.end.values[:-1]) > 1e6])
for k in np.unique(seg):
    s = d[seg == k]
    ax.fill_between(s.mid, s.CI95LO, s.CI95HI, color=col, alpha=0.22, lw=0,
                    label="95% bootstrap interval" if k == 0 else None)
    ax.plot(s.mid, s.recombRate, color=col, lw=0.4, alpha=0.9,
            label="ReLERNN estimate (3-6 kb windows)" if k == 0 else None)
    ax.plot(s.mid, s.recombRate.rolling(100, center=True, min_periods=40).median(), color=C["ink"], lw=1.4,
            label="rolling median (≈0.6 Mb)" if k == 0 else None)
broad = d.recombRate.rolling(333, center=True, min_periods=150).median()
ax.axhline(med, color=C["muted"], ls=":", lw=0.9)
ax.text(np.mean(holes[0]), 0.015, "masked\n(π-elevated\nblock)", ha="center", va="bottom", fontsize=7, color=C["ink"])
ax.text(np.mean(CEN), 0.015, "centro-\nmere", ha="center", va="bottom", fontsize=7, color=C["ink"])
ax.set_ylim(0, 0.8); ax.set_ylabel("recombination\n(cM/Mb)")
ax.set_title(f"Chromosome 1 recombination map  (median {med:.3f} cM/Mb)", loc="left", fontweight="bold", fontsize=10)
ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.2), frameon=False, fontsize=7.5, ncol=3); despine(ax)
ax.set_xlabel("chromosome 1 position (Mb)")
ax.set_xlim(0, d.end.max() / 1e6)
fig.savefig(OUT, dpi=200, bbox_inches="tight")
f6 = d.recombRate.rolling(100, center=True, min_periods=40).median()
print("wrote", OUT, f"| median {med:.3f}, 0.6Mb roll {f6.min():.3f}-{f6.max():.3f}, 2Mb roll {broad.min():.3f}-{broad.max():.3f}, "
      f"holes {[(round(a,2),round(b,2)) for a,b in holes]}, flank {fl:+.3f}")
