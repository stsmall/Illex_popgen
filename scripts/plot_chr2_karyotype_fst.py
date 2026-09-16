"""Supp: per-SNP F_ST between the chr2 arrangements (AA vs BB karyotypes) across the
WHOLE of chromosome 2, aggregated into windows (ratio-of-averages Hudson estimator).
Unlike the population F_ST scan (which asks about geography), this asks how differentiated
the two arrangements are along the chromosome: collinear chr2 sits at ~0 (the karyotypes
are one panmictic population outside the inversion), while the inverted body (60.54-79.50 Mb)
is strongly differentiated -- the F_ST~0.37 'structure' is between karyotypes and confined
to the inversion. AA n=254, BB n=95 arrangement homozygotes."""
import os
os.environ.setdefault("MPLCONFIGDIR", "/dev/shm/mplcache")
import sys
_HERE = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.dirname(_HERE) if os.path.basename(_HERE) == "scripts" else _HERE
sys.path.insert(0, _HERE)
from figstyle import apply, C, despine
apply()
import numpy as np, pandas as pd, matplotlib.pyplot as plt

O = "/sietch_colab/data_share/illex/popgen_data/analysis/steps/03_karyotype/inversion_content"
OUT = os.path.join(_ROOT, "supp_figures", "figS20_chr2_karyotype_fst.png")
INV_S, INV_E = 60_540_000, 79_500_000
N1, N2 = 2 * 254, 2 * 95          # allele counts (diploid homozygotes)
WIN = 100_000

# Windowed Hudson FST precomputed by a streaming awk pass (ratio-of-averages, N1=508
# AA + N2=190 BB alleles): columns win_start, sum_num, sum_den, n_snps.
w = pd.read_csv(f"{O}/chr2_karyo_fst_windows.tsv", sep="\t", header=None,
                names=["win", "num", "den", "n"])
w = w[(w.n >= 10) & (w.den > 0)]
w["fst"] = (w.num / w.den).clip(lower=0)
w["mid"] = w.win + WIN / 2
w["inv"] = (w.mid >= INV_S) & (w.mid <= INV_E)

coll_med = float(w.loc[~w.inv, "fst"].median())
inv_med = float(w.loc[w.inv, "fst"].median())

fig, ax = plt.subplots(figsize=(11, 3.3))
ax.axvspan(INV_S / 1e6, INV_E / 1e6, color=C["shade"], lw=0, zorder=0)
ax.scatter(w.loc[~w.inv, "mid"] / 1e6, w.loc[~w.inv, "fst"], s=5, color=C["muted"],
           linewidths=0, alpha=0.7, zorder=2, rasterized=True, label=f"collinear (median {coll_med:.3f})")
ax.scatter(w.loc[w.inv, "mid"] / 1e6, w.loc[w.inv, "fst"], s=7, color=C["warm"],
           linewidths=0, alpha=0.85, zorder=3, label=f"inverted body (median {inv_med:.2f})")
ax.axhline(coll_med, color=C["ink"], ls=":", lw=0.9, zorder=1)
ax.set_xlabel("chromosome 2 position (Mb)")
ax.set_ylabel(r"$F_{ST}$ between arrangements  (AA vs BB)")
ax.set_ylim(-0.02, min(1.0, w.fst.max() * 1.08))
ax.set_xlim(0, w.mid.max() / 1e6)
ax.set_title("Between-karyotype differentiation is confined to the chromosome-2 inversion",
             loc="left", fontweight="bold", fontsize=11, pad=8)
ax.legend(loc="upper left", frameon=False, fontsize=9, handletextpad=0.4)
despine(ax)
fig.savefig(OUT, dpi=200, bbox_inches="tight")
print("wrote", OUT)
print(f"windows: {len(w)} | collinear median FST={coll_med:.4f} | inverted median FST={inv_med:.4f} "
      f"| inverted max={w.loc[w.inv,'fst'].max():.3f}")
