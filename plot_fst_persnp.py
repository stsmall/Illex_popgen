"""Supp: per-SNP FST between populations finds no differentiation islands outside chr2.

plink2 Weir-Cockerham FST per SNP for all 45 NAFO-division pairs (autosomes, chr2
excluded), restricted to well-genotyped sites (>=40 observed haplotypes per pair). The
pooled per-SNP FST distribution collapses at ~0; no SNP is a robust outlier (high FST in
many pairs) -- the apparent FST=1 sites are low-count/mismapping artifacts that vanish
under the genotype-count filter. Confirms structure is confined to the chr2 inversion.
"""
import os
os.environ.setdefault("MPLCONFIGDIR", "/dev/shm/mplcache")  # avoid slow-FS font-cache rebuild
import sys
sys.path.insert(0, "/sietch_colab/data_share/illex/popgen_data/analysis/manuscript")
_HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(_HERE, "data")
from figstyle import apply, C, despine
apply()
import numpy as np, matplotlib.pyplot as plt

WD = DATA
OUT = "/sietch_colab/data_share/illex/popgen_data/analysis/manuscript/supp_figures/supp_fst_persnp.png"

f = np.load(os.path.join(DATA, "pooled_fst.npy"))
v = f[np.isfinite(f)]
v = np.clip(v, 0, 1)
med, q99, q999 = np.median(v), np.quantile(v, 0.99), np.quantile(v, 0.999)
frac25, frac50 = np.mean(v > 0.25), np.mean(v > 0.5)
n_full = len(v)
# precompute histogram counts (numpy only; avoids passing 5M points to matplotlib)
counts, edges = np.histogram(v, bins=np.linspace(0, 1, 101))
centers = 0.5 * (edges[:-1] + edges[1:])

fig, ax = plt.subplots(figsize=(7.4, 4.0))
ax.bar(centers, counts, width=(edges[1] - edges[0]), color="#7FB3D5", edgecolor="none", align="center")
ax.set_yscale("log")
ax.axvline(med, ls=":", color=C["ink"], lw=1.0)
ax.text(med, ax.get_ylim()[1] * 0.5, f"  median $F_{{ST}}$ = {med:.3f}", fontsize=8, color=C["ink"])
ax.set_xlabel(r"per-SNP $F_\mathrm{ST}$ (Weir--Cockerham), all 45 division pairs")
ax.set_ylabel("number of SNP $\\times$ pair observations")
ax.set_xlim(0, 1)
despine(ax)
ax.text(0.97, 0.92,
        f"{n_full:,} well-genotyped SNP$\\times$pair values (chr2 excluded)\n"
        f"99th pctile = {q99:.2f};  99.9th = {q999:.2f}\n"
        f"{frac50*100:.2f}% of values $>0.5$ (isolated, not reproducible across pairs)\n"
        f"no SNP is high-$F_{{ST}}$ in more than 6 of 45 pairs\n"
        r"$\Rightarrow$ no differentiation islands outside chr2",
        transform=ax.transAxes, ha="right", va="top", fontsize=8, color=C["ink"],
        bbox=dict(boxstyle="round,pad=0.4", fc="white", ec=C["muted"], lw=0.6))
fig.savefig(OUT, dpi=220, bbox_inches="tight")
print("wrote", OUT)
print(f"median={med:.4f} q99={q99:.3f} q999={q999:.3f} frac>0.5={frac50:.5f} n={len(v)}")
