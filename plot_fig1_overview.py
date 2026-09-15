import os
#!/usr/bin/env python
"""Fig 1: no autosomal structure (excl chr2 inversion) + ZO sex assignment (no W chromosome; female is Z-hemizygous).
Run with bioinfo-buddy python.
"""
import sys; sys.path.insert(0, "/sietch_colab/data_share/illex/popgen_data/analysis/manuscript")
_HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(_HERE, "data")
from figstyle import apply, C, SEXES, despine
apply()
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

SC = DATA
SEXDIR = "/sietch_colab/data_share/illex/popgen_data/analysis/steps/12_sex"
OUT = "/sietch_colab/data_share/illex/popgen_data/analysis/manuscript/figures/fig1_overview.png"

AUTO_PCA_PREFIX = "auto_pca4"   # set by driver once final PCA is ready
AUTO_N_SNPS = None              # filled below from log if present
Z_N_SNPS = 6199

# --- load sex assignment (chrZ:autosome coverage ratio -> GMM -> sex) ---
sexdf = pd.read_csv(f"{SEXDIR}/sex_assignment.tsv", sep="\t")
sex_map = dict(zip(sexdf["sample"], sexdf["sex"]))
n_male = (sexdf.sex == "male").sum()
n_female = (sexdf.sex == "female").sum()
# chr42-corrected Z:autosome coverage ratio: chr42 is a known-autosomal chromosome
# that shares chrZ's small-chromosome mapping/composition bias (elevated read density
# in ALL samples regardless of sex, ~1.27-1.30x); dividing it out isolates the
# sex-linked ploidy signal from that shared technical bias.
sexdf["Z_auto_corrected"] = sexdf["Z_ratio"] / sexdf["chr42_ratio"]

# --- panel a: autosomal PCA excluding chr2 ---
auto = pd.read_csv(f"{SC}/{AUTO_PCA_PREFIX}.eigenvec", sep="\t")
auto_eig = np.loadtxt(f"{SC}/{AUTO_PCA_PREFIX}.eigenval")
auto_pve = 100 * auto_eig / auto_eig.sum()
with open(f"{SC}/auto_nsnps.txt") as fh:
    auto_nsnps = int(fh.read().strip())

# --- panel b: chrZ PCA ---
z = pd.read_csv(f"{SC}/z_pca.eigenvec", sep="\t")
z_eig = np.loadtxt(f"{SC}/z_pca.eigenval")
z_pve = 100 * z_eig / z_eig.sum()
# coverage-based sex assignment (chrZ:chr1 read ratio; NOT the prior labels), with
# an unknown band in the valley -- shared by panels b and c.
_cov = pd.read_csv("/sietch_colab/data_share/illex/popgen_data/analysis/steps/12_sex/chrZ_chr1_cov.tsv", sep="\t")
LOTHR, HITHR = 0.97, 1.12
def _covsex(r):
    if not np.isfinite(r) or r <= 0.3: return "unknown"
    return "ZO" if r < LOTHR else ("ZZ" if r > HITHR else "unknown")
cov_sex = {s: _covsex(r) for s, r in zip(_cov["sample"], _cov["Z_over_1"])}
z["covsex"] = z["IID"].map(cov_sex).fillna("unknown")

fig = plt.figure(figsize=(11.5, 4.0))
gs = fig.add_gridspec(1, 3, wspace=0.42)
axA = fig.add_subplot(gs[0, 0])
axB = fig.add_subplot(gs[0, 1])
axC = fig.add_subplot(gs[0, 2])

# panel a
# axis clipped to the robust core (median +/- 8*MAD on PC2) to show the single
# undifferentiated cloud; all n=621 samples are plotted, a handful of
# high-leverage individuals (elevated relatedness/quality, ~1.3%) fall outside
# the view -- see report for identities.
axA.scatter(auto["PC1"], auto["PC2"], s=10, color=C["ink"], alpha=0.3, edgecolors="none")
med2 = auto["PC2"].median(); mad2 = (auto["PC2"] - med2).abs().median() * 1.4826
axA.set_ylim(med2 - 8 * mad2, med2 + 8 * mad2)
axA.set_xlabel(f"PC1 ({auto_pve[0]:.1f}%)")
axA.set_ylabel(f"PC2 ({auto_pve[1]:.1f}%)")
axA.set_title(f"autosomal PCA, chr2 excluded\n(n={len(auto)}, {auto_nsnps:,} SNPs)", fontsize=9, fontweight="normal")
despine(axA)

# panel b -- coloured by the coverage-based assignment (ZO / ZZ / unknown)
COVCOL = {"ZO": SEXES["female"], "ZZ": SEXES["male"], "unknown": "#9E9E9E"}
COVLAB = {"ZO": "ZO female", "ZZ": "ZZ male", "unknown": "unknown"}
for grp in ["ZZ", "ZO", "unknown"]:
    m = z["covsex"] == grp
    if m.sum():
        axB.scatter(z.loc[m, "PC1"], z.loc[m, "PC2"], s=10, color=COVCOL[grp],
                    alpha=(0.5 if grp != "unknown" else 0.9),
                    edgecolors="none", zorder=(3 if grp == "unknown" else 2),
                    label=f"{COVLAB[grp]} (n={m.sum()})")
axB.set_xlabel(f"PC1 ({z_pve[0]:.1f}%)")
axB.set_ylabel(f"PC2 ({z_pve[1]:.1f}%)")
axB.set_title(f"chrZ PCA\n(n={len(z)}, {Z_N_SNPS:,} SNPs)", fontsize=9, fontweight="normal")
axB.legend(loc="best", markerscale=1.6, handletextpad=0.3)
despine(axB)

# panel c: DIRECT chrZ:chr1 read-coverage ratio per sample, assigned BY COVERAGE
# ALONE (not from any prior label -- avoids circularity). chrZ and chr1 read densities
# from samtools idxstats. Two clean modes -> ZO (1 Z) / ZZ (2 Z); the valley between
# them is left UNKNOWN.
cov = pd.read_csv("/sietch_colab/data_share/illex/popgen_data/analysis/steps/12_sex/chrZ_chr1_cov.tsv", sep="\t")
r = cov["Z_over_1"].to_numpy(float)
r = r[np.isfinite(r) & (r > 0.3)]                # drop zero-coverage failed samples
LOTHR, HITHR = 0.97, 1.12                        # valley edges (data-driven)
zo, unk, zz = r[r < LOTHR], r[(r >= LOTHR) & (r <= HITHR)], r[r > HITHR]
bins = np.linspace(0.6, 1.5, 46)
axC.hist(zo, bins=bins, color=SEXES["female"], alpha=0.85, label=f"ZO, 1 Z  (n={len(zo)})")
axC.hist(zz, bins=bins, color=SEXES["male"], alpha=0.85, label=f"ZZ, 2 Z  (n={len(zz)})")
axC.hist(unk, bins=bins, color="#9E9E9E", alpha=0.85, label=f"unknown  (n={len(unk)})")
for t in (LOTHR, HITHR):
    axC.axvline(t, ls="--", lw=0.8, color=C["muted"], zorder=0)
axC.set_xlabel("chrZ : chr1 read-coverage ratio\n(direct, per sample)")
axC.set_ylabel("number of samples")
axC.set_title("sex from coverage alone", fontsize=9, fontweight="normal")
axC.legend(fontsize=7.0, loc="upper center", frameon=False)
mzo, mzz = np.median(zo), np.median(zz)
axC.text(0.5, 0.52, f"modes {mzo:.2f} / {mzz:.2f}\nratio {mzz/mzo:.2f} ($<$2:1: chrZ\nover-covered $\\sim$1.26$\\times$\nin all samples)",
         transform=axC.transAxes, fontsize=6.6, va="top", ha="center", color=C["muted"])
despine(axC)

for ax, letter in zip([axA, axB, axC], "abc"):
    ax.text(-0.06, 1.14, letter, transform=ax.transAxes, fontsize=13, fontweight="bold", va="top", ha="left")

fig.savefig(OUT, dpi=200, bbox_inches="tight")
print("wrote", OUT)
med_m = sexdf.loc[sexdf.sex=="male","Z_auto_corrected"].median()
med_f = sexdf.loc[sexdf.sex=="female","Z_auto_corrected"].median()
print(f"panel a: autosomal PCA excl chr2, n_samples={len(auto)}, n_snps={auto_nsnps}")
print(f"panel b: chrZ PCA, n_samples={len(z)}, n_snps={Z_N_SNPS}; male={(z.sex=='male').sum()} female={(z.sex=='female').sum()}")
print(f"panel c: median chr42-corrected Z:auto ratio -- male={med_m:.3f} female={med_f:.3f} (n male={n_male} female={n_female})")
