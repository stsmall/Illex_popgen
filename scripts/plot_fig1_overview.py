import os
#!/usr/bin/env python
"""Fig 1: no autosomal structure (excl chr2 inversion) + ZO sex assignment (no W chromosome; female is Z-hemizygous).
Run with bioinfo-buddy python.
"""
import sys; sys.path.insert(0, "/sietch_colab/data_share/illex/popgen_data/analysis/manuscript")
_HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(_HERE, "..", "data")
from figstyle import apply, C, SEXES, despine
apply()
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

SC = DATA
SEXDIR = "/sietch_colab/data_share/illex/popgen_data/analysis/steps/12_sex"
OUT = os.path.join(_HERE, "..", "figures", "fig2_overview.png")

AUTO_PCA_PREFIX = "auto_pca4"   # set by driver once final PCA is ready
AUTO_N_SNPS = None              # filled below from log if present
Z_N_SNPS = 6199

# 12 BAMs with essentially no data (data/qc/empty_bams_12.tsv) are excluded from every panel (n=621).
EMPTY = set(pd.read_csv(os.path.join(DATA, "qc", "empty_bams_12.tsv"), sep="\t", comment="#")["sample"])

# --- load sex assignment (chrZ:autosome coverage ratio -> 2-component GMM -> sex; assign_sex.py) ---
# 'confident' = GMM posterior >= 0.9. Panels b and c both use these calls.
sexdf = pd.read_csv(f"{SEXDIR}/sex_assignment.tsv", sep="\t")
sexdf = sexdf[~sexdf["sample"].isin(EMPTY)].reset_index(drop=True)
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
auto = auto[~auto["IID"].isin(EMPTY)]
# auto_pca4.eigenval and auto_nsnps.txt are no longer on disk (searched analysis/, illex_manuscript/,
# home and scratch, 2026-11-11). Panel a is unchanged (same 621-sample eigenvec), so when they are
# missing the values are taken from the previous render of this panel (figures/fig2_overview.png of
# 2026-09-16, made from those files): PC1 29.2%, PC2 10.0%, 170,907 SNPs.
if os.path.exists(f"{SC}/{AUTO_PCA_PREFIX}.eigenval") and os.path.exists(f"{SC}/auto_nsnps.txt"):
    auto_eig = np.loadtxt(f"{SC}/{AUTO_PCA_PREFIX}.eigenval")
    auto_pve = 100 * auto_eig / auto_eig.sum()
    with open(f"{SC}/auto_nsnps.txt") as fh:
        auto_nsnps = int(fh.read().strip())
else:
    print("WARNING: auto_pca4.eigenval / auto_nsnps.txt missing -> panel a labels from the 2026-09-16 render")
    auto_pve = np.array([29.2, 10.0])
    auto_nsnps = 170907

# --- panel b: chrZ PCA ---
z = pd.read_csv(f"{SC}/z_pca.eigenvec", sep="\t")
z_eig = np.loadtxt(f"{SC}/z_pca.eigenval")
z_pve = 100 * z_eig / z_eig.sum()
z = z[~z["IID"].isin(EMPTY)]
# mixture-model sex calls (sex_assignment.tsv) shared by panels b and c
z = z.merge(sexdf[["sample", "sex", "confident"]], left_on="IID", right_on="sample", how="left")
assert z["sex"].notna().all()

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

# panel b -- coloured by the mixture-model sex call; filled = confident (posterior >= 0.9),
# open = not confident
SEXLAB = {"male": "ZZ male", "female": "ZO female"}
for sx in ["male", "female"]:
    ms = z["sex"] == sx
    mc = ms & z["confident"]; ml = ms & ~z["confident"]
    axB.scatter(z.loc[mc, "PC1"], z.loc[mc, "PC2"], s=10, color=SEXES[sx], alpha=0.5,
                edgecolors="none", zorder=2,
                label=f"{SEXLAB[sx]} (n={ms.sum()}, {mc.sum()} confident)")
    if ml.sum():
        axB.scatter(z.loc[ml, "PC1"], z.loc[ml, "PC2"], s=16, facecolors="none",
                    edgecolors=SEXES[sx], linewidths=0.9, zorder=3,
                    label=f"{SEXLAB[sx]}, not confident (n={ml.sum()})")
axB.set_xlabel(f"PC1 ({z_pve[0]:.1f}%)")
axB.set_ylabel(f"PC2 ({z_pve[1]:.1f}%)")
axB.set_title(f"chrZ PCA\n(n={len(z)}, {Z_N_SNPS:,} SNPs)", fontsize=9, fontweight="normal")
axB.legend(loc="best", markerscale=1.6, handletextpad=0.3, fontsize=7.0)
despine(axB)

# panel c: chrZ:autosome read-density ratio per sample (samtools idxstats; sex_from_cov.sh), the
# quantity the two-component Gaussian mixture of assign_sex.py is fitted to. Bars are coloured by the
# mixture-model call; hatched = not confident (posterior < 0.9).
r_all = sexdf["Z_ratio"].to_numpy(float)
lo = np.floor(np.nanmin(r_all) * 20) / 20; hi = np.ceil(np.nanmax(r_all) * 20) / 20
bins = np.arange(lo, hi + 1e-9, 0.02)
for sx, zlab in [("female", "ZO, 1 Z"), ("male", "ZZ, 2 Z")]:
    ms = sexdf["sex"] == sx
    rc = sexdf.loc[ms & sexdf["confident"], "Z_ratio"]; rl = sexdf.loc[ms & ~sexdf["confident"], "Z_ratio"]
    axC.hist(rc, bins=bins, color=SEXES[sx], alpha=0.85,
             label=f"{zlab}  (n={ms.sum()}, {len(rc)} confident)")
    if len(rl):
        axC.hist(rl, bins=bins, facecolor="white", edgecolor=SEXES[sx], hatch="////", lw=0.6,
                 label=f"{zlab}, not confident  (n={len(rl)})")
axC.set_xlabel("chrZ : autosome read-coverage ratio\n(per sample)")
axC.set_ylabel("number of samples")
axC.set_title("sex from coverage (Gaussian mixture)", fontsize=9, fontweight="normal")
axC.legend(fontsize=6.6, loc="upper center", frameon=False)
mzo = sexdf.loc[(sexdf.sex == "female") & sexdf.confident, "Z_ratio"].median()
mzz = sexdf.loc[(sexdf.sex == "male") & sexdf.confident, "Z_ratio"].median()
axC.text(0.84, 0.6, f"medians {mzo:.2f} / {mzz:.2f}\nratio {mzz/mzo:.2f} ($<$2:1)",
         transform=axC.transAxes, fontsize=6.6, va="top", ha="center", color=C["muted"])
despine(axC)

for ax, letter in zip([axA, axB, axC], "abc"):
    ax.text(-0.06, 1.14, letter, transform=ax.transAxes, fontsize=13, fontweight="bold", va="top", ha="left")

fig.savefig(OUT, dpi=200, bbox_inches="tight")
print("wrote", OUT)
med_m = sexdf.loc[sexdf.sex=="male","Z_auto_corrected"].median()
med_f = sexdf.loc[sexdf.sex=="female","Z_auto_corrected"].median()
print(f"panel a: autosomal PCA excl chr2, n_samples={len(auto)}, n_snps={auto_nsnps}")
print(f"panel b: chrZ PCA, n_samples={len(z)}, n_snps={Z_N_SNPS}; male={(z.sex=='male').sum()} ({(z.confident & (z.sex=='male')).sum()} confident) female={(z.sex=='female').sum()} ({(z.confident & (z.sex=='female')).sum()} confident)")
print(f"panel c: n={len(sexdf)}; male={n_male} ({(sexdf.confident & (sexdf.sex=='male')).sum()} confident) female={n_female} ({(sexdf.confident & (sexdf.sex=='female')).sum()} confident); medians confident ZO={mzo:.3f} ZZ={mzz:.3f}")
print(f"panel c: median chr42-corrected Z:auto ratio -- male={med_m:.3f} female={med_f:.3f} (n male={n_male} female={n_female})")
