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

# PCAs from pipeline/steps/05_karyo_pca_arms/build_pca621.sh (621 individuals, LD-thinned SNPs, MAF >= 0.05).
# Variance explained = eigenvalue / trace of the full relationship matrix (in *_pca621_summary.tsv).

# 12 BAMs with essentially no data (data/qc/empty_bams_12.tsv) are excluded from every panel (n=621).
EMPTY = set(pd.read_csv(os.path.join(DATA, "qc", "empty_bams_12.tsv"), sep="\t", comment="#")["sample"])

# --- load sex assignment (chrZ:autosome coverage ratio -> 2-component GMM -> sex; assign_sex.py) ---
# 'confident' = GMM posterior >= 0.9. Panels b and c both use these calls.
sexdf = pd.read_csv(f"{SEXDIR}/sex_assignment.tsv", sep="\t")
sexdf = sexdf[~sexdf["sample"].isin(EMPTY)].reset_index(drop=True)
# manual overrides (data/qc/sex_overrides.tsv): five calls with male-like chrZ heterozygosity are set to not
# confident, and the two with chrZ coverage above the male range are recalled as male
_ov = pd.read_csv(os.path.join(DATA, "qc", "sex_overrides.tsv"), sep="\t", comment="#").set_index("sample")
_i = sexdf["sample"].isin(_ov.index)
sexdf.loc[_i, "sex"] = sexdf.loc[_i, "sample"].map(_ov["sex"])
sexdf.loc[_i, "confident"] = sexdf.loc[_i, "sample"].map(_ov["confident"]).astype(bool)
sex_map = dict(zip(sexdf["sample"], sexdf["sex"]))
n_male = (sexdf.sex == "male").sum()
n_female = (sexdf.sex == "female").sum()
# chr42-corrected Z:autosome coverage ratio: chr42 is a known-autosomal chromosome
# that shares chrZ's small-chromosome mapping/composition bias (elevated read density
# in ALL samples regardless of sex, ~1.27-1.30x); dividing it out isolates the
# sex-linked ploidy signal from that shared technical bias.
sexdf["Z_auto_corrected"] = sexdf["Z_ratio"] / sexdf["chr42_ratio"]

# --- panel a: autosomal PCA excluding chr2 ---
def load_pca(prefix):
    ev = pd.read_csv(f"{SC}/{prefix}.eigenvec", sep="\t")
    sm = pd.read_csv(f"{SC}/{prefix}_summary.tsv", sep="\t").iloc[0]
    assert not ev["IID"].isin(EMPTY).any()
    return ev, np.array([sm.pve_pc1, sm.pve_pc2]), int(sm.n_snps)


auto, auto_pve, auto_nsnps = load_pca("auto_pca621")

# --- panel b: chrZ PCA ---
z, z_pve, z_nsnps = load_pca("z_pca621")
# mixture-model sex calls (sex_assignment.tsv) shared by panels b and c
z = z.merge(sexdf[["sample", "sex", "confident"]], left_on="IID", right_on="sample", how="left")
assert z["sex"].notna().all()

fig = plt.figure(figsize=(11.5, 4.6))
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
# open = not confident; the five overridden calls (data/qc/sex_overrides.tsv) in grey
SEXLAB = {"male": "ZZ male", "female": "ZO female"}
flag = z["sample"].isin(_ov.index)
for sx in ["male", "female"]:
    ms = (z["sex"] == sx) & ~flag
    mc = ms & z["confident"]; ml = ms & ~z["confident"]
    axB.scatter(z.loc[mc, "PC1"], z.loc[mc, "PC2"], s=10, color=SEXES[sx], alpha=0.5,
                edgecolors="none", zorder=2,
                label=f"{SEXLAB[sx]} (n={ms.sum()}, {mc.sum()} confident)")
    if ml.sum():
        axB.scatter(z.loc[ml, "PC1"], z.loc[ml, "PC2"], s=16, facecolors="none",
                    edgecolors=SEXES[sx], linewidths=0.9, zorder=3,
                    label=f"{SEXLAB[sx]}, not confident (n={ml.sum()})")
axB.scatter(z.loc[flag, "PC1"], z.loc[flag, "PC2"], s=18, color="#7f7f7f", edgecolors="white", linewidths=0.5,
            zorder=4, label=f"ambiguous, male-like chrZ heterozygosity (n={flag.sum()})")
axB.set_xlabel(f"PC1 ({z_pve[0]:.1f}%)")
axB.set_ylabel(f"PC2 ({z_pve[1]:.1f}%)")
axB.set_title(f"chrZ PCA\n(n={len(z)}, {z_nsnps:,} SNPs)", fontsize=9, fontweight="normal")
axB.legend(loc="upper center", bbox_to_anchor=(0.5, -0.2), ncol=2, markerscale=1.4, handletextpad=0.3, columnspacing=0.8, fontsize=6.5, frameon=False)
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
print(f"panel b: chrZ PCA, n_samples={len(z)}, n_snps={z_nsnps}; male={(z.sex=='male').sum()} ({(z.confident & (z.sex=='male')).sum()} confident) female={(z.sex=='female').sum()} ({(z.confident & (z.sex=='female')).sum()} confident)")
print(f"panel c: n={len(sexdf)}; male={n_male} ({(sexdf.confident & (sexdf.sex=='male')).sum()} confident) female={n_female} ({(sexdf.confident & (sexdf.sex=='female')).sum()} confident); medians confident ZO={mzo:.3f} ZZ={mzz:.3f}")
print(f"panel c: median chr42-corrected Z:auto ratio -- male={med_m:.3f} female={med_f:.3f} (n male={n_male} female={n_female})")
