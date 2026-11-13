"""Among-karyotype variance (eta^2) of PC1 and the per-chromosome Gaussian-mixture screen, on PCAs fitted on the
619 analysed individuals (refit_pca_619.sh; final chr2 karyotypes from data/qc/karyotypes_619.tsv, karyotype_final).

eta^2 = between-karyotype sum of squares / total sum of squares of PC1 (one-way ANOVA). PC1 sign does not matter.
1) windowed PCA along chr2 (refit_619/chr2_winpca_619_coords.npy, 200-kb windows, 100-kb step): mean over windows
   inside the breakpoints (60.54-79.50 Mb) and over collinear windows >1 Mb from either breakpoint.
2) per-chromosome PCA of thinned SNPs (refit_619/karyo_pca_coords_619.tsv): eta^2 per chromosome, and, independent of
   the karyotype labels, a GMM screen of PC1: Delta BIC = BIC(k=1) - min(BIC(k=2), BIC(k=3)), cluster separation =
   (max - min of the 3-component means) / sqrt(mean component variance), n_init 4, random_state 0.
"""
import numpy as np
import pandas as pd
from sklearn.mixture import GaussianMixture

R = "/sietch_colab/data_share/illex/popgen_data/analysis/steps/05_karyo_pca_arms/refit_619"
M = "/sietch_colab/data_share/illex/popgen_data/illex_manuscript"
kar = pd.read_csv(f"{M}/data/qc/karyotypes_619.tsv", sep="\t", comment="#").set_index("sample")["karyotype_final"]


def eta2(x, lab):
    ok = np.isin(lab, ["AA", "AB", "BB"]) & np.isfinite(x)
    x, l = x[ok], lab[ok]
    sb = sum((l == k).sum() * (x[l == k].mean() - x.mean()) ** 2 for k in ["AA", "AB", "BB"])
    return sb / ((x - x.mean()) ** 2).sum()


# per-chromosome PCA: eta^2 and GMM screen
d = pd.read_csv(f"{R}/karyo_pca_coords_619.tsv", sep="\t", dtype={"chrom": str})
r = []
for c, g in d.groupby("chrom"):
    assert len(g) == 619
    x = g.PC1.to_numpy()
    lab = g["sample"].map(kar).to_numpy()
    X = x.reshape(-1, 1)
    bic = {k: GaussianMixture(k, n_init=4, random_state=0).fit(X).bic(X) for k in (1, 2, 3)}
    g3 = GaussianMixture(3, n_init=4, random_state=0).fit(X)
    sep = (g3.means_.max() - g3.means_.min()) / np.sqrt(g3.covariances_.ravel().mean())
    best = min((2, 3), key=lambda k: bic[k])
    r.append((c, eta2(x, lab), bic[1] - bic[best], best, sep))
r = pd.DataFrame(r, columns=["chrom", "eta2", "dBIC", "best_k", "sep"]).sort_values("dBIC", ascending=False)
r.to_csv(f"{R}/perchrom_eta2_gmm_619.tsv", sep="\t", index=False)
print(r.head(6).round(3).to_string())
aut = r[~r.chrom.isin(["2", "Z"])]
print("chr2 eta2", round(r[r.chrom == "2"].eta2.iloc[0], 4), "| max eta2 other", round(r[r.chrom != "2"].eta2.max(), 4),
      "| max dBIC other autosome", round(aut.dBIC.max(), 1), aut.loc[aut.dBIC.idxmax(), "chrom"], "| n chroms", len(r))

# windowed PCA chr2
names = open(f"{R}/chr2_619_samples.txt").read().split()
C = np.load(f"{R}/chr2_winpca_619_coords.npy")[:, :, 0]
assert C.shape[1] == len(names) == 619
w = pd.read_csv(f"{R}/chr2_winpca_619_windows.tsv", sep="\t")
ln = np.array([kar[i] for i in names])
en = np.array([eta2(C[i], ln) for i in range(len(w))])
inv = (w.start >= 60.54e6) & (w.end <= 79.50e6)
col = (w.end < 60.54e6 - 1e6) | (w.start > 79.50e6 + 1e6)
print("windowed inv mean/median", np.nanmean(en[inv]).round(4), np.nanmedian(en[inv]).round(4),
      "col mean/median", np.nanmean(en[col]).round(4), np.nanmedian(en[col]).round(4),
      "windows", int(inv.sum()), int(col.sum()))
