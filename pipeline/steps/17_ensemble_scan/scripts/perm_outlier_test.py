#!/usr/bin/env python
"""Permutation-calibrated per-SNP geographic Hudson FST outlier test (MAF > 0.05).

Null: per-SNP FST from division labels shuffled among the 632 individuals (perm/01..NN),
pooled genome-wide. Because every permutation scores exactly the same SNPs, the pooled
null is matched to the empirical data in allele frequency, sample size, missingness and
linkage. Thresholds: pooled null quantiles 99.9, 99.99, 99.999 and the pooled maximum.
Observed counts (real labels, perm/00) are compared with expected = n_SNPs x null tail rate,
separately for the autosomes outside the chr2 inversion (60.54-79.50 Mb), the inversion,
and chrZ. SNPs above the null 99.99th percentile are annotated as genic if they fall in a
gene interval (genes.bed). Usage: perm_outlier_test.py PERM_DIR GENES_BED OUT_DIR
"""
import sys, os, glob, json, numpy as np, pandas as pd
P, GENES, OUT = sys.argv[1:4]
os.makedirs(OUT, exist_ok=True)
chroms = [str(c) for c in range(1, 46)] + ["Z"]
perms = sorted(d for d in os.listdir(P) if d.isdigit() and d != "00")
real, pos, maf, chrom = [], [], [], []
for c in chroms:
    r = np.load(f"{P}/00/{c}.fst.npy"); real.append(r); pos.append(np.load(f"{P}/{c}.pos.npy"))
    maf.append(np.load(f"{P}/{c}.maf.npy")); chrom.append(np.full(len(r), c, dtype=object))
real = np.concatenate(real); pos = np.concatenate(pos); maf = np.concatenate(maf); chrom = np.concatenate(chrom)
inv = (chrom == "2") & (pos >= 60_540_000) & (pos <= 79_500_000)
groups = {"autosomes_outside_inversion": (chrom != "Z") & ~inv, "chr2_inversion": inv, "chrZ": chrom == "Z"}
# pooled null, autosomes outside the inversion (the reference for thresholds)
nullv = np.concatenate([np.concatenate([np.load(f"{P}/{k}/{c}.fst.npy") for c in chroms])[groups["autosomes_outside_inversion"]]
                        for k in perms])
nullv = nullv[np.isfinite(nullv)]
q = {"q99.9": .999, "q99.99": .9999, "q99.999": .99999}
thr = {k: float(np.quantile(nullv, v)) for k, v in q.items()}; thr["null_max"] = float(nullv.max())
rate = {k: float((nullv > t).mean()) for k, t in thr.items()}
res = {"n_permutations": len(perms), "n_null_values": int(len(nullv)), "null_median": float(np.median(nullv)),
       "null_sd": float(nullv.std()), "thresholds": thr, "null_tail_rate": rate, "groups": {}}
for g, m in groups.items():
    x = real[m & np.isfinite(real)]
    e = {"n_snps": int(len(x)), "median": float(np.median(x)), "sd": float(x.std()), "max": float(x.max())}
    for k, t in thr.items():
        e[k] = {"observed": int((x > t).sum()), "expected": float(len(x) * rate[k])}
    res["groups"][g] = e
# annotate SNPs above the null 99.99th percentile (autosomes outside the inversion and the inversion)
genes = pd.read_csv(GENES, sep="\t", names=["chrom", "start", "end", "gene"], dtype={"chrom": str})
hi = np.where(np.isfinite(real) & (real > thr["q99.99"]) & (chrom != "Z"))[0]
rows = []
for i in hi:
    gc = genes[(genes.chrom == chrom[i]) & (genes.start <= pos[i]) & (genes.end >= pos[i])]
    rows.append({"CHROM": chrom[i], "POS": int(pos[i]), "fst": round(float(real[i]), 5), "maf": round(float(maf[i]), 4),
                 "above_null_max": bool(real[i] > thr["null_max"]),
                 "overlap": "genic" if len(gc) else "intergenic", "gene": ";".join(gc.gene) if len(gc) else "."})
out = pd.DataFrame(rows, columns=["CHROM", "POS", "fst", "maf", "above_null_max", "overlap", "gene"])
out.to_csv(f"{OUT}/fst_outliers_q9999.tsv", sep="\t", index=False)
res["outliers_q9999"] = {"n": int(len(out)), "genic": int((out.overlap == "genic").sum()),
                         "above_null_max": int(out.above_null_max.sum()),
                         "genic_genes": sorted(set(";".join(out.gene[out.overlap == "genic"]).split(";")) - {""})}
json.dump(res, open(f"{OUT}/fst_perm_summary.json", "w"), indent=1)
# compact plotting data: every SNP above null q99.9, plus a 2% random subsample of the rest
rng = np.random.default_rng(1); ok = np.isfinite(real)
sel = ok & ((real > thr["q99.9"]) | (rng.random(len(real)) < 0.02))
pd.DataFrame({"CHROM": chrom[sel], "POS": pos[sel], "fst": real[sel].round(5)}).to_csv(f"{OUT}/fst_plot_points.tsv.gz", sep="\t", index=False)
qs = np.linspace(0, 1, 2001)[1:-1]
np.savez(f"{OUT}/fst_quantiles.npz", qs=qs, real=np.quantile(real[groups["autosomes_outside_inversion"] & ok], qs),
         null=np.quantile(nullv, qs), inv=np.quantile(real[inv & ok], qs),
         chr2_collinear=np.quantile(real[(chrom == "2") & ~inv & ok], qs))
h_edges = np.linspace(-0.06, 0.4, 231)
np.savez(f"{OUT}/fst_hist.npz", edges=h_edges, real=np.histogram(real[groups["autosomes_outside_inversion"] & ok], h_edges)[0],
         null=np.histogram(nullv, h_edges)[0], inv=np.histogram(real[inv & ok], h_edges)[0],
         chr2_collinear=np.histogram(real[(chrom == "2") & ~inv & ok], h_edges)[0])
print(json.dumps(res, indent=1))
