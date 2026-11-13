#!/usr/bin/env python
"""Refit the chr2 karyotyping PCA on the 619 analysed individuals only.

The original fits (recluster_karyo.py on chr2_karyo_coords.tsv, 653 individuals of the chr2 callset;
recluster_downsampled.py on the depth-equalised VCF, 633) included the 12 empty BAMs, the two
individuals excluded as contaminated and, at full depth, 20 small_2026 individuals. This script repeats
both fits with only the 619 in data/qc/karyotypes_619.tsv: bcftools subsets chr2:60-80 Mb to the 619
and keeps biallelic SNPs with MAF > 0.05 among them, pg_gpu randomized PCA (6 components, haplotype
pairs averaged per individual), then a 3-component GMM on PC1 (n_init 10, random_state 0), PC1
low->high = AA/AB/BB after orienting PC1 so that the old BB calls are positive.

Usage:
  pca   (pg_gpu env):       refit_pca_619.py pca OUTDIR
  calls (sklearn env):      refit_pca_619.py calls OUTDIR
"""
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd

A = Path("/sietch_colab/data_share/illex/popgen_data/analysis/steps")
M = Path(__file__).resolve().parents[3]
KAR = M / "data/qc/karyotypes_619.tsv"
BCF = "/home/ssmall/miniforge3/envs/bioinfo-buddy/bin/bcftools"
REGION = "2:60000000-80000000"
VCFS = {"full": A / "00_callset/filtered/2/variants_filt.vcf.gz",
        "ds": A / "03_karyotype/downsample/chr2inv_ds.vcf.gz"}


def ds_name(s):
    return s.split("/")[-1].replace(".ds.bam", "").replace(".bam", "")


def pca(out):
    from pg_gpu import HaplotypeMatrix, decomposition
    keep = pd.read_csv(KAR, sep="\t", comment="#")["sample"].tolist()
    for tag, vcf in VCFS.items():
        names = subprocess.check_output([BCF, "query", "-l", str(vcf)]).decode().split()
        raw = {ds_name(n): n for n in names}
        sel = [raw[s] for s in keep]
        (out / f"{tag}_samples.txt").write_text("\n".join(sel) + "\n")
        sub = out / f"{tag}_619.vcf.gz"
        subprocess.run(f"{BCF} view -r {REGION} -S {out}/{tag}_samples.txt -m2 -M2 -v snps {vcf} -Ou | "
                       f"{BCF} view -q 0.05:minor -Oz -o {sub} && {BCF} index -t -f {sub}",
                       shell=True, check=True)
        h = HaplotypeMatrix.from_vcf(str(sub), region=REGION)
        n = h.num_haplotypes // 2
        c, _ = decomposition.randomized_pca(h, n_components=6)
        c = c.get() if hasattr(c, "get") else np.asarray(c)
        if c.shape[0] == 2 * n:
            c = (c[:n] + c[n:2 * n]) / 2.0
        order = [ds_name(s) for s in subprocess.check_output([BCF, "query", "-l", str(sub)]).decode().split()]
        pd.DataFrame({"sample": order, "PC1": c[:, 0], "PC2": c[:, 1]}).to_csv(out / f"{tag}_coords.tsv", sep="\t", index=False)
        print(f"[{tag}] {h.num_variants} SNPs, {n} individuals", flush=True)


def calls(out):
    from sklearn.mixture import GaussianMixture
    kar = pd.read_csv(KAR, sep="\t", comment="#").set_index("sample")
    res = {}
    for tag, col in (("full", "karyotype"), ("ds", "karyotype_ds")):
        d = pd.read_csv(out / f"{tag}_coords.tsv", sep="\t").set_index("sample").loc[kar.index]
        x = d.PC1.to_numpy()
        if np.mean(x[kar[col] == "BB"]) < 0:
            x = -x
        gm = GaussianMixture(n_components=3, n_init=10, random_state=0).fit(x.reshape(-1, 1))
        o = np.argsort(gm.means_.ravel())
        name = {o[0]: "AA", o[1]: "AB", o[2]: "BB"}
        lab = pd.Series([name[l] for l in gm.predict(x.reshape(-1, 1))], index=kar.index)
        res[tag] = lab
        sizes = lab.value_counts().reindex(["AA", "AB", "BB"]).tolist()
        print(f"[{tag}] refit sizes AA/AB/BB {sizes}, agree with old {col}: {(lab == kar[col]).sum()}/619, "
              f"with karyotype_final: {(lab == kar.karyotype_final).sum()}/619")
    final = res["full"].where(res["full"] == res["ds"], res["ds"])
    dis = (res["full"] != res["ds"]).sum()
    sizes = final.value_counts().reindex(["AA", "AB", "BB"]).tolist()
    print(f"refit full vs ds disagree: {dis}; refit final (ds where they disagree) AA/AB/BB {sizes}; "
          f"agree with karyotype_final: {(final == kar.karyotype_final).sum()}/619")
    diff = kar.index[final != kar.karyotype_final]
    if len(diff):
        print(pd.DataFrame({"old_final": kar.karyotype_final[diff], "refit_full": res["full"][diff],
                            "refit_ds": res["ds"][diff]}).to_string())
    pd.DataFrame({"refit_full": res["full"], "refit_ds": res["ds"], "refit_final": final,
                  "karyotype_final": kar.karyotype_final}).to_csv(out / "refit_calls_619.tsv", sep="\t")


if __name__ == "__main__":
    out = Path(sys.argv[2]); out.mkdir(parents=True, exist_ok=True)
    {"pca": pca, "calls": calls}[sys.argv[1]](out)
