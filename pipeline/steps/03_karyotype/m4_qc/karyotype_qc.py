#!/usr/bin/env python
"""Genotype quality and technical differences between chr2 karyotype groups (review M4).

Steps (run from analysis/steps/03_karyotype/m4_qc; plink2 and bcftools outputs are produced by karyotype_qc.sh):
  het      pooled and per-individual called Ho/He and call rate on chr20 (biallelic PASS SNPs of the 619), all SNPs
           and MAF >= 0.05. He per site is the unbiased 2pq*2n/(2n-1) from called genotypes, pooled as
           sum(Ho)/sum(He), the estimator of msinv illex/abc/error_model.genotype_error_stats. Writes
           perind_chr20_maf05.tsv and reports Ho/He and call rate against read depth and by karyotype.
  match    97 2a/2a individuals matched to the 97 2+/2+ on autosomal read count (greedy nearest neighbour) ->
           pops_depthmatched.tsv, for compute_flank_div_snps.py (POPS=..., OUT=chr2_40_100_windows_depthmatched.csv).
  ratios   pi_I/pi_S and dxy/pi_S over the 189 differentiated windows of the full-group table (FST > 0.15,
           chr2:60-80 Mb), jackknife over 1-Mb blocks, for the full and the depth-matched groups.
  glpost   per-individual karyotype posterior from ANGSD genotype likelihoods (gl_nearfixed.beagle.gz, -doMajorMinor 3
           with REF/ALT from the callset) at the biallelic near-fixed sites (|AF_2a - AF_2+| >= 0.9 between the
           breakpoints, nearfixed_af.tsv). P(data | K) = prod_sites sum_g GL(g) P(g | K), with P(g | K) from the
           arrangement allele frequencies (bounded to [0.01, 0.99]); prior from Hardy-Weinberg at p(2a) = 0.602.
Usage: karyotype_qc.py {het,match,ratios,sites,glpost}
"""
import gzip
import sys

import numpy as np
import pandas as pd
from scipy.stats import chi2_contingency, kruskal, mannwhitneyu, spearmanr

A = "/sietch_colab/data_share/illex/popgen_data/analysis/steps"
K = pd.read_csv("/sietch_colab/data_share/illex/popgen_data/illex_manuscript/data/qc/karyotypes_619.tsv", sep="\t",
                comment="#").set_index("sample").karyotype_final
READS = pd.read_csv(f"{A}/12_sex/sex_cov.tsv", sep="\t").set_index("sample").autosome_reads
DIV = pd.read_csv("/sietch_colab/data_share/illex/popgen_data/pggpu_illex/popstats/pops.tsv", sep="\t", header=None,
                  index_col=0)[1]


def pooled(h, m):
    n = h.HOM_A1_CT + h.HET_A1_CT + h.TWO_AX_CT
    p = (2 * h.HOM_A1_CT + h.HET_A1_CT) / (2 * n)
    he = (2 * p * (1 - p) * (2 * n) / (2 * n - 1) * n)[m].sum()
    return h.HET_A1_CT[m].sum() / he, (n[m] / 619).mean()


def het():
    h = pd.read_csv("chr20_all.hardy", sep="\t")
    n = h.HOM_A1_CT + h.HET_A1_CT + h.TWO_AX_CT
    p = (2 * h.HOM_A1_CT + h.HET_A1_CT) / (2 * n)
    maf = np.minimum(p, 1 - p)
    for lab, m in (("all", n >= 2), ("MAF>=0.05", (n >= 2) & (maf >= 0.05))):
        r, c = pooled(h, m)
        print(f"{lab}: {int(m.sum())} SNPs, pooled Ho/He {r:.3f}, mean call rate {c:.3f}")
    hz = pd.read_csv("chr20_maf05.het", sep="\t").set_index("#IID")
    ms = pd.read_csv("chr20_maf05.smiss", sep="\t").set_index("#IID").F_MISS
    d = pd.DataFrame({"K": K, "hohe": (hz.OBS_CT - hz["O(HOM)"]) / (hz.OBS_CT - hz["E(HOM)"]), "call": 1 - ms,
                      "reads": READS / 1e6, "div": DIV}).dropna(subset=["hohe"])
    d.to_csv("perind_chr20_maf05.tsv", sep="\t")
    print(f"per-individual Ho/He median {d.hohe.median():.3f}; Spearman with reads {spearmanr(d.hohe, d.reads)[0]:.2f}; "
          f"call rate with reads {spearmanr(d.call, d.reads)[0]:.2f}")
    print("Ho/He by read quintile", d.groupby(pd.qcut(d.reads, 5), observed=True).hohe.median().round(3).tolist())
    for v in ("reads", "call", "hohe"):
        g = [d[d.K == k][v] for k in ("AA", "AB", "BB")]
        print(f"{v}: medians {[round(x.median(), 3) for x in g]}, Kruskal p = {kruskal(*g)[1]:.3f}")
    ct = pd.crosstab(d["div"], d.K)
    ct = ct[ct.sum(1) >= 5]
    print(f"division x karyotype chi2 p = {chi2_contingency(ct)[1]:.3f} ({len(ct)} divisions with >= 5)")


def match():
    P = pd.read_csv(f"{A}/03_karyotype/karyofix/pops_karyofix.tsv", sep="\t", header=None, names=["s", "k"])
    aa, bb = P[P.k == "AA"].s.tolist(), P[P.k == "BB"].s.tolist()
    pool, pick = {s: READS[s] for s in aa}, []
    for s in sorted(bb, key=lambda x: READS[x]):
        t = min(pool, key=lambda a: abs(pool[a] - READS[s]))
        pick.append(t)
        pool.pop(t)
    pd.DataFrame({"s": pick + bb, "k": ["AA"] * len(pick) + ["BB"] * len(bb)}).to_csv(
        "pops_depthmatched.tsv", sep="\t", header=False, index=False)
    print("matched reads medians", np.median(READS[pick]) / 1e6, np.median(READS[bb]) / 1e6,
          "MW p", mannwhitneyu(READS[pick], READS[bb])[1])


def ratios():
    full = pd.read_csv(f"{A}/03_karyotype/karyofix/chr2_40_100_karyo_windows_snps_karyofix.csv")
    sel = full[(full.fst > 0.15) & (full.window_start >= 60_000_000) & (full.window_stop <= 80_000_000)].window_start

    def est(b):
        return b.pi_AA.sum() / b.pi_BB.sum(), b.dxy.sum() / b.pi_BB.sum()
    for lab, f in (("full 223/97", f"{A}/03_karyotype/karyofix/chr2_40_100_karyo_windows_snps_karyofix.csv"),
                   ("depth-matched 97/97", "chr2_40_100_windows_depthmatched.csv")):
        d = pd.read_csv(f)
        b = d[d.window_start.isin(sel)].copy()
        b["blk"] = b.window_start // 1_000_000
        th, bl = np.array(est(b)), b.blk.unique()
        jk = np.array([est(b[b.blk != k]) for k in bl])
        se = np.sqrt((len(bl) - 1) / len(bl) * ((jk - jk.mean(0)) ** 2).sum(0))
        c = d[(d.window_stop < 59_540_000) | (d.window_start > 80_500_000)]
        print(f"{lab}: pi_I/pi_S {th[0]:.3f} +- {se[0]:.3f}, dxy/pi_S {th[1]:.3f} +- {se[1]:.3f}, "
              f"collinear pi ratio {c.pi_AA.sum() / c.pi_BB.sum():.3f}")


def glpost():
    af = pd.read_csv("nearfixed_af.tsv", sep="\t").set_index("pos")
    with gzip.open("gl_nearfixed.beagle.gz", "rt") as fh:
        hdr = fh.readline().split()
        rows = [l.split() for l in fh]
    names = [s.split("/")[-1].replace(".bam", "") for s in open("bams619.txt").read().split()]
    assert len(hdr) == 3 + 3 * len(names)
    pos = np.array([int(r[0].split("_")[1]) for r in rows])
    gl = np.array([r[3:] for r in rows], float).reshape(len(rows), len(names), 3)   # REF/REF, REF/ALT, ALT/ALT
    keep = np.isin(pos, af.index)
    pos, gl = pos[keep], gl[keep]
    fa = np.clip(af.loc[pos, "fA"].to_numpy(), 0.01, 0.99)   # ALT frequency in 2a
    fb = np.clip(af.loc[pos, "fB"].to_numpy(), 0.01, 0.99)   # ALT frequency in 2+

    def pg(x, y):   # genotype probabilities (0, 1, 2 ALT copies) for one haplotype from each class
        return np.stack([(1 - x) * (1 - y), x * (1 - y) + (1 - x) * y, x * y], 1)
    ll = {}
    for k, (x, y) in {"AA": (fa, fa), "AB": (fa, fb), "BB": (fb, fb)}.items():
        ll[k] = np.log(np.einsum("sig,sg->si", gl, pg(x, y)) + 1e-300).sum(0)
    p = 0.602
    prior = {"AA": p * p, "AB": 2 * p * (1 - p), "BB": (1 - p) ** 2}
    L = np.stack([ll[k] + np.log(prior[k]) for k in ("AA", "AB", "BB")], 1)
    L = np.exp(L - L.max(1, keepdims=True))
    post = L / L.sum(1, keepdims=True)
    informative = (gl.max(2) - gl.min(2) > 1e-6).sum(0)
    out = pd.DataFrame(post, index=names, columns=["P_AA", "P_AB", "P_BB"])
    out["call_gl"] = out.idxmax(1).str[2:]
    out["P_max"] = out[["P_AA", "P_AB", "P_BB"]].max(1)
    out["sites_with_reads"] = informative
    out["karyotype_final"] = K.reindex(out.index)
    out.to_csv("gl_karyotype_posterior.tsv", sep="\t")
    agree = (out.call_gl == out.karyotype_final)
    print(f"{len(pos)} sites; agreement with final calls {agree.sum()}/{len(out)}; "
          f"P_max >= 0.99 for {(out.P_max >= 0.99).sum()}, median sites with reads {int(np.median(informative))}")
    print(out[~agree].round(3).to_string())


def sites():
    """Biallelic near-fixed sites between the breakpoints: |AF_2a - AF_2+| >= 0.9 from the karyofix allele-frequency
    tables (positions duplicated by multiallelic records dropped), alleles from inv_biallelic.tsv."""
    kf = f"{A}/03_karyotype/karyofix/"

    def af(f):
        d = pd.read_csv(kf + f, sep="\t", header=None, names=["c", "pos", "af"])
        d = d[(d.pos >= 60_540_000) & (d.pos <= 79_500_000)]
        return d.drop_duplicates("pos", keep=False).set_index("pos").af
    v = pd.read_csv("inv_biallelic.tsv", sep="\t", header=None, names=["pos", "ref", "alt"])
    v = v.drop_duplicates("pos", keep=False).set_index("pos")
    d = v.join(af("AA223_chr2_10_80_af.txt").rename("fA"), how="inner").join(af("BB97_chr2_10_80_af.txt").rename("fB"), how="inner")
    d = d[(d.fA - d.fB).abs() >= 0.9]
    d.reset_index().assign(c="2")[["c", "pos", "ref", "alt"]].to_csv("nearfixed_sites.txt", sep="\t", header=False, index=False)
    d.reset_index()[["pos", "ref", "alt", "fA", "fB"]].to_csv("nearfixed_af.tsv", sep="\t", index=False)
    with open("nearfixed_regions.txt", "w") as f:
        f.writelines(f"2:{p}-{p}\n" for p in d.index)
    print("near-fixed biallelic sites", len(d))


if __name__ == "__main__":
    {"het": het, "match": match, "ratios": ratios, "sites": sites, "glpost": glpost}[sys.argv[1]]()
