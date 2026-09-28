#!/usr/bin/env python
"""
Downstream of the FST simulation null:
  - compare empirical per-SNP mean geographic FST to the neutral sim null
  - call outliers exceeding the sim threshold (99.9th pct / max)
  - annotate genic outliers against genes.bed, GO-enrich vs SNP-able background
  - build supplement figure + result tables

Run with the mkado-vcf env (pandas/numpy/scipy/matplotlib).
"""
import os
os.environ['MPLCONFIGDIR'] = '/dev/shm/mplcache'
os.makedirs('/dev/shm/mplcache', exist_ok=True)
import sys, json
import numpy as np, pandas as pd
from scipy.stats import fisher_exact
sys.path.insert(0, "/sietch_colab/data_share/illex/popgen_data/analysis/manuscript")
from figstyle import apply, C, despine, save
import matplotlib.pyplot as plt

KING = "/tmp/claude-1003/-sietch-colab-data-share-illex-popgen-data-mkado-illex/6d7afefb-6654-47e0-96f2-bc56ef699fbb/scratchpad/king"
RES  = "/sietch_colab/data_share/illex/popgen_data/analysis/steps/17_ensemble_scan/results/fst_persnp"
STEP = "/sietch_colab/data_share/illex/popgen_data/analysis/steps/17_ensemble_scan/results"
FIG  = "/sietch_colab/data_share/illex/popgen_data/analysis/manuscript/supp_figures/supp_fst_expanded.png"
INV0, INV1 = 60_540_000, 79_500_000

# ---------------- load ----------------
# Recompute empirical per-SNP mean FST from the 45 plink2 .fst.var files with the
# SAME filters/clip as the simulation (obs>=40 per pair; per-pair FST clipped to
# [0,1]; per-SNP mean across pairs; require >=20 informative pairs).
import glob
def emp_persnp():
    files = sorted(glob.glob(f"{KING}/fst_snp.*.fst.var"))
    frames = []
    for f in files:
        d = pd.read_csv(f, sep="\t", usecols=[0,1,3,4],
                        names=["CHROM","POS","obs","fst"], header=0)
        d = d[d["obs"] >= 40].dropna(subset=["fst"])
        d["fst"] = np.clip(d["fst"].to_numpy(), 0, 1.0)
        frames.append(d[["CHROM","POS","fst"]])
    allf = pd.concat(frames, ignore_index=True)
    g = allf.groupby(["CHROM","POS"], sort=False).agg(
        meanf=("fst","mean"), maxf=("fst","max"), npair=("fst","size")).reset_index()
    return g[g["npair"] >= 20].reset_index(drop=True)
emp = emp_persnp()
emp["CHROM"] = emp["CHROM"].astype(int)
sim = np.load(f"{RES}/sim_persnp_meanfst.npy")                      # per-SNP mean FST null
sim_pool = np.load(f"{RES}/sim_pooled_fst.npy")
emp_pool = np.load(f"{KING}/pooled_fst_finite.npy")
chr2 = pd.read_csv(f"{RES}/chr2_persnp_fst.tsv", sep="\t")
summ = json.load(open(f"{RES}/sim_null_summary.json"))

thr_999 = float(np.quantile(sim, .999))
thr_max = float(sim.max())
print(f"[thresh] sim per-SNP mean FST: median={np.median(sim):.4f} q99={np.quantile(sim,.99):.4f} "
      f"q999={thr_999:.4f} max={thr_max:.4f} (n_sim={len(sim):,})")
print(f"[emp]  per-SNP mean FST: median={emp.meanf.median():.4f} q99={emp.meanf.quantile(.99):.4f} "
      f"q999={emp.meanf.quantile(.999):.4f} max={emp.meanf.max():.4f} (n_emp={len(emp):,})")

# ---------------- outliers ----------------
out_max = emp[emp.meanf > thr_max].copy().sort_values("meanf", ascending=False)
out_999 = emp[emp.meanf > thr_999].copy()
print(f"[outliers] empirical SNPs > sim-max({thr_max:.4f}): {len(out_max)}")
print(f"[outliers] empirical SNPs > sim-q999({thr_999:.4f}): {len(out_999)}")

# ---------------- gene annotation ----------------
genes = pd.read_csv(f"{RES}/genes.bed", sep="\t", names=["chrom","start","end","gene"])
genes["chrom"] = genes["chrom"].astype(str)
by_chrom = {c: g.sort_values("start").reset_index(drop=True) for c, g in genes.groupby("chrom")}

def genes_hit(df):
    """set of gene IDs overlapped by SNP positions in df (CHROM,POS)."""
    hit = set()
    for c, sub in df.groupby("CHROM"):
        gc = by_chrom.get(str(int(c)))
        if gc is None: continue
        starts = gc["start"].to_numpy(); ends = gc["end"].to_numpy(); gids = gc["gene"].to_numpy()
        pos = sub["POS"].to_numpy()
        # for each gene, count SNPs inside; simpler: for each SNP find candidate genes
        idx = np.searchsorted(starts, pos, side="right") - 1
        for k, p in zip(idx, pos):
            # walk back a few genes in case of overlap
            j = k
            while j >= 0 and starts[j] >= p - 2_000_000:
                if starts[j] <= p <= ends[j]:
                    hit.add(gids[j])
                j -= 1
    return hit

def annotate(df):
    rows = []
    for c, sub in df.groupby("CHROM"):
        gc = by_chrom.get(str(int(c)))
        for _, r in sub.iterrows():
            gene, ov = ".", "intergenic"
            if gc is not None:
                m = (gc["start"] <= r["POS"]) & (gc["end"] >= r["POS"])
                if m.any():
                    gene = ";".join(gc.loc[m, "gene"].tolist()); ov = "genic"
            rows.append(dict(CHROM=int(c), POS=int(r["POS"]), meanf=r["meanf"],
                             maxf=r["maxf"], npair=int(r["npair"]), overlap=ov, gene=gene))
    return pd.DataFrame(rows).sort_values("meanf", ascending=False)

ann = annotate(out_max) if len(out_max) else pd.DataFrame(
        columns=["CHROM","POS","meanf","maxf","npair","overlap","gene"])
ann.to_csv(f"{STEP}/fst_outliers_annotated.tsv", sep="\t", index=False)
n_genic = int((ann["overlap"] == "genic").sum()) if len(ann) else 0
print(f"[annot] genic outliers: {n_genic} of {len(ann)}")

# ---------------- GO enrichment ----------------
def go_enrich(fg_genes, bg_genes):
    g2go, term, cat = {}, {}, {}
    with open(f"{RES}/gene_go.tsv") as fh:
        for ln in fh:
            g, go = ln.rstrip("\n").split("\t")
            g2go.setdefault(g, set()).add(go)
    with open(f"{RES}/go_terms.tsv") as fh:
        for ln in fh:
            go, t, c = ln.rstrip("\n").split("\t")
            term[go] = t; cat[go] = c
    fg = [g for g in fg_genes if g in g2go]
    bg = [g for g in bg_genes if g in g2go]
    bgset = set(bg) - set(fg)
    Nfg, Nbg = len(fg), len(bgset)
    if Nfg == 0:
        return pd.DataFrame(columns=["go_id","term","category","n_fg","n_bg",
                                     "N_fg","N_bg","fold","odds_ratio","p","FDR"])
    from collections import Counter
    cf, cb = Counter(), Counter()
    for g in fg:
        for go in g2go[g]: cf[go] += 1
    for g in bgset:
        for go in g2go[g]: cb[go] += 1
    rows = []
    for go, a in cf.items():
        if a < 2: continue
        c_ = cb.get(go, 0)
        b = Nfg - a; d = Nbg - c_
        OR, p = fisher_exact([[a, b], [c_, d]], alternative="greater")
        fold = (a/Nfg) / max(c_/Nbg, 1e-12)
        rows.append([go, term.get(go, go), cat.get(go, "."), a, c_, Nfg, Nbg,
                     round(fold, 3), round(float(OR), 3), p])
    if not rows:
        return pd.DataFrame(columns=["go_id","term","category","n_fg","n_bg",
                                     "N_fg","N_bg","fold","odds_ratio","p","FDR"])
    R = pd.DataFrame(rows, columns=["go_id","term","category","n_fg","n_bg",
                                    "N_fg","N_bg","fold","odds_ratio","p"])
    R = R.sort_values("p").reset_index(drop=True)
    m = len(R); R["FDR"] = (R["p"] * m / (np.arange(m)+1)).clip(upper=1.0)
    R["FDR"] = R["FDR"][::-1].cummin()[::-1]
    return R

bg_genes = genes_hit(emp)                                   # SNP-able background
fg_genes = genes_hit(out_max) if len(out_max) else set()
print(f"[GO] background genes (SNP-able): {len(bg_genes)}; foreground outlier genes: {len(fg_genes)}")
GO = go_enrich(fg_genes, bg_genes)
GO.to_csv(f"{STEP}/fst_outlier_GO.tsv", sep="\t", index=False)
n_sig = int((GO["FDR"] < 0.05).sum()) if len(GO) else 0
print(f"[GO] terms tested: {len(GO)}; FDR<0.05: {n_sig}")

# ---------------- sim-null summary table ----------------
def pct(a):
    return dict(n=int(len(a)), median=float(np.median(a)), mean=float(np.mean(a)),
                q99=float(np.quantile(a,.99)), q999=float(np.quantile(a,.999)),
                max=float(np.max(a)))
snull = pd.DataFrame([
    dict(statistic="per-SNP mean FST (across pairs)", source="simulation_null", **pct(sim)),
    dict(statistic="per-SNP mean FST (across pairs)", source="empirical",        **pct(emp.meanf.to_numpy())),
    dict(statistic="per-SNP mean FST (across pairs)", source="empirical_chr2",    **pct(chr2.meanf.to_numpy())),
    dict(statistic="per-SNP mean FST (chr2 inversion)", source="empirical_chr2_inv",
         **pct(chr2.loc[chr2.in_inversion,"meanf"].to_numpy())),
    dict(statistic="per-SNP-per-pair FST", source="simulation_null", **pct(sim_pool)),
    dict(statistic="per-SNP-per-pair FST", source="empirical",        **pct(emp_pool[np.isfinite(emp_pool)])),
])
snull["outlier_thr_q999"] = thr_999
snull["outlier_thr_max"]  = thr_max
snull["n_emp_gt_q999"] = len(out_999)
snull["n_emp_gt_max"]  = len(out_max)
snull.to_csv(f"{STEP}/persnp_fst_sim_null.tsv", sep="\t", index=False)
print(f"[table] wrote persnp_fst_sim_null.tsv")

# ================= FIGURE =================
apply()
fig, axes = plt.subplots(1, 3, figsize=(13.5, 4.0))

# ---- (a) empirical vs sim null, per-SNP mean FST ----
ax = axes[0]
bins = np.linspace(0, max(emp.meanf.max(), thr_max)*1.05, 70)
ax.hist(sim, bins=bins, density=True, color=C["muted"], alpha=0.55,
        label=f"neutral sim null (n={len(sim):,})")
ax.hist(emp.meanf, bins=bins, density=True, histtype="step", color=C["accent"],
        lw=1.6, label=f"empirical (n={len(emp):,})")
ax.axvline(thr_999, color=C["warm"], lw=1.2, ls="--",
           label=f"sim 99.9th = {thr_999:.3f}")
ax.axvline(thr_max, color=C["warm"], lw=1.2, ls=":",
           label=f"sim max = {thr_max:.3f}")
ax.set_xlabel("per-SNP mean $F_{ST}$ across 45 division pairs")
ax.set_ylabel("density")
ax.set_title("a  Geographic $F_{ST}$ vs neutral expectation")
ax.legend(fontsize=7, loc="upper right")
despine(ax)

# ---- (b) chr2 / inversion ----
ax = axes[1]
inv = chr2.loc[chr2.in_inversion, "meanf"].to_numpy()
non = chr2.loc[~chr2.in_inversion, "meanf"].to_numpy()
bins2 = np.linspace(0, max(chr2.meanf.max(), thr_max)*1.05, 60)
ax.hist(non, bins=bins2, density=True, color=C["faint"], alpha=0.9,
        label=f"chr2 non-inversion (n={len(non):,})")
ax.hist(inv, bins=bins2, density=True, histtype="step", color="#009E73", lw=1.8,
        label=f"chr2 inversion 60.5–79.5 Mb (n={len(inv):,})")
ax.axvline(thr_999, color=C["warm"], lw=1.2, ls="--", label="sim 99.9th")
ax.axvline(np.median(inv), color="#009E73", lw=1.0, ls=":", alpha=0.8)
ax.set_xlabel("per-SNP mean $F_{ST}$ across 45 division pairs")
ax.set_ylabel("density")
ax.set_title("b  chr2 inversion is non-clinal")
ax.legend(fontsize=7, loc="upper right")
despine(ax)
ax.text(0.97, 0.55, f"inversion median = {np.median(inv):.4f}\ninversion max = {inv.max():.3f}",
        transform=ax.transAxes, ha="right", va="top", fontsize=7.5, color="#222")

# ---- (c) outliers / annotation ----
ax = axes[2]
if len(out_max) == 0:
    ax.axis("off")
    ax.set_title("c  Genic outliers")
    msg = ("No geographic $F_{ST}$ outlier survives\n"
           "the neutral simulation threshold.\n\n"
           f"Empirical max per-SNP mean $F_{{ST}}$ = {emp.meanf.max():.3f}\n"
           f"Neutral sim max = {thr_max:.3f}\n"
           f"Neutral sim 99.9th pct = {thr_999:.3f}\n\n"
           f"{len(out_999)} SNPs exceed the sim 99.9th pct\n"
           f"(expected ≈ {0.001*len(emp):.0f} by chance).\n\n"
           "All geographic differentiation is\nconsistent with drift + finite\n"
           "sampling under the fitted growth\ndemography — no local adaptation\nsignal, no functional enrichment.")
    ax.text(0.02, 0.92, msg, transform=ax.transAxes, ha="left", va="top",
            fontsize=9, color="#222", linespacing=1.5)
else:
    top = ann.head(20)
    col = np.where(top["overlap"].to_numpy() == "genic", C["warm"], C["muted"])
    y = np.arange(len(top))[::-1]
    ax.scatter(top["meanf"], y, c=col, s=22, zorder=3)
    ax.axvline(thr_max, color=C["warm"], lw=1.0, ls=":")
    ax.set_yticks(y)
    ax.set_yticklabels([f"chr{c}:{p//1000}kb" for c, p in zip(top["CHROM"], top["POS"])],
                       fontsize=6)
    ax.set_xlabel("per-SNP mean $F_{ST}$")
    ax.set_title(f"c  {len(out_max)} outliers > sim max ({n_genic} genic)")
    despine(ax)

fig.tight_layout(rect=(0, 0.06, 1, 1))
foot = ("Neutral null: single panmictic population under the fitted moments 2-epoch growth demography "
        "($N_0{=}6.8{\\times}10^6$, $N_{ANC}{=}5.5{\\times}10^5$, $T{=}7.7{\\times}10^5$ gen, $\\mu{=}3{\\times}10^{-9}$), "
        "at the empirical per-division sample sizes with matched 44% per-sample missingness; individuals randomly split into the same 10 "
        "\"divisions\" (no real structure). Per-SNP WC $F_{ST}$ across 45 pairs, obs≥40/pair, mean over ≥20 informative pairs, "
        "per-pair $F_{ST}$ clipped to [0,1].  (a) autosomes (chr2 excluded); (b) chr2 full callset (thinned).")
fig.text(0.5, 0.015, foot, ha="center", va="bottom", fontsize=6.3, color="#555", wrap=True)
np.save(f"{RES}/emp_persnp_meanfst.npy", emp["meanf"].to_numpy().astype(np.float32))
save(fig, FIG)
print("[fig] wrote", FIG)

# ---------------- console report ----------------
print("\n================= REPORT =================")
print(json.dumps({
  "sim_null_persnp_mean": {k: round(v,4) for k,v in pct(sim).items()},
  "empirical_persnp_mean": {k: round(v,4) for k,v in pct(emp.meanf.to_numpy()).items()},
  "chr2_all": {k: round(v,4) for k,v in pct(chr2.meanf.to_numpy()).items()},
  "chr2_inversion": {k: round(v,4) for k,v in pct(chr2.loc[chr2.in_inversion,'meanf'].to_numpy()).items()},
  "threshold_q999": round(thr_999,4), "threshold_max": round(thr_max,4),
  "n_emp_gt_q999": len(out_999), "n_emp_gt_max": len(out_max),
  "n_genic_outliers": n_genic, "n_GO_FDR05": n_sig,
}, indent=2))
