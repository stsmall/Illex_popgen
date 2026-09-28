#!/usr/bin/env python
"""Genome-wide FST downstream: aggregate per-chrom per-SNP mean FST, compare to the
already-computed coalescent neutral null, call/annotate genic outliers, GO-enrich,
build supp figure + result tables.  Run in mkado-vcf env."""
import os
os.environ['MPLCONFIGDIR'] = '/dev/shm/mplcache'
os.makedirs('/dev/shm/mplcache', exist_ok=True)
import sys, json, glob, subprocess
import numpy as np, pandas as pd
from collections import Counter
from scipy.stats import fisher_exact
sys.path.insert(0, "/sietch_colab/data_share/illex/popgen_data/analysis/manuscript")
from figstyle import apply, C, despine, save
import matplotlib.pyplot as plt

GW   = "/tmp/claude-1003/-sietch-colab-data-share-illex-popgen-data-mkado-illex/6d7afefb-6654-47e0-96f2-bc56ef699fbb/scratchpad/gwfst"
RESP = "/sietch_colab/data_share/illex/popgen_data/analysis/steps/17_ensemble_scan/results/fst_persnp"
STEP = "/sietch_colab/data_share/illex/popgen_data/analysis/steps/17_ensemble_scan/results"
FIG  = "/sietch_colab/data_share/illex/popgen_data/analysis/manuscript/supp_figures/supp_fst_genomewide.png"
CDSBED = "/sietch_colab/data_share/illex/popgen_data/degenotate_illex/illex.cds_sites.merge.bed"
BEDTOOLS = "/home/ssmall/miniforge3/envs/annotation-expression-buddy-barrnap_trnascan/bin/bedtools"

# ---------------- thresholds from the already-computed neutral null ----------------
summ    = json.load(open(f"{RESP}/sim_null_summary.json"))
thr_999 = summ["persnp_mean"]["q999"]      # 0.10499...
thr_max = summ["persnp_mean"]["max"]       # 0.32201...
sim     = np.load(f"{RESP}/sim_persnp_meanfst.npy")   # neutral per-SNP mean FST null
print(f"[null] per-SNP mean FST: q999={thr_999:.4f} max={thr_max:.4f} (n_sim={len(sim):,})")

# ---------------- load genome-wide empirical per-SNP mean FST ----------------
def chrom_key(c):
    return 100 if str(c) == "Z" else int(c)
files = sorted(glob.glob(f"{GW}/perchrom_*.tsv"), key=lambda f: chrom_key(f.split("perchrom_")[1].split(".tsv")[0]))
frames = [pd.read_csv(f, sep="\t", dtype={"CHROM": str}) for f in files]
emp = pd.concat(frames, ignore_index=True)
emp["POS"] = emp["POS"].astype(int)
emp["is_auto"] = emp["CHROM"] != "Z"
print(f"[emp] genome-wide per-SNP records: {len(emp):,} across {emp.CHROM.nunique()} chroms")
auto = emp[emp.is_auto]

# ---------------- outliers ----------------
out_max = emp[emp.meanf > thr_max].copy().sort_values("meanf", ascending=False)
out_999 = emp[emp.meanf > thr_999].copy().sort_values("meanf", ascending=False)
print(f"[outliers] > sim-max ({thr_max:.4f}): {len(out_max)}  (auto {int((out_max.is_auto).sum())})")
print(f"[outliers] > sim-q999 ({thr_999:.4f}): {len(out_999)}  (auto {int((out_999.is_auto).sum())})")

# ---------------- gene / CDS annotation ----------------
genes = pd.read_csv(f"{RESP}/genes.bed", sep="\t", names=["chrom","start","end","gene"], dtype={"chrom": str})
by_chrom = {c: g.sort_values("start").reset_index(drop=True) for c, g in genes.groupby("chrom")}

def genes_hit(df):
    hit = set()
    for c, sub in df.groupby("CHROM"):
        gc = by_chrom.get(str(c))
        if gc is None: continue
        starts, ends, gids = gc["start"].to_numpy(), gc["end"].to_numpy(), gc["gene"].to_numpy()
        for p in sub["POS"].to_numpy():
            k = np.searchsorted(starts, p, side="right") - 1
            j = k
            while j >= 0 and starts[j] >= p - 2_000_000:
                if starts[j] <= p <= ends[j]:
                    hit.add(gids[j])
                j -= 1
    return hit

_CDS = pd.read_csv(CDSBED, sep="\t", names=["chrom","start","end"], dtype={"chrom": str})
_CDS_BY = {c: (g["start"].to_numpy(), g["end"].to_numpy())
           for c, g in _CDS.sort_values(["chrom","start"]).groupby("chrom")}
def cds_distance(df):
    """bp from each SNP to nearest CDS interval (0 = inside a coding exon); -1 if chrom has no CDS."""
    dist = {}
    for _, r in df.iterrows():
        c, p = str(r["CHROM"]), int(r["POS"])
        se = _CDS_BY.get(c)
        if se is None:
            dist[(c, p)] = -1; continue
        starts, ends = se
        inside = (starts <= p) & (p <= ends)
        if inside.any():
            dist[(c, p)] = 0
        else:
            dist[(c, p)] = int(min(np.min(np.abs(starts - p)), np.min(np.abs(ends - p))))
    return dist

def annotate(df):
    rows = []
    dist = cds_distance(df)
    for c, sub in df.groupby("CHROM"):
        gc = by_chrom.get(str(c))
        for _, r in sub.iterrows():
            gene, ov = ".", "intergenic"
            if gc is not None:
                m = (gc["start"] <= r["POS"]) & (gc["end"] >= r["POS"])
                if m.any():
                    gene = ";".join(gc.loc[m, "gene"].tolist()); ov = "genic"
            rows.append(dict(CHROM=str(c), POS=int(r["POS"]), meanf=round(float(r["meanf"]),5),
                             maxf=round(float(r["maxf"]),5), npair=int(r["npair"]), overlap=ov,
                             dist_to_cds=dist.get((str(c), int(r["POS"])), -1), gene=gene))
    if not rows:
        return pd.DataFrame(columns=["CHROM","POS","meanf","maxf","npair","overlap","dist_to_cds","gene"])
    return pd.DataFrame(rows).sort_values("meanf", ascending=False)

# annotate lenient (q999) outliers -- superset of strict; mark which pass max
ann = annotate(out_999)
if len(ann):
    ann["exceeds_sim_max"] = ann["meanf"] > thr_max
ann.to_csv(f"{STEP}/fst_genomewide_outliers.tsv", sep="\t", index=False)
n_genic_999 = int((ann["overlap"] == "genic").sum()) if len(ann) else 0
n_genic_max = int(((ann["overlap"] == "genic") & (ann.get("exceeds_sim_max", False))).sum()) if len(ann) else 0
print(f"[annot] q999 outliers genic: {n_genic_999}/{len(ann)}; of the >max set genic: {n_genic_max}")

# ---------------- GO enrichment (genic outlier genes vs SNP-able background) ----------------
def go_enrich(fg_genes, bg_genes):
    g2go = {}
    with open(f"{RESP}/gene_go.tsv") as fh:
        for ln in fh:
            g, go = ln.rstrip("\n").split("\t"); g2go.setdefault(g, set()).add(go)
    term, cat = {}, {}
    with open(f"{RESP}/go_terms.tsv") as fh:
        for ln in fh:
            go, t, c = ln.rstrip("\n").split("\t"); term[go] = t; cat[go] = c
    fg = [g for g in fg_genes if g in g2go]
    bgset = set(g for g in bg_genes if g in g2go) - set(fg)
    Nfg, Nbg = len(fg), len(bgset)
    cols = ["go_id","term","category","n_fg","n_bg","N_fg","N_bg","fold","odds_ratio","p","FDR"]
    if Nfg == 0:
        return pd.DataFrame(columns=cols), Nfg, Nbg
    cf, cb = Counter(), Counter()
    for g in fg:
        for go in g2go[g]: cf[go] += 1
    for g in bgset:
        for go in g2go[g]: cb[go] += 1
    rows = []
    for go, a in cf.items():
        if a < 2: continue
        c_ = cb.get(go, 0); b = Nfg - a; d = Nbg - c_
        OR, p = fisher_exact([[a, b],[c_, d]], alternative="greater")
        fold = (a/Nfg) / max(c_/Nbg, 1e-12)
        rows.append([go, term.get(go,go), cat.get(go,"."), a, c_, Nfg, Nbg, round(fold,3), round(float(OR),3), p])
    if not rows:
        return pd.DataFrame(columns=cols), Nfg, Nbg
    R = pd.DataFrame(rows, columns=cols[:-1]).sort_values("p").reset_index(drop=True)
    m = len(R); R["FDR"] = (R["p"] * m / (np.arange(m)+1)).clip(upper=1.0)
    R["FDR"] = R["FDR"][::-1].cummin()[::-1]
    return R, Nfg, Nbg

bg_genes = genes_hit(emp)                                  # all SNP-able (scanned) genes = background
fg_genes = genes_hit(out_999) if len(out_999) else set()   # genes overlapped by outlier SNPs
GO, Nfg, Nbg = go_enrich(fg_genes, bg_genes)
GO.to_csv(f"{STEP}/fst_genomewide_GO.tsv", sep="\t", index=False)
n_sig = int((GO["FDR"] < 0.05).sum()) if len(GO) else 0
print(f"[GO] background SNP-able genes: {len(bg_genes)}; foreground outlier genes: {len(fg_genes)} "
      f"(annotated fg={Nfg}); terms tested: {len(GO)}; FDR<0.05: {n_sig}")

# ---------------- per-chrom + overall summary table ----------------
def pct_row(label, a, nsnp=None):
    a = np.asarray(a, float)
    return dict(chrom=label, n_snp=int(len(a)) if nsnp is None else nsnp,
                median=round(float(np.median(a)),5), q99=round(float(np.quantile(a,.99)),5),
                q999=round(float(np.quantile(a,.999)),5), max=round(float(np.max(a)),5),
                n_gt_q999=int((a > thr_999).sum()), n_gt_max=int((a > thr_max).sum()))
rows = []
for c in sorted(emp.CHROM.unique(), key=chrom_key):
    rows.append(pct_row(c, emp.loc[emp.CHROM==c, "meanf"].to_numpy()))
rows.append(pct_row("AUTOSOMES", auto["meanf"].to_numpy()))
rows.append(pct_row("GENOME_ALL", emp["meanf"].to_numpy()))
summ_df = pd.DataFrame(rows)
summ_df["sim_thr_q999"] = round(thr_999,5); summ_df["sim_thr_max"] = round(thr_max,5)
summ_df.to_csv(f"{STEP}/fst_genomewide_summary.tsv", sep="\t", index=False)
print(f"[table] wrote fst_genomewide_summary.tsv ({len(summ_df)} rows)")

# ================= FIGURE =================
apply()
fig, axes = plt.subplots(1, 2, figsize=(12.0, 4.4))

# (a) empirical vs neutral null
ax = axes[0]
hi = max(auto.meanf.max(), thr_max) * 1.05
bins = np.linspace(0, hi, 80)
ax.hist(sim, bins=bins, density=True, color=C["muted"], alpha=0.55,
        label=f"neutral sim null (n={len(sim):,})")
ax.hist(auto.meanf, bins=bins, density=True, histtype="step", color=C["accent"], lw=1.7,
        label=f"empirical autosomes (n={len(auto):,})")
ax.axvline(thr_999, color=C["warm"], lw=1.2, ls="--", label=f"sim 99.9th = {thr_999:.3f}")
ax.axvline(thr_max, color=C["warm"], lw=1.2, ls=":",  label=f"sim max = {thr_max:.3f}")
ax.set_xlabel("per-SNP mean $F_{ST}$ across 45 division pairs")
ax.set_ylabel("density")
ax.set_title("a  Genome-wide geographic $F_{ST}$ vs neutral null")
ax.legend(fontsize=7.5, loc="upper right")
despine(ax)

# (b) manhattan
ax = axes[1]
order = sorted(emp.CHROM.unique(), key=chrom_key)
xoff, ticks, tlab, x0 = {}, [], [], 0
for i, c in enumerate(order):
    n = int((emp.CHROM == c).sum())
    xoff[c] = x0; ticks.append(x0 + n/2); tlab.append(c); x0 += n + max(1, n//50)
cols = {c: (C["ink"] if i % 2 == 0 else C["muted"]) for i, c in enumerate(order)}
for c in order:
    sub = emp[emp.CHROM == c]
    xs = xoff[c] + np.arange(len(sub))
    ax.scatter(xs, sub["meanf"], s=1.4, c=(C["accent"] if c=="Z" else cols[c]), rasterized=True, linewidths=0)
if len(out_999):
    for c in out_999.CHROM.unique():
        sub = emp[emp.CHROM == c].reset_index(drop=True)
        mask = sub["meanf"] > thr_999
        ax.scatter(xoff[c] + np.where(mask)[0], sub.loc[mask, "meanf"], s=10, c=C["warm"], zorder=5, linewidths=0)
ax.axhline(thr_999, color=C["warm"], lw=1.0, ls="--", label=f"sim 99.9th = {thr_999:.3f}")
ax.axhline(thr_max, color=C["warm"], lw=1.0, ls=":",  label=f"sim max = {thr_max:.3f}")
ax.set_xticks(ticks[::2]); ax.set_xticklabels(tlab[::2], fontsize=5.5, rotation=90)
ax.set_xlabel("chromosome"); ax.set_ylabel("per-SNP mean $F_{ST}$")
ax.set_title("b  Per-SNP $F_{ST}$ across the genome")
ax.set_ylim(0, hi)
ax.legend(fontsize=7.5, loc="upper left")
despine(ax)

fig.tight_layout(rect=(0, 0.07, 1, 1))
foot = ("Genome-wide per-SNP Weir & Cockerham $F_{ST}$ across 45 NAFO-division pairs on the FULL baker_2025 callsets "
        "(45 autosomes + chrZ; each chrom randomly thinned to 2% of SNPs, seed 20240915), obs$\\geq$40/pair, mean over "
        "$\\geq$20 informative pairs, per-pair $F_{ST}$ clipped to [0,1].  Neutral null: single panmictic population under "
        "the fitted 2-epoch growth demography ($N_0{=}6.8{\\times}10^6$, $N_{ANC}{=}5.5{\\times}10^5$, $T{=}7.7{\\times}10^5$ gen, "
        "$\\mu{=}3{\\times}10^{-9}$) at empirical per-division $n$ with matched missingness, individuals randomly split into 10 "
        "\"divisions\".  chrZ (blue in b) is a sex chromosome and differentiation there is confounded by division sex ratio.")
fig.text(0.5, 0.012, foot, ha="center", va="bottom", fontsize=6.2, color="#555", wrap=True)
save(fig, FIG)

# ---------------- console report ----------------
print("\n================= GENOME-WIDE FST REPORT =================")
print(json.dumps({
  "n_snp_genome": int(len(emp)), "n_snp_auto": int(len(auto)), "n_chrom": int(emp.CHROM.nunique()),
  "n_genes_covered_background": len(bg_genes),
  "empirical_auto": {"median": round(float(auto.meanf.median()),4), "q99": round(float(auto.meanf.quantile(.99)),4),
                     "q999": round(float(auto.meanf.quantile(.999)),4), "max": round(float(auto.meanf.max()),4)},
  "empirical_genome": {"median": round(float(emp.meanf.median()),4), "max": round(float(emp.meanf.max()),4)},
  "threshold_q999": round(thr_999,4), "threshold_max": round(thr_max,4),
  "n_gt_q999_genome": int(len(out_999)), "n_gt_q999_auto": int((out_999.is_auto).sum()) if len(out_999) else 0,
  "n_gt_max_genome": int(len(out_max)), "n_gt_max_auto": int((out_max.is_auto).sum()) if len(out_max) else 0,
  "expected_gt_q999_by_chance": round(0.001*len(emp),1),
  "n_genic_outliers_q999": n_genic_999, "n_GO_FDR05": n_sig,
}, indent=2))
