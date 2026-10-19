"""Supplementary Figs S2 and S19: per-SNP geographic Hudson FST across the 45 NAFO-division
pairs, every biallelic SNP with pooled MAF > 0.05, against a permutation null (division labels
shuffled among individuals). Inputs from pipeline/steps/17_ensemble_scan/scripts/perm_outlier_test.py."""
import os, sys, json
os.environ.setdefault("MPLCONFIGDIR", "/dev/shm/mplcache")
_HERE = os.path.dirname(os.path.abspath(__file__)); _ROOT = os.path.dirname(_HERE)
sys.path.insert(0, _HERE)
from figstyle import apply, C, despine
apply()
import numpy as np, pandas as pd, matplotlib.pyplot as plt

R = sys.argv[1] if len(sys.argv) > 1 else "/sietch_colab/data_share/illex/popgen_data/analysis/steps/17_ensemble_scan/fst_maf05/results"
OUT = os.path.join(_ROOT, "supp_figures")
S = json.load(open(f"{R}/fst_perm_summary.json")); T = S["thresholds"]
H = np.load(f"{R}/fst_hist.npz"); Q = np.load(f"{R}/fst_quantiles.npz")
mid = (H["edges"][:-1] + H["edges"][1:]) / 2; w = np.diff(H["edges"])
dens = lambda c: c / c.sum() / w
TL = {"q99.9": "99.9th pct", "q99.99": "99.99th pct", "q99.999": "99.999th pct", "null_max": "maximum"}

# ---------------- S2 ----------------
fig, ax = plt.subplots(1, 3, figsize=(10.8, 3.4))
a = ax[0]
a.step(mid, dens(H["null"]), where="mid", color=C["muted"], lw=1.3, label=f"permutation null ({S['n_permutations']} shuffles)")
a.step(mid, dens(H["real"]), where="mid", color=C["accent"], lw=1.3, label="observed")
a.set_yscale("log"); a.set_xlabel(r"per-SNP $F_{ST}$ (Hudson)"); a.set_ylabel("density")
for k in ("q99.99", "null_max"):
    a.axvline(T[k], color=C["ink"], ls=":", lw=0.8)
a.text(T["null_max"], a.get_ylim()[1] * 0.3, " null max", fontsize=6.5, ha="left")
a.legend(frameon=False, fontsize=6.5, loc="center right"); a.set_title("a   autosomes outside the inversion", loc="left", fontsize=8, fontweight="bold")
despine(a)
b = ax[1]
b.step(mid, dens(H["chr2_collinear"]), where="mid", color=C["muted"], lw=1.3, label="chr2 collinear")
b.step(mid, dens(H["inv"]), where="mid", color="#009E73", lw=1.3, label="chr2 inversion")
b.set_yscale("log"); b.set_xlabel(r"per-SNP $F_{ST}$ (Hudson)"); b.legend(frameon=False, fontsize=6.5)
b.set_title("b   chromosome 2", loc="left", fontsize=8, fontweight="bold"); despine(b)
c = ax[2]
g = S["groups"]["autosomes_outside_inversion"]; ks = ["q99.9", "q99.99", "q99.999"]
xo = np.arange(len(ks))
c.bar(xo - 0.18, [g[k]["expected"] for k in ks], 0.36, color=C["faint"], edgecolor=C["muted"], label="expected (null)")
c.bar(xo + 0.18, [g[k]["observed"] for k in ks], 0.36, color=C["accent"], label="observed")
c.set_yscale("log"); c.set_xticks(xo); c.set_xticklabels([TL[k] for k in ks], fontsize=6.5, rotation=20)
c.set_ylabel("SNPs above null threshold"); c.legend(frameon=False, fontsize=6.5)
c.text(0.98, 0.62, f"above null maximum ({T['null_max']:.3f}): {g['null_max']['observed']}", transform=c.transAxes, ha="right", fontsize=6.5)
c.set_title("c   observed vs expected", loc="left", fontsize=8, fontweight="bold"); despine(c)
fig.tight_layout(); fig.savefig(f"{OUT}/figS2_fst_expanded.png", dpi=220, bbox_inches="tight"); plt.close(fig)

# ---------------- S19 ----------------
P = os.path.join(os.path.dirname(R), "perm")   # observed per-SNP FST (perm 00 = real labels), every SNP
order = [str(i) for i in range(1, 46)] + ["Z"]
fr = []
for ch in order:
    f = np.load(f"{P}/00/{ch}.fst.npy"); ok = np.isfinite(f)
    fr.append(pd.DataFrame({"CHROM": ch, "POS": np.load(f"{P}/{ch}.pos.npy")[ok], "fst": f[ok]}))
pts = pd.concat(fr, ignore_index=True)
off, ticks, lab, cur = {}, [], [], 0
for ch in order:
    m = pts.CHROM == ch
    if not m.any(): continue
    mx = pts.POS[m].max(); off[ch] = cur; ticks.append(cur + mx / 2); lab.append(ch); cur += mx + 5e6
pts["x"] = pts.POS + pts.CHROM.map(off)
fig, ax = plt.subplots(1, 2, figsize=(11.5, 3.4), gridspec_kw={"width_ratios": [3.2, 1]})
a = ax[0]
for i, ch in enumerate(order):
    m = pts.CHROM == ch
    col = "#009E73" if False else (C["accent"] if ch == "Z" else (C["ink"] if i % 2 == 0 else C["muted"]))
    a.scatter(pts.x[m] / 1e6, pts.fst[m], s=1.2, color=col, lw=0, rasterized=True)
inv = (pts.CHROM == "2") & (pts.POS >= 60.54e6) & (pts.POS <= 79.5e6)
a.scatter(pts.x[inv] / 1e6, pts.fst[inv], s=1.2, color="#009E73", lw=0, rasterized=True)
for k, ls in (("q99.99", ":"), ("null_max", "--")):
    a.axhline(T[k], color=C["warm"], ls=ls, lw=0.8)
    a.text(cur / 1e6, T[k], f"null {TL[k]}", fontsize=6.5, va="bottom", ha="left", color=C["warm"])
a.set_xlim(-10, cur / 1e6 + 160)
tk = [i for i in range(len(lab)) if i % 2 == 0 or lab[i] == "Z"]
a.set_xticks(np.array([ticks[i] for i in tk]) / 1e6); a.set_xticklabels([lab[i] for i in tk], fontsize=5.5, rotation=90)
a.set_xlabel("chromosome"); a.set_ylabel(r"per-SNP $F_{ST}$ (Hudson)")
a.set_title("a   genome-wide per-SNP $F_{ST}$", loc="left", fontsize=8, fontweight="bold"); despine(a)
b = ax[1]
qs = Q["qs"]; t = qs >= 0.9
b.plot(Q["null"][t], Q["real"][t], color=C["accent"], lw=1.3, label="autosomes outside inversion")
b.plot(Q["null"][t], Q["inv"][t], color="#009E73", lw=1.3, label="chr2 inversion")
lim = [min(Q["null"][t].min(), Q["real"][t].min()), max(Q["null"][t].max(), Q["real"][t].max())]
b.plot(lim, lim, color=C["muted"], ls="--", lw=0.8)
b.set_xlabel("permutation null quantile"); b.set_ylabel("observed quantile"); b.legend(frameon=False, fontsize=6.5)
b.set_title("b   upper tail (90th-99.95th pct)", loc="left", fontsize=8, fontweight="bold"); despine(b)
fig.tight_layout(); fig.savefig(f"{OUT}/figS19_fst_genomewide.png", dpi=220, bbox_inches="tight"); plt.close(fig)
print("wrote figS2_fst_expanded.png and figS19_fst_genomewide.png")
