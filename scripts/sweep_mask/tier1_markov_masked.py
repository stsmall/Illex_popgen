#!/usr/bin/env python3
"""Tier-1 = Markov-corrected diploSHIC hard/soft calls, + enrichment.

The HMM (design-informed params) decodes each 100 kb window as N / linked (LL,LR) /
center C. In this 91%-linked genome the value of the Markov step is spatial coherence:
it retains raw hard/soft calls that sit in coherent sweep centers and drops spatially
isolated ones. Tier-1 = the retained C-state windows, merged into localized calls.

Ranks Tier-1 by the pre-computed empirical-outlier percentile S_pct (rank within
callable x cds strata) -- the project's empirical-FDR statistic -- and tests whether
Tier-1 is enriched at the Tier-2 cross-method-concordant loci.
"""
import sys, pandas as pd, numpy as np
from scipy.stats import fisher_exact
D = "/sietch_colab/data_share/illex/popgen_data/analysis/steps/14_sweep_seqmodel"
R = f"{D}/results/empirical_scan_fullsfs"
OUT = f"{R}/hmm_decode_masked"

win = pd.read_csv(f"{R}/outlier_scan_45_masked/windows.tsv", sep="\t")
win["chrom"] = win["chrom"].astype(str)
sp = pd.read_csv(f"{OUT}/state_path.tsv.gz", sep="\t")
sp["chrom"] = sp["chrom"].astype(str)
reg2 = pd.read_csv(f"{R}/outlier_scan_45_masked/regions.tsv", sep="\t")   # 234 Tier-2 candidates
reg2["chrom"] = reg2["chrom"].astype(str)

# raw diploSHIC argmax call
P = win[["p_neutral", "p_linked_hard", "p_linked_soft", "p_hard", "p_soft"]].to_numpy()
lab = np.array(["neutral", "linkedHard", "linkedSoft", "hard", "soft"])
win["raw_call"] = lab[P.argmax(1)]
win["pC"] = win["p_hard"] + win["p_soft"]
win["raw_hs"] = win["raw_call"].isin(["hard", "soft"])

# join HMM state
w = win.merge(sp[["chrom", "start", "state"]], on=["chrom", "start"], how="left")
w["markovC"] = w["state"] == "C"

# ---- Tier-2 concordant subsets ----
# concordant-34 = diploSHIC + >=1 cross-method (SF2 or RAiSD) = max_nmethods==2
reg2["concordant"] = (reg2["sf2_support"] > 0) | (reg2["raisd_support"] > 0)
n_conc = int(reg2["concordant"].sum())

# label each window by whether it lies in a Tier-2 candidate / concordant region
def in_region(w, regs):
    flag = np.zeros(len(w), bool)
    for _, r in regs.iterrows():
        m = (w["chrom"] == r["chrom"]) & (w["start"] < r["end"]) & (w["end"] > r["start"])
        flag |= m.to_numpy()
    return flag
w["in_tier2"] = in_region(w, reg2)
w["in_conc"] = in_region(w, reg2[reg2["concordant"]])

# ---- summary ----
raw_hs = int(w["raw_hs"].sum())
mk_hs = int(w["markovC"].sum())
print("========== Markov correction of diploSHIC hard/soft calls ==========")
print(f"  raw diploSHIC hard/soft windows : {raw_hs}  ({100*raw_hs/len(w):.2f}% of {len(w)})")
print(f"  Markov-retained (C-state)       : {mk_hs}  ({100*mk_hs/len(w):.2f}%)")
print(f"  removed as spatially isolated   : {raw_hs-mk_hs}  ({100*(raw_hs-mk_hs)/raw_hs:.0f}% of raw)")
kept = w[w["markovC"]]
print(f"  Markov hard : soft (peak pC)    : {int((kept['p_hard']>=kept['p_soft']).sum())} : {int((kept['p_hard']<kept['p_soft']).sum())}")

# merge C windows into localized Tier-1 calls (gap tolerance 1 window = 100kb)
tier1 = []
for c, g in w[w["markovC"]].sort_values(["chrom", "start"]).groupby("chrom"):
    starts = g["start"].to_numpy(); ends = g["end"].to_numpy()
    i = 0
    while i < len(starts):
        j = i
        while j + 1 < len(starts) and starts[j + 1] <= ends[j] + 100000:
            j += 1
        sub = g.iloc[i:j+1]
        tier1.append(dict(chrom=c, start=int(starts[i]), end=int(ends[j]), n_win=j-i+1,
                          peak_pC=float(sub["pC"].max()), max_Spct=float(sub["S_pct"].max()),
                          evidence="hard" if sub["p_hard"].sum() >= sub["p_soft"].sum() else "soft",
                          in_tier2=bool(sub["in_tier2"].any()), in_conc=bool(sub["in_conc"].any()),
                          chr2_42=c in ("2", "42")))
        i = j + 1
t1 = pd.DataFrame(tier1).sort_values("max_Spct", ascending=False).reset_index(drop=True)
t1.to_csv(f"{OUT}/tier1_markov_calls.tsv", sep="\t", index=False)
print(f"\n  Tier-1 Markov calls (merged C-runs): {len(t1)} localized regions")
print(f"    hard : soft = {int((t1.evidence=='hard').sum())} : {int((t1.evidence=='soft').sum())}")
print(f"    on chr2/chr42 : {int(t1.chr2_42.sum())}")
print(f"    in top-1% empirical stratum (max_Spct>=0.99): {int((t1.max_Spct>=0.99).sum())}")
print(f"    overlapping a Tier-2 candidate  : {int(t1.in_tier2.sum())}")
print(f"    overlapping a Tier-2 CONCORDANT : {int(t1.in_conc.sum())}  (of {n_conc} concordant regions)")

# ---- enrichment: is Markov C over-represented at concordant loci? (window-level Fisher) ----
a = int((w["markovC"] & w["in_conc"]).sum())      # C & concordant
b = int((w["markovC"] & ~w["in_conc"]).sum())     # C & not
c_ = int((~w["markovC"] & w["in_conc"]).sum())    # notC & concordant
d_ = int((~w["markovC"] & ~w["in_conc"]).sum())
OR, p = fisher_exact([[a, b], [c_, d_]], alternative="greater")
print(f"\n  Enrichment of Markov-C at concordant loci: OR={OR:.2f}, p={p:.2e}")
print(f"    ({a} of {a+c_} concordant-region windows are Markov-C = {100*a/max(a+c_,1):.0f}%, "
      f"vs {100*b/max(b+d_,1):.1f}% genome background)")
# S_pct enrichment: Markov C vs raw-only (dropped) calls
drop = w[w["raw_hs"] & ~w["markovC"]]
print(f"\n  median S_pct: Markov-kept C = {kept['S_pct'].median():.3f}  vs  raw-but-dropped = {drop['S_pct'].median():.3f}")
print(f"  (Markov retention favors higher empirical-outlier windows: keeps the real centers)")
print(f"\n  wrote {OUT}/tier1_markov_calls.tsv")
