#!/usr/bin/env python3
"""Apply the concordant33 per-region BGS test to all 658 Tier-1 Markov calls."""
import importlib.util, pandas as pd, numpy as np
spec=importlib.util.spec_from_file_location("bgs","/sietch_colab/data_share/illex/popgen_data/analysis/steps/14_sweep_seqmodel/calibration/bgs_screen_model.py")
m=importlib.util.module_from_spec(spec); spec.loader.exec_module(m)

A="/sietch_colab/data_share/illex/popgen_data/analysis"
OUT=f"{A}/steps/14_sweep_seqmodel/results/empirical_scan_fullsfs/hmm_decode"
t1=pd.read_csv(f"{OUT}/tier1_markov_calls.tsv",sep="\t",dtype={"chrom":str})
regs=[(r.chrom,int(r.start),int(r.end),r.evidence) for r in t1.itertuples()]
res=m.test_regions(regs)
bt=pd.DataFrame(res)[["chrom","start","end","sweep","pi","pi_pct","flank_cds","fc_pct","resid_vs_BGS","verdict"]]
# carry chr2_42 for reporting
bt=bt.merge(t1[["chrom","start","end","chr2_42"]],on=["chrom","start","end"],how="left")
bt.to_csv(f"{OUT}/tier1_BGS_test.tsv",sep="\t",index=False,
          columns=["chrom","start","end","sweep","pi","pi_pct","flank_cds","fc_pct","resid_vs_BGS","verdict"])

robust=bt[bt.verdict.str.startswith("SWEEP-like")].copy()
robust.to_csv(f"{OUT}/tier1_BGS_robust_calls.tsv",sep="\t",index=False,
              columns=["chrom","start","end","sweep","pi","pi_pct","flank_cds","fc_pct","resid_vs_BGS","verdict"])

print("==== Tier-1 per-region BGS screen ====")
print(f"total Tier-1 calls scored : {len(bt)}   (NA/no-pi: {int(bt.pi.isna().sum())})")
print("\nverdict breakdown:")
for v,n in bt.verdict.value_counts().items(): print(f"  {n:4d}  {v}")
print(f"\nBGS-robust (SWEEP-like, resid<=-0.3): {len(robust)}")
sub1=int((robust.verdict=='SWEEP-like (pi well below BGS expectation)').sum())
sub2=int((robust.verdict=='SWEEP-like (low pi, low flank-CDS: BGS unlikely)').sum())
print(f"    of which 'pi well below BGS expectation'      : {sub1}")
print(f"    of which 'low pi, low flank-CDS: BGS unlikely': {sub2}")
print(f"  hard : soft = {int((robust.sweep=='hard').sum())} : {int((robust.sweep=='soft').sum())}")
print(f"  on chr2 : {int((robust.chrom=='2').sum())}   on chr42 : {int((robust.chrom=='42').sum())}   chr2_or_42 flag: {int(robust.chr2_42.fillna(False).sum())}")
print(f"\nBGS-plausible (excluded) : {int((bt.verdict.str.startswith('BGS-plausible')).sum())}")
print(f"weak (excluded)          : {int((bt.verdict.str.startswith('weak')).sum())}")
print(f"\nwrote {OUT}/tier1_BGS_test.tsv  ({len(bt)} rows)")
print(f"wrote {OUT}/tier1_BGS_robust_calls.tsv  ({len(robust)} rows)")
