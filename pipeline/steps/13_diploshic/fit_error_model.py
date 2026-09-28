import json, numpy as np, allel, pandas as pd
WD="/sietch_colab/data_share/illex/popgen_data/analysis/steps/13_diploshic"
ANA="/sietch_colab/data_share/illex/popgen_data/analysis"
# per-sample mean depth (proxy) for the 633 baker samples (callset has GT only, no DP)
dp=pd.read_csv(f"{ANA}/steps/03_karyotype/depth.tsv", sep="\t")
dmap=dict(zip(dp["sample"], dp["depth_proxy"]))
samples=[s.strip() for s in open(f"{WD}/all_samples.txt")]
means=np.clip(np.array([dmap[s] for s in samples]), 0, 30)   # (633,) mean depth/sample
# per-genotype depth ~ Poisson(sample mean); pool many draws per sample -> depth histogram
rng=np.random.default_rng(0)
draws=rng.poisson(np.repeat(means, 3000)).clip(0, 60)
hist=np.bincount(draws, minlength=61).astype(float); hist/=hist.sum()
depth_hist=[[int(i), float(hist[i])] for i in range(61) if hist[i]>0]
# missing rate measured from the real combined callset GT (no DP needed)
miss=0; tot=0
for c in ["1","10","20","30","40"]:
    gt=allel.read_vcf(f"{WD}/combined/all.{c}.vcf.gz", fields=["calldata/GT"])["calldata/GT"]
    miss += int((gt[:,:,0]<0).sum()); tot += gt.shape[0]*gt.shape[1]
# het-dropout: at depth d, P(all reads same allele | true het) = 2*0.5^d (both-hom outcomes) -> called hom
het_drop=[[int(d), float(2*0.5**d) if d>0 else 1.0] for d in range(0,21)]
out=dict(n_samples=len(samples), miss_rate=float(miss/tot), depth_hist=depth_hist,
         het_dropout_by_depth=het_drop,
         filters=dict(biallelic_snp=True, max_f_missing=0.5, min_maf=0.01))
json.dump(out, open(f"{WD}/ascertainment.json","w"), indent=1)
print("n", len(samples), "miss_rate", out["miss_rate"], "mean_depth", float((draws).mean()))
