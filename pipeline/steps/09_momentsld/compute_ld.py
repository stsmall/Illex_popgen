#!/usr/bin/env python
# momentsLD LD stats via pg_gpu over SPACED 2Mb windows (bootstrap regions), chr13-18, baker-633,
# accessible, flat 1cM/Mb, unphased. Dense data (24% segregating) -> 2Mb windows keep pair-count tractable.
import pickle, time, subprocess, os
from pg_gpu.moments_ld import compute_ld_statistics
A="/sietch_colab/data_share/illex/popgen_data/analysis"
MD=f"{A}/steps/09_momentsld"
ACCESS="/sietch_colab/data_share/illex/popgen_data/degenotate_illex/accessible_sites.bed"
BB="/home/ssmall/miniforge3/envs/bioinfo-buddy/bin"
R_BINS=[0,1e-6,2e-6,5e-6,1e-5,2e-5,5e-5,1e-4,2e-4,5e-4,1e-3]
WIN=2_000_000
STARTS=[5,20,35,50,65,80]        # Mb window-starts (spaced ~15Mb -> independent); 6/chrom
CHRS=[13,14,15,16,17,18]
region={}; i=0
for c in CHRS:
    vcf=f"{A}/steps/00_callset/filtered/{c}/variants_filt.vcf.gz"
    for s in STARTS:
        start=s*1_000_000; end=start+WIN
        tmp=f"{MD}/w_{c}_{s}.vcf.gz"
        subprocess.run(f"{BB}/bcftools view -r {c}:{start}-{end} -Oz -o {tmp} {vcf} && {BB}/tabix -f -p vcf {tmp}",
                       shell=True, stderr=subprocess.DEVNULL)
        nv=int(subprocess.run(f"{BB}/bcftools index -n {tmp}",shell=True,capture_output=True,text=True).stdout or 0)
        if nv<500: os.remove(tmp); print(f"skip {c}:{s}Mb (nv={nv})",flush=True); continue
        t=time.time()
        try:
            d=compute_ld_statistics(vcf_file=tmp, rec_map_file=f"{MD}/recmap.flat.{c}.txt",
                pop_file=f"{MD}/pops.txt", pops=["illex"], r_bins=R_BINS,
                use_genotypes=True, accessible_bed=ACCESS, report=False)
        except Exception as e:
            print(f"FAIL {c}:{s}Mb nv={nv}: {e}",flush=True); os.remove(tmp); continue
        try:
            import cupy as cp
            for k,v in list(d.items()):
                if isinstance(v,cp.ndarray): d[k]=cp.asnumpy(v)
                elif isinstance(v,list): d[k]=[cp.asnumpy(x) if isinstance(x,cp.ndarray) else x for x in v]
        except Exception: pass
        region[i]=d; i+=1
        os.remove(tmp)
        print(f"region {i-1} = chr{c}:{s}Mb nv={nv} {time.time()-t:.0f}s",flush=True)
with open(f"{MD}/ld_regions.pkl","wb") as fh:
    pickle.dump({"regions":region,"r_bins":R_BINS,"mu":3e-9,"rec_rate_per_bp":1e-8},fh)
print(f"SAVED {len(region)} regions -> ld_regions.pkl",flush=True)
