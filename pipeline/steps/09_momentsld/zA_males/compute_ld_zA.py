#!/usr/bin/env python
"""momentsLD two-locus statistics (pg_gpu, unphased genotypes) for the 330 ZZ males on chrZ and on the
autosomes, at fine distance bins (0-5 kb), for a Z vs autosome comparison of rho = 4 Ne r. Flat map
of 1 cM/Mb (1e-8 per bp), so r bins are bp x 1e-8. Autosomes: 2-Mb windows at 10 and 40 Mb on all 45
autosomes, skipping windows that overlap the chr2 inversion (60.54-79.50 Mb) or the 22 anomalous-
diversity blocks. Z: a 2-Mb window every 4 Mb. Accessible sites only. Writes ld_regions_{A,Z}.pkl."""
import pickle, time, subprocess, os, sys
import numpy as np
from pg_gpu.moments_ld import compute_ld_statistics
A = "/sietch_colab/data_share/illex/popgen_data/analysis"
MD = f"{A}/steps/09_momentsld/zA_males"
ACCESS = "/sietch_colab/data_share/illex/popgen_data/degenotate_illex/accessible_sites.bed"
BLK = f"{A}/steps/14_sweep_seqmodel/config/multicopy_block_mask.bed"
BB = "/home/ssmall/miniforge3/envs/bioinfo-buddy/bin"
BP = [0, 10, 25, 50, 100, 200, 500, 1000, 2000, 5000]
R_BINS = [b * 1e-8 for b in BP]
WIN = 2_000_000
fai = {l.split()[0]: int(l.split()[1]) for l in open("/sietch_colab/data_share/illex/alignments/anchorwave_dir/ill_coin_alignment/Illex_F24.primary.clean.fa.fai")}
blocks = [(l.split()[0], int(l.split()[1]), int(l.split()[2])) for l in open(BLK)]
def bad(c, s, e):
    if c == "2" and s < 79_500_000 and e > 60_540_000: return True
    return any(bc == c and s < be and e > bs for bc, bs, be in blocks)
group = sys.argv[1]
if group == "A":
    wins = [(str(c), s) for c in range(1, 46) for s in (10_000_000, 40_000_000)]
else:
    wins = [("Z", s) for s in range(3_000_000, fai["Z"] - WIN, 4_000_000)]
region = {}; i = 0
for c, start in wins:
    end = start + WIN
    if end > fai[c] or bad(c, start, end): print(f"skip {c}:{start}", flush=True); continue
    vcf = f"{A}/steps/00_callset/filtered/{c}/variants_filt.vcf.gz"
    tmp = f"{MD}/w_{c}_{start//1_000_000}.vcf.gz"
    subprocess.run(f"{BB}/bcftools view -S {MD}/../../12_sex/males.txt --force-samples -r {c}:{start}-{end} -m2 -M2 -v snps -Ou {vcf} "
                   f"| {BB}/bcftools view -c 1:minor -Oz -o {tmp} && {BB}/tabix -f -p vcf {tmp}", shell=True, stderr=subprocess.DEVNULL)
    nv = int(subprocess.run(f"{BB}/bcftools index -n {tmp}", shell=True, capture_output=True, text=True).stdout or 0)
    if nv < 500: os.remove(tmp); print(f"skip {c}:{start} nv={nv}", flush=True); continue
    t = time.time()
    try:
        d = compute_ld_statistics(vcf_file=tmp, rec_map_file=f"{MD}/recmap.flat.{c}.txt", pop_file=f"{MD}/pops_males.txt",
                                  pops=["male"], r_bins=R_BINS, use_genotypes=True, accessible_bed=ACCESS, report=False)
    except Exception as e:
        print(f"FAIL {c}:{start} nv={nv}: {e}", flush=True); os.remove(tmp); continue
    try:
        import cupy as cp
        for k, v in list(d.items()):
            if isinstance(v, cp.ndarray): d[k] = cp.asnumpy(v)
            elif isinstance(v, list): d[k] = [cp.asnumpy(x) if isinstance(x, cp.ndarray) else x for x in v]
    except Exception: pass
    region[i] = d; i += 1; os.remove(tmp); os.remove(tmp + ".tbi")
    print(f"region {i-1} = {c}:{start} nv={nv} {time.time()-t:.0f}s", flush=True)
pickle.dump({"regions": region, "r_bins": R_BINS, "bp_bins": BP, "rec_rate_per_bp": 1e-8, "mu": 3e-9},
            open(f"{MD}/ld_regions_{group}.pkl", "wb"))
print(f"SAVED {len(region)} regions ({group})", flush=True)
