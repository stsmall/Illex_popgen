#!/usr/bin/env python
# Genome-wide diversity summary for I. argentinus (10 low-cov dedup BAMs -> Illex_F24 ref),
# from per-chrom thetaStat avg.pestPG (folded, accessible sites). Mirrors the illecebrosus
# summarize_theta.py: per-chrom pi/thetaW/Tajima's D + pooled genome-wide (autosomes, excl
# chr2 inversion + chr42 sex; chrZ reported separately). Writes summary.tsv.
import os
D="/sietch_colab/data_share/illex/popgen_data/analysis/steps/08_demography/argentinus/thetas"
OUT="/sietch_colab/data_share/illex/popgen_data/analysis/steps/08_demography/argentinus/summary.tsv"
CHRS=[str(c) for c in range(1,46)]+["Z"]
ILL_PI=0.009302   # illecebrosus baker-633 genome-wide reference (same pipeline)
rows={}
for c in CHRS:
    f=f"{D}/arg.{c}.avg.pestPG"
    if not os.path.exists(f): continue
    with open(f) as fh:
        fh.readline()
        parts=fh.readline().split()
        if len(parts)<14: continue
        tW=float(parts[3]); tP=float(parts[4]); tajd=float(parts[8]); n=int(parts[13])
        if n==0: continue
        rows[c]=dict(tW=tW, tP=tP, tajd=tajd, n=n)
print(f"{'chr':>4} {'pi':>9} {'thetaW':>9} {'TajD':>7} {'nSites':>12}")
for c in CHRS:
    if c not in rows: continue
    r=rows[c]; pi=r['tP']/r['n']; tw=r['tW']/r['n']
    tag=" <-inv" if c=="2" else (" <-sex" if c=="42" else "")
    print(f"{c:>4} {pi:9.6f} {tw:9.6f} {r['tajd']:7.3f} {r['n']:12d}{tag}")
pool=[c for c in rows if c not in ("2","42","Z")]
sN=sum(rows[c]['n'] for c in pool); sP=sum(rows[c]['tP'] for c in pool); sW=sum(rows[c]['tW'] for c in pool)
mTajD=sum(rows[c]['tajd']*rows[c]['n'] for c in pool)/sN
gw_pi=sP/sN; gw_tw=sW/sN
print(f"\n=== GENOME-WIDE (autosomes, excl chr2 inversion + chr42 sex; {len(pool)} chroms) ===")
print(f"pi      = {gw_pi:.6f}")
print(f"thetaW  = {gw_tw:.6f}")
print(f"TajD    = {mTajD:.4f} (nSites-weighted mean)")
print(f"nSites  = {sN:,} (sites with data)")
print(f"\n vs I. illecebrosus (same pipeline): pi={ILL_PI:.6f} -> argentinus/illecebrosus = {gw_pi/ILL_PI:.3f}x")
zc=rows.get('Z')
if zc: print(f"chrZ (separate): pi={zc['tP']/zc['n']:.6f} thetaW={zc['tW']/zc['n']:.6f} TajD={zc['tajd']:.3f}")
c2=rows.get('2')
if c2: print(f"chr2 (inversion, excl): pi={c2['tP']/c2['n']:.6f} TajD={c2['tajd']:.3f}")
# write summary.tsv
with open(OUT,"w") as o:
    o.write("metric\tI_argentinus\tI_illecebrosus_ref\n")
    o.write(f"pi_genomewide\t{gw_pi:.6f}\t{ILL_PI:.6f}\n")
    o.write(f"thetaW_genomewide\t{gw_tw:.6f}\t0.031613\n")
    o.write(f"TajimasD_genomewide\t{mTajD:.4f}\t-2.0658\n")
    o.write(f"nSites_pooled\t{sN}\t1340915867\n")
    o.write(f"n_chroms_pooled\t{len(pool)}\t43\n")
    o.write(f"n_samples\t10\t633\n")
    if zc: o.write(f"chrZ_pi\t{zc['tP']/zc['n']:.6f}\t0.005577\n")
    if c2: o.write(f"chr2_pi_inv\t{c2['tP']/c2['n']:.6f}\t0.009066\n")
    o.write(f"pi_ratio_arg_over_ill\t{gw_pi/ILL_PI:.3f}\t1.000\n")
print(f"\nwrote {OUT}")
