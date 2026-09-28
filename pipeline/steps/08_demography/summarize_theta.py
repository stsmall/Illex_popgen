#!/usr/bin/env python
# Genome-wide diversity summary from per-chrom thetaStat avg.pestPG (46 chroms, baker-633, folded, accessible).
# Per-chrom pi/thetaW/Tajima's D + pooled genome-wide (autosomes, excl chr2 inversion + chr42 sex when POOLING).
import glob, os
D="/sietch_colab/data_share/illex/popgen_data/analysis/steps/08_demography/thetas"
CHRS=[str(c) for c in list(range(1,46))]+["Z"]  # 1..45 + Z (no chr0); chr2=inversion, chr42=sex
rows={}
for c in CHRS:
    f=f"{D}/baker.{c}.avg.pestPG"
    if not os.path.exists(f): continue
    with open(f) as fh:
        fh.readline()
        parts=fh.readline().split()
        if len(parts)<14: continue
        tW=float(parts[3]); tP=float(parts[4]); tajd=float(parts[8]); n=int(parts[13])
        rows[c]=dict(tW=tW, tP=tP, tajd=tajd, n=n)
print(f"{'chr':>4} {'pi':>9} {'thetaW':>9} {'TajD':>7} {'nSites':>12}")
for c in CHRS:
    if c not in rows: continue
    r=rows[c]; pi=r['tP']/r['n']; tw=r['tW']/r['n']
    tag=" <-inv" if c=="2" else (" <-sex" if c=="42" else "")
    print(f"{c:>4} {pi:9.6f} {tw:9.6f} {r['tajd']:7.3f} {r['n']:12d}{tag}")
# pooled genome-wide = autosomes excluding chr2 (inversion) and chr42 (sex)
pool=[c for c in rows if c not in ("2","42","Z")]
sN=sum(rows[c]['n'] for c in pool); sP=sum(rows[c]['tP'] for c in pool); sW=sum(rows[c]['tW'] for c in pool)
# nSites-weighted mean Tajima's D
mTajD=sum(rows[c]['tajd']*rows[c]['n'] for c in pool)/sN
print(f"\n=== GENOME-WIDE (autosomes, excl chr2 inversion + chr42 sex; {len(pool)} chroms) ===")
print(f"pi      = {sP/sN:.6f}")
print(f"thetaW  = {sW/sN:.6f}")
print(f"TajD    = {mTajD:.4f} (nSites-weighted mean)")
print(f"nSites  = {sN:,} (sites with data)")
zc=rows.get('Z');
if zc: print(f"\nchrZ (reported separately): pi={zc['tP']/zc['n']:.6f} thetaW={zc['tW']/zc['n']:.6f} TajD={zc['tajd']:.3f}")
c2=rows.get('2')
if c2: print(f"chr2 (inversion, excluded from pool): pi={c2['tP']/c2['n']:.6f} TajD={c2['tajd']:.3f}")
