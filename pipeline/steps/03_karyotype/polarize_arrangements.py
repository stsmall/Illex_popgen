#!/usr/bin/env python
# Tally which arrangement (AA vs BB) is ancestral at diagnostic SNPs.
import pandas as pd, numpy as np
O="/sietch_colab/data_share/illex/popgen_data/analysis/steps/03_karyotype/polarize"

aa=pd.read_csv(f"{O}/aa_af.txt",sep="\t",names=["chrom","pos","ref","alt","af_aa"],na_values=".")
bb=pd.read_csv(f"{O}/bb_af.txt",sep="\t",names=["chrom","pos","ref","alt","af_bb"],na_values=".")
df=aa.merge(bb[["pos","af_bb"]],on="pos")
df=df.dropna(subset=["af_aa","af_bb"])
# arrangement-diagnostic: near-fixed opposite alleles
df["ddiff"]=(df.af_aa-df.af_bb).abs()
diag=df[(df.ddiff>0.8)&((df.af_aa<0.1)|(df.af_aa>0.9))&((df.af_bb<0.1)|(df.af_bb>0.9))].copy()
print(f"SNPs in region: {len(df)}  arrangement-diagnostic (|dAF|>0.8, both near-fixed): {len(diag)}")

# each arrangement's fixed base
diag["aa_base"]=np.where(diag.af_aa>0.5, diag.alt, diag.ref)   # AA-hom fixed allele
diag["bb_base"]=np.where(diag.af_bb>0.5, diag.alt, diag.ref)

def tally(anc_map, label):
    a=diag.copy(); a["anc"]=a.pos.map(anc_map)
    a=a[a.anc.isin(["A","C","G","T"])]
    a=a[(a.anc==a.aa_base)|(a.anc==a.bb_base)]     # ancestral must match one arrangement's fixed allele
    aa_anc=(a.anc==a.aa_base).sum(); bb_anc=(a.anc==a.bb_base).sum(); n=len(a)
    print(f"\n[{label}] usable diagnostic sites: {n}")
    print(f"  AA arrangement ancestral: {aa_anc} ({100*aa_anc/n:.1f}%)")
    print(f"  BB arrangement ancestral: {bb_anc} ({100*bb_anc/n:.1f}%)")
    print(f"  => ANCESTRAL (older) arrangement: {'AA' if aa_anc>bb_anc else 'BB'}  (derived/younger: {'BB' if aa_anc>bb_anc else 'AA'})")
    return aa_anc,bb_anc,n

# est-sfs ancestral
anc=pd.read_csv(f"{O}/anc_estsfs.txt",sep="\t",names=["chrom","pos","aa"],na_values=".")
anc_map=dict(zip(anc.pos, anc.aa.astype(str).str.upper()))
tally(anc_map,"est-sfs ancestral")

# coindetii allele (outgroup): GT 0->ref, 1->alt; use if homozygous
try:
    co=pd.read_csv(f"{O}/coindetii.txt",sep="\t",header=None,na_values=".")
    co.columns=["chrom","pos","ref","alt"]+[f"gt{i}" for i in range(co.shape[1]-4)]
    def coin_allele(r):
        g=str(r.gt0)
        if g in ("0/0","0|0"): return r.ref
        if g in ("1/1","1|1"): return r.alt
        return None
    co["cb"]=co.apply(coin_allele,axis=1)
    cmap=dict(zip(co.pos, co.cb))
    tally(cmap,"coindetii outgroup")
except Exception as e:
    print(f"[coindetii cross-check skipped: {e}]")
