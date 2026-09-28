#!/usr/bin/env python
# SF3B4 substitution inspection: illecebrosus vs coindetii (AnchorWave) codon comparison
# -> Dn/Ds + nonsyn positions; map onto RRM core vs C-term linker via OG0006508 conservation;
# polarize (is the illecebrosus residue derived among cephalopod orthologs?).
import gzip, numpy as np
from collections import Counter

D="/sietch_colab/data_share/illex/popgen_data/analysis/steps/07_sf3b4"
GFF="/sietch_colab/data_share/illex/popgen_data/degenotate_illex/Illex_F24.gene_lnc_pseudo.func.fix.sq3.FINAL.v2.fixID.gff3"
MSA="/sietch_colab/data_share/illex/alignments/Orthofinder/OrthoFinder/Results_Feb27/MultipleSequenceAlignments/OG0006508.fa"
CODON={}
bases="TCAG"
aas="FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
for i,aa in enumerate(aas): CODON[bases[i>>4]+bases[(i>>2)&3]+bases[i&3]]=aa
def tr(c): return CODON.get(c.upper(),"X")

# illecebrosus CDS (from gffread, + strand, in transcript order)
seq=""
for l in open(f"{D}/sf3b4_illex.cds.fa"):
    if not l.startswith(">"): seq+=l.strip()
illex=seq.upper()
# CDS genomic positions in transcript order (+ strand -> ascending exons)
exons=[]
for l in open(GFF):
    f=l.split("\t")
    if len(f)>4 and f[0]=="24" and f[2]=="CDS" and "LOC_00013210" in l: exons.append((int(f[3]),int(f[4])))
exons.sort()
positions=[p for s,e in exons for p in range(s,e+1)]
assert len(positions)==len(illex), f"{len(positions)} vs {len(illex)}"

# coindetii base per chr24 position
coin={}
with gzip.open(f"{D}/maf_bases/24.out2_bases.tsv.gz","rt") as fh:
    for line in fh:
        c,p,b=line.split("\t");
        p=int(p)
        if p>=exons[0][0] and p<=exons[-1][1]: coin[p]=b.strip().upper()
coin_cds="".join(coin.get(p,"N") for p in positions)

# codon compare
Dn=Ds=0; nonsyn=[]  # (aa_pos 1-based, illex_aa, coin_aa)
for i in range(0,len(illex)-2,3):
    ic=illex[i:i+3]; cc=coin_cds[i:i+3]; aap=i//3+1
    if "N" in cc or "-" in cc or len(cc)<3: continue
    ia,ca=tr(ic),tr(cc)
    if ia=="X" or ca=="X" or ia=="*" or ca=="*": continue
    if ic!=cc:
        if ia==ca: Ds+=1
        else: Dn+=1; nonsyn.append((aap,ia,ca))
print(f"SF3B4 illecebrosus-vs-coindetii (AnchorWave): Dn={Dn} Ds={Ds}  (MK reported Dn29/Ds9)")

# OG MSA -> per-illecebrosus-residue conservation (functional core vs variable linker) + polarization
msa={}; name=None
for l in open(MSA):
    if l.startswith(">"): name=l[1:].strip().split()[0]; msa[name]=""
    else: msa[name]+=l.strip()
illkey=[k for k in msa if "LOC_00013210" in k][0]
illrow=msa[illkey]
others=[k for k in msa if k!=illkey]
# map ungapped illex residue index -> MSA column
col_of=[]; ung=0
for col,ch in enumerate(illrow):
    if ch!="-": col_of.append(col)  # col_of[residue_i-1]=column
def conservation(col):
    aa=[msa[o][col] for o in others if msa[o][col]!="-"]
    if not aa: return 0.0,None
    m,cnt=Counter(aa).most_common(1)[0]
    return cnt/len(aa), m
core=linker=0; pol=0
rows=[]
for (aap,ia,ca) in nonsyn:
    if aap-1>=len(col_of): continue
    col=col_of[aap-1]; cons,modal=conservation(col)
    region="CORE(RRM?)" if cons>=0.6 else "linker/variable"
    if cons>=0.6: core+=1
    else: linker+=1
    derived = (modal is not None and ia!=modal)   # illex residue differs from ortholog consensus -> illex-derived
    if derived: pol+=1
    rows.append((aap,ia,ca,round(cons,2),modal,region,"illex-derived" if derived else "shared/ambig"))
print(f"\nnonsyn substitutions in conserved core (cons>=0.6): {core}  vs variable/linker: {linker}")
print(f"polarized as illecebrosus-derived (illex aa != cephalopod ortholog consensus): {pol}/{len(nonsyn)}")
print("\npos illexAA coinAA  cons modalOrtho  region  polarity")
for r in sorted(rows): print(f"{r[0]:>4} {r[1]}    {r[2]}     {r[3]:.2f}  {r[4]}     {r[5]:<16} {r[6]}")
