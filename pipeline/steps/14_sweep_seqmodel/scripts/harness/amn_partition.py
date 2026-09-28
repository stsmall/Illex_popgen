#!/usr/bin/env python3
"""Partition illex polymorphism in the Amn/Zfest region by feature class + codon position,
to localize the chrZ sweep target (Amn coding protein vs antisense lncRNA vs regulatory/intron)."""
import sys
import numpy as np

GFF = "/sietch_colab/data_share/illex/andy_links/Illex_F24.gene_lnc_pseudo.func.fix.sq3.FINAL.v2.gff3"
SNPS = sys.argv[1]
AMN = "LOC_00003844"; LNC1 = "XLOC_014991U_lncRNA"; LNC2 = "lncrna-gene-8"
REG0, REG1 = 21460000, 21560000   # focused Amn/Zfest window

def attr(s, key):
    for f in s.split(";"):
        if f.strip().startswith(key + "="): return f.strip()[len(key)+1:]
    return ""

cds=[]; utr=[]; amn_exon=[]; lnc_exon=[]; amn_mrna=set()
for line in open(GFF):
    if line[0]=="#": continue
    f=line.rstrip("\n").split("\t")
    if f[0]!="Z": continue
    s,e=int(f[3]),int(f[4])
    if e<REG0 or s>REG1: continue
    typ=f[2]; par=attr(f[8],"Parent"); ident=attr(f[8],"ID")
    if typ=="mRNA" and attr(f[8],"Parent")==AMN or (typ=="mRNA" and AMN in ident): amn_mrna.add(ident)
# 2nd pass now that we know Amn mRNA ids (also accept Parent chain by AMN gene)
for line in open(GFF):
    if line[0]=="#": continue
    f=line.rstrip("\n").split("\t")
    if f[0]!="Z": continue
    s,e=int(f[3]),int(f[4]); typ=f[2]; par=attr(f[8],"Parent"); strand=f[6]; phase=f[7]
    if e<REG0 or s>REG1: continue
    is_amn = (par in amn_mrna) or (AMN in par) or (AMN in attr(f[8],"ID"))
    is_lnc = (LNC1 in par or LNC2 in par or LNC1 in attr(f[8],"ID") or LNC2 in attr(f[8],"ID"))
    if typ=="CDS" and is_amn: cds.append((s,e,strand,int(phase) if phase.isdigit() else 0))
    elif typ in ("five_prime_UTR","three_prime_UTR") and is_amn: utr.append((s,e))
    elif typ=="exon" and is_amn: amn_exon.append((s,e))
    elif typ=="exon" and is_lnc: lnc_exon.append((s,e))

# base-level codon position for Amn CDS (+ strand: order ascending, running index from 0)
cds.sort()
codonpos={}  # genomic pos -> 1/2/3
idx=0
for s,e,strand,ph in cds:
    for p in range(s,e+1):
        codonpos[p]= (idx%3)+1; idx+=1
def spanset(ivs):
    ss=set()
    for s,e in ivs:
        ss.update(range(s,e+1))
    return ss
utr_s=spanset(utr); amn_exon_s=spanset(amn_exon); lnc_exon_s=spanset(lnc_exon)
cds_s=set(codonpos)
# Amn gene span for intron
amn_span=set(range(21471373,21541796+1))
intron_s=amn_span - amn_exon_s   # intronic = in gene, not in any Amn exon

# category base counts (bp)
bp={ "CDS_pos1":sum(1 for v in codonpos.values() if v==1),
     "CDS_pos2":sum(1 for v in codonpos.values() if v==2),
     "CDS_pos3":sum(1 for v in codonpos.values() if v==3),
     "Amn_5UTR":len(utr_s - cds_s),
     "Amn_intron":len(intron_s),
     "lncRNA_exon":len(lnc_exon_s) }
# classify SNPs (priority: CDS pos -> UTR -> lncRNA_exon -> intron)
snp=[int(x) for x in open(SNPS) if x.strip()]
cnt={k:0 for k in bp}; cnt["lncRNA_exon_antisense_in_Amn"]=0
for p in snp:
    if p in codonpos: cnt["CDS_pos%d"%codonpos[p]]+=1
    elif p in utr_s: cnt["Amn_5UTR"]+=1
    elif p in lnc_exon_s: cnt["lncRNA_exon"]+=1
    elif p in intron_s: cnt["Amn_intron"]+=1

print(f"{'category':<14}{'SNPs':>6}{'bp':>7}{'SNP/kb':>9}")
order=["CDS_pos1","CDS_pos2","CDS_pos3","Amn_5UTR","Amn_intron","lncRNA_exon"]
for k in order:
    b=bp.get(k,0); c=cnt.get(k,0)
    print(f"{k:<14}{c:>6}{b:>7}{(1000*c/b if b else 0):>9.2f}")
# note lncRNA overlaps Amn 3' region (antisense); how many lncRNA-exon SNPs are also in Amn intron/3'?
