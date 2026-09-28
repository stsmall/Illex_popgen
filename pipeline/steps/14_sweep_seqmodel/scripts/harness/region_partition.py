#!/usr/bin/env python3
"""Partition illex polymorphism in a region by feature class + codon position (strand-aware),
to localize a sweep: CDS pos1/2/3, UTR, intron (focal gene), and INTERGENIC.
Args: snps_file gene_id gene_start gene_end strand region_start region_end"""
import sys
SNPS, GID, GS, GE, STR, R0, R1 = sys.argv[1], sys.argv[2], int(sys.argv[3]), int(sys.argv[4]), sys.argv[5], int(sys.argv[6]), int(sys.argv[7])
GFF = "/sietch_colab/data_share/illex/andy_links/Illex_F24.gene_lnc_pseudo.func.fix.sq3.FINAL.v2.gff3"
def attr(s,k):
    for f in s.split(";"):
        if f.strip().startswith(k+"="): return f.strip()[len(k)+1:]
    return ""
cds=[]; exon=[]; utr=[]; gene_spans=[]; focal_mrna=set()
for line in open(GFF):
    if line[0]=="#": continue
    f=line.rstrip("\n").split("\t")
    if f[0]!="Z": continue
    s,e=int(f[3]),int(f[4])
    if e<R0 or s>R1: continue
    if f[2]=="gene": gene_spans.append((s,e))
    if f[2]=="mRNA" and (GID in attr(f[8],"Parent") or GID in attr(f[8],"ID")): focal_mrna.add(attr(f[8],"ID"))
for line in open(GFF):
    if line[0]=="#": continue
    f=line.rstrip("\n").split("\t")
    if f[0]!="Z": continue
    s,e=int(f[3]),int(f[4])
    if e<R0 or s>R1: continue
    isf = (attr(f[8],"Parent") in focal_mrna) or (GID in attr(f[8],"Parent")) or (GID in attr(f[8],"ID"))
    if not isf: continue
    if f[2]=="CDS": cds.append((s,e))
    elif f[2]=="exon": exon.append((s,e))
    elif f[2] in ("five_prime_UTR","three_prime_UTR"): utr.append((s,e))
def sset(ivs):
    ss=set()
    for s,e in ivs: ss.update(range(s,e+1))
    return ss
# codon position, strand-aware; use unique CDS bases in transcript order
cds_bases=sorted(sset(cds))
if STR=="-": cds_bases=cds_bases[::-1]
codonpos={p:(i%3)+1 for i,p in enumerate(cds_bases)}
utr_s=sset(utr); exon_s=sset(exon); gene_s=set(range(GS,GE+1))
intron_s=gene_s-exon_s
allgene_s=sset(gene_spans)
region_s=set(range(R0,R1+1))
intergenic_s=region_s-allgene_s
bp={"CDS_pos1":sum(v==1 for v in codonpos.values()),"CDS_pos2":sum(v==2 for v in codonpos.values()),
    "CDS_pos3":sum(v==3 for v in codonpos.values()),"UTR":len(utr_s-set(codonpos)),
    "intron":len(intron_s),"intergenic":len(intergenic_s)}
cnt={k:0 for k in bp}
for line in open(SNPS):
    p=int(line);
    if p in codonpos: cnt["CDS_pos%d"%codonpos[p]]+=1
    elif p in utr_s: cnt["UTR"]+=1
    elif p in intron_s: cnt["intron"]+=1
    elif p in intergenic_s: cnt["intergenic"]+=1
print(f"{'category':<12}{'SNPs':>6}{'bp':>8}{'SNP/kb':>9}")
for k in ["CDS_pos1","CDS_pos2","CDS_pos3","UTR","intron","intergenic"]:
    b=bp[k]; print(f"{k:<12}{cnt[k]:>6}{b:>8}{(1000*cnt[k]/b if b else 0):>9.2f}")
