"""Tier-1 gene-ontology enrichment (reconstructed; reproduces the published run exactly).
Background: genes with EnTAP GO terms that are in the annotation GFF, chrZ excluded. Tier-1 genes: genes in
tier1_genes.tsv (overlapping a Tier-1 call; nearest-gene assignments excluded) that are in that background, optionally restricted
to a set of calls. One-sided Fisher's exact test per term with >=1 Tier-1 gene and >=3 background genes; Benjamini-Hochberg FDR.
Usage: tier1_go.py CALLS.tsv OUT.tsv [MASK.bed]   (MASK removes background genes inside masked regions)"""
import sys, pandas as pd, numpy as np
from scipy.stats import fisher_exact
from statsmodels.stats.multitest import multipletests
R="/sietch_colab/data_share/illex/popgen_data/analysis/steps/14_sweep_seqmodel/results/empirical_scan_fullsfs/hmm_decode"
E="/sietch_colab/data_share/illex/annotations/gene_annotation_dir/project/functional_entap/entap_outfiles/final_results/annotated_gene_ontology_terms.tsv"
GFF="/sietch_colab/data_share/illex/andy_links/Illex_F24.gene_lnc_pseudo.func.fix.sq3.FINAL.v2.gff3"
calls=pd.read_csv(sys.argv[1],sep="\t"); calls["chrom"]=calls.chrom.astype(str)
e=pd.read_csv(E,sep="\t"); e["gene"]=e.query_sequence.str.replace(r"-mRNA-\d+$","",regex=True); e=e.drop_duplicates(["gene","go_id"])
g=[]
for ln in open(GFF):
    if ln[0]=="#": continue
    f=ln.split("\t")
    if len(f)>8 and f[2]=="gene":
        gid=[kv.split("=",1)[1] for kv in f[8].split(";") if kv.startswith("ID=")]
        if gid: g.append((f[0],int(f[3]),int(f[4]),gid[0].replace("gene-","")))
G=pd.DataFrame(g,columns=["chrom","s","e","gene"]); G["chrom"]=G.chrom.astype(str); G=G[G.chrom!="Z"]
if len(sys.argv)>3:
    m=pd.read_csv(sys.argv[3],sep="\t",header=None,usecols=[0,1,2],names=["c","s","e"]); m["c"]=m.c.astype(str)
    inm=np.zeros(len(G),bool)
    for r in m.itertuples(): inm|=(G.chrom.values==r.c)&(G.s.values<r.e)&(G.e.values>r.s)
    G=G[~inm]
bg=set(G.gene)&set(e.gene)
t=pd.read_csv(f"{R}/tier1_genes.tsv",sep="\t"); t["chrom"]=t.chrom.astype(str)
t=t[t.overlap=="overlap"]; t=t.merge(calls[["chrom","start"]].rename(columns={"start":"call_start"}),on=["chrom","call_start"])
T=set(t.gene)&bg
e=e[e.gene.isin(bg)]; go2=e.groupby("go_id").gene.apply(set); name=e.drop_duplicates("go_id").set_index("go_id")
rows=[]; NT,NB=len(T),len(bg)
for go,gs in go2.items():
    a=len(gs&T); b=len(gs)
    if a<1 or b<3: continue
    _,p=fisher_exact([[a,NT-a],[b-a,NB-NT-(b-a)]],alternative="greater")
    rows.append(dict(go_id=go,term=name.loc[go,"go_term"],category=name.loc[go,"category"],n_tier1=a,n_background=b,N_tier1=NT,N_background=NB,fold=round(a/NT/(b/NB),3),p=p))
d=pd.DataFrame(rows); d["FDR"]=multipletests(d.p,method="fdr_bh")[1]; d=d.sort_values("p")
d.to_csv(sys.argv[2],sep="\t",index=False)
print(f"N_tier1={NT} N_background={NB} terms FDR<0.05: {(d.FDR<0.05).sum()} (BP {((d.FDR<0.05)&(d.category=='biological_process')).sum()})")
for go in ["GO:0009416","GO:0009314","GO:0048812","GO:0000902","GO:0048699","GO:0048666","GO:0006366","GO:0009583","GO:0007602"]:
    x=d[d.go_id==go]
    if len(x): print(f"  {x.term.iloc[0]}: fold {x.fold.iloc[0]:.2f} FDR {x.FDR.iloc[0]:.1e} (n={x.n_tier1.iloc[0]})")
