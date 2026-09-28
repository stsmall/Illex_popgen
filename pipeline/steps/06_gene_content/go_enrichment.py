#!/usr/bin/env python
# GO enrichment of chr2:60-80Mb inversion genes (102) vs genome background,
# using EnTAP GO terms. Hypergeometric test, BH-corrected.
import pandas as pd, numpy as np
from scipy.stats import hypergeom
from statsmodels.stats.multitest import multipletests

D="/sietch_colab/data_share/illex/popgen_data/analysis/steps/06_gene_content"
GO=("/sietch_colab/data_share/illex/annotations/gene_annotation_dir/project/"
    "functional_entap/entap_outfiles/final_results/annotated_without_contam_gene_ontology_terms.tsv")

def to_gene(q): return str(q).split("-mRNA")[0].split("-RA")[0]
go=pd.read_csv(GO,sep="\t")
go["gene"]=go["query_sequence"].map(to_gene)
# drop root terms (uninformative)
roots={"GO:0005575","GO:0003674","GO:0008150"}
go=go[~go["go_id"].isin(roots)]
gene2go=go.groupby("gene")["go_id"].apply(set).to_dict()
term_name=dict(zip(go["go_id"],go["go_term"])); term_cat=dict(zip(go["go_id"],go["category"]))
bg_genes=set(gene2go)                                  # background = genes with >=1 GO
N=len(bg_genes)
inv=set(pd.read_csv(f"{D}/inv_genes.tsv",sep="\t",header=None)[0])
inv_bg=inv & bg_genes
K=len(inv_bg)
print(f"background genes w/ GO: {N}; inversion genes: {len(inv)} (w/ GO: {K})", flush=True)

# term -> gene sets
from collections import defaultdict
term_genes=defaultdict(set)
for g,ts in gene2go.items():
    for t in ts: term_genes[t].add(g)

rows=[]
for t,genes in term_genes.items():
    n=len(genes); k=len(genes & inv_bg)
    if k<2 or n<3: continue                            # require >=2 inversion genes, term not ultra-rare
    p=hypergeom.sf(k-1,N,n,K)                           # P(X>=k)
    rows.append((t,term_cat.get(t,""),term_name.get(t,""),k,n,K*n/N,p))
res=pd.DataFrame(rows,columns=["go_id","cat","term","inv_k","bg_n","expected","p"])
res["fdr"]=multipletests(res["p"],method="fdr_bh")[1]
res=res.sort_values("p")
res.to_csv(f"{D}/go_enrichment.tsv",sep="\t",index=False)
print("\n=== top enriched GO terms (inversion vs background) ===", flush=True)
print(res.head(20)[["cat","term","inv_k","bg_n","expected","p","fdr"]].to_string(index=False), flush=True)
sig=res[res.fdr<0.1]
print(f"\nFDR<0.1: {len(sig)} terms", flush=True)
