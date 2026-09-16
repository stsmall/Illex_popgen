"""Characterise the chr2 inversion's content, independent of any sweep signal:
  (1) genes captured by the inversion body + their GO enrichment vs the genome,
  (2) the architecture of fixed differences between arrangements (AA vs BB) -- how
      many, how many fall in exons/CDS, and whether they cluster in a few genes.
Inversion body: chr2:60.54-79.50 Mb (~19 Mb). AA/BB per-SNP allele freqs from step 03.
"""
import os
import numpy as np, pandas as pd
from scipy.stats import fisher_exact
from statsmodels.stats.multitest import multipletests

A = "/sietch_colab/data_share/illex/popgen_data"
GFF = "/sietch_colab/data_share/illex/andy_links/Illex_F24.gene_lnc_pseudo.func.fix.sq3.FINAL.v2.gff3"
GENEGO = f"{A}/analysis/steps/17_ensemble_scan/results/fst_persnp/gene_go.tsv"
GOTERMS = "/sietch_colab/data_share/illex/annotations/gene_annotation_dir/project/functional_entap/entap_outfiles/final_results/annotated_gene_ontology_terms.tsv"
CDS = f"{A}/mkado_illex/inputs/cds.bed"
AA = f"{A}/analysis/steps/03_karyotype/polarize/aa_af.txt"
BB = f"{A}/analysis/steps/03_karyotype/polarize/bb_af.txt"
OUT = f"{A}/analysis/steps/03_karyotype/inversion_content"
os.makedirs(OUT, exist_ok=True)
INV_S, INV_E = 60_540_000, 79_500_000

# ---------- genes in the inversion ----------
genes = []
for ln in open(GFF):
    if ln.startswith("#"):
        continue
    f = ln.rstrip("\n").split("\t")
    if len(f) < 9 or f[2] != "gene":
        continue
    chrom, s, e = f[0], int(f[3]), int(f[4])
    gid = ""
    for kv in f[8].split(";"):
        if kv.startswith("ID=") or kv.startswith("geneID="):
            gid = kv.split("=", 1)[1]
    genes.append((chrom, s, e, gid))
g = pd.DataFrame(genes, columns=["chrom", "start", "end", "gene"])
g["chrom"] = g["chrom"].astype(str)
inv_genes = g[(g.chrom == "2") & (g.start < INV_E) & (g.end > INV_S)]["gene"].tolist()
inv_genes = sorted(set(x.replace("gene-", "") for x in inv_genes if x))
print(f"[genes] inversion body carries {len(inv_genes)} genes; total annotated genes {g.gene.nunique()}")

# ---------- GO enrichment (BP), inversion genes vs genome background ----------
gg = pd.read_csv(GENEGO, sep="\t", header=None, names=["gene", "go"])
gt = pd.read_csv(GOTERMS, sep="\t")
gt["gene"] = gt["query_sequence"].str.replace(r"-mRNA-\d+$", "", regex=True)
id2name = dict(zip(gt.go_id, gt.go_term))
id2cat = dict(zip(gt.go_id, gt.category))
bp_ids = {i for i, c in id2cat.items() if c == "biological_process"}

all_genes = set(gg.gene.unique())
inv_set = set(inv_genes) & all_genes
N_inv, N_bg = len(inv_set), len(all_genes)
go2genes = gg.groupby("go")["gene"].apply(set)
rows = []
for go, gene_set in go2genes.items():
    if go not in bp_ids:
        continue
    a = len(gene_set & inv_set)
    if a < 2:
        continue
    b = N_inv - a
    c = len(gene_set) - a
    d = N_bg - N_inv - c
    OR, p = fisher_exact([[a, b], [c, d]], alternative="greater")
    fold = (a / N_inv) / (len(gene_set) / N_bg)
    rows.append((go, id2name.get(go, go), a, len(gene_set), fold, OR, p))
enr = pd.DataFrame(rows, columns=["go_id", "term", "n_inv", "n_bg", "fold", "odds_ratio", "p"])
enr["FDR"] = multipletests(enr.p, method="fdr_bh")[1] if len(enr) else []
enr = enr.sort_values("p")
enr.to_csv(f"{OUT}/inversion_GO_enrichment.tsv", sep="\t", index=False)
sig = enr[enr.FDR < 0.05]
print(f"[GO] {N_inv} inversion genes with GO; {len(sig)} BP terms FDR<0.05, {(enr.p<0.05).sum()} nominal")
print(sig.head(12)[["term", "n_inv", "fold", "FDR"]].to_string(index=False))

# ---------- fixed differences between arrangements ----------
aa = pd.read_csv(AA, sep="\t", header=None, names=["chrom", "pos", "ref", "alt", "aa"])
bb = pd.read_csv(BB, sep="\t", header=None, names=["chrom", "pos", "ref", "alt", "bb"])
m = aa.merge(bb[["pos", "bb"]], on="pos")
m["delta"] = (m.aa - m.bb).abs()
inv = m[(m.pos >= INV_S) & (m.pos <= INV_E)].copy()
print(f"\n[fixed] chr2 SNPs {len(m)}; in inversion body {len(inv)}")
for thr in (0.8, 0.9, 0.95, 0.99):
    print(f"  |dAF|>={thr}: whole-chr2={int((m.delta>=thr).sum())}  inversion={int((inv.delta>=thr).sum())}")

# exonic vs non-exonic among near-fixed differences in the inversion
cds = pd.read_csv(CDS, sep="\t", header=None, names=["chrom", "s", "e"])
cds = cds[cds.chrom.astype(str) == "2"].sort_values("s")
starts, ends = cds.s.values, cds.e.values


def in_cds(pos):
    i = np.searchsorted(starts, pos, side="right") - 1
    return i >= 0 and pos < ends[i]


FIX = 0.9
fx = inv[inv.delta >= FIX].copy()
fx["exonic"] = fx.pos.apply(in_cds)
n_ex, n_all = int(fx.exonic.sum()), len(fx)
# background exon fraction of the inversion body (bp in CDS / total)
inv_cds_bp = int(np.clip(ends, INV_S, INV_E).clip(min=INV_S).sum() - np.clip(starts, INV_S, INV_E).clip(max=INV_E).sum())
exon_frac = inv_cds_bp / (INV_E - INV_S)
print(f"\n[fixed>={FIX}] in inversion: {n_all} fixed diffs, {n_ex} exonic ({100*n_ex/max(n_all,1):.1f}%); "
      f"inversion CDS covers {100*exon_frac:.1f}% of bp")
# enrichment of fixed diffs in exons vs expectation
all_ex = inv.pos.apply(in_cds)
base_ex = all_ex.mean()
print(f"  exon fraction of ALL inversion SNPs: {100*base_ex:.1f}%  -> fixed-diff exon enrichment "
      f"{ (n_ex/max(n_all,1)) / max(base_ex,1e-9):.2f}x")

# clustering by gene: assign each exonic fixed diff to its gene
gi = g[(g.chrom == "2") & (g.start < INV_E) & (g.end > INV_S)].copy()
gi["gene"] = gi.gene.str.replace("gene-", "")
def gene_of(pos):
    h = gi[(gi.start <= pos) & (gi.end >= pos)]
    return h.gene.iloc[0] if len(h) else None
fxe = fx[fx.exonic].copy()
fxe["gene"] = fxe.pos.apply(gene_of)
per_gene = fxe.dropna(subset=["gene"]).groupby("gene").size().sort_values(ascending=False)
ngenes_hit = per_gene.size
top = per_gene.head(10)
print(f"\n[cluster] {n_ex} exonic fixed diffs fall in {ngenes_hit} genes; "
      f"top {min(10,ngenes_hit)} genes carry {int(top.sum())} ({100*top.sum()/max(n_ex,1):.0f}%):")
print(top.to_string())
fxe.to_csv(f"{OUT}/inversion_exonic_fixed_diffs.tsv", sep="\t", index=False)
per_gene.to_csv(f"{OUT}/inversion_fixed_per_gene.tsv", sep="\t", header=["n_fixed_exonic"])
print("\nwrote", OUT)
