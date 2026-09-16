"""Fast, vectorized fixed-difference architecture in the chr2 inversion."""
import numpy as np, pandas as pd
A = "/sietch_colab/data_share/illex/popgen_data"
K = f"{A}/analysis/steps/03_karyotype/polarize"
CDS = f"{A}/mkado_illex/inputs/cds.bed"
GFF = "/sietch_colab/data_share/illex/andy_links/Illex_F24.gene_lnc_pseudo.func.fix.sq3.FINAL.v2.gff3"
OUT = f"{A}/analysis/steps/03_karyotype/inversion_content"
INV_S, INV_E = 60_540_000, 79_500_000
LOG = open(f"{OUT}/fixed_diffs_summary.txt", "w")
def say(*a): print(*a); print(*a, file=LOG); LOG.flush()

aa = pd.read_csv(f"{K}/aa_af.txt", sep="\t", header=None, names=["c","pos","ref","alt","aa"])
bb = pd.read_csv(f"{K}/bb_af.txt", sep="\t", header=None, names=["c","pos","ref","alt","bb"])
m = aa[["pos","aa"]].copy(); m["bb"] = bb["bb"].values          # positions aligned (verified)
m["delta"] = (m.aa - m.bb).abs()
inv = m[(m.pos >= INV_S) & (m.pos <= INV_E)].reset_index(drop=True)
say(f"chr2:60-80Mb SNPs {len(m)}; inversion body (60.54-79.50Mb) {len(inv)}")

# CDS on chr2 -> interval arrays
cds = pd.read_csv(CDS, sep="\t", header=None, names=["c","s","e"])
cds = cds[cds.c.astype(str) == "2"].sort_values("s").reset_index(drop=True)
cs, ce = cds.s.values, cds.e.values
def exonic_mask(pos):
    i = np.searchsorted(cs, pos, side="right") - 1
    ok = i >= 0
    out = np.zeros(len(pos), bool)
    out[ok] = pos[ok] < ce[i[ok]]
    return out
inv_exon = exonic_mask(inv.pos.values)
say(f"exonic fraction of inversion SNPs: {100*inv_exon.mean():.1f}%  (CDS bp coverage baseline)")

for thr in (0.8, 0.9, 0.95, 0.99):
    fx = inv.delta.values >= thr
    ex = fx & inv_exon
    n, ne = int(fx.sum()), int(ex.sum())
    enr = (ne/max(n,1)) / max(inv_exon.mean(),1e-9)
    say(f"|dAF|>={thr}: fixed-diffs={n}  exonic={ne} ({100*ne/max(n,1):.1f}%)  exon-enrichment={enr:.2f}x")

# gene assignment for exonic fixed diffs at 0.9
genes = []
for ln in open(GFF):
    if ln.startswith("#"): continue
    f = ln.split("\t")
    if len(f) < 9 or f[2] != "gene": continue
    if f[0] != "2": continue
    s, e = int(f[3]), int(f[4])
    if e < INV_S or s > INV_E: continue
    gid = next((kv.split("=",1)[1] for kv in f[8].split(";") if kv.startswith("ID=")), "")
    name = next((kv.split("=",1)[1] for kv in f[8].split(";") if kv.startswith("Note=") or kv.startswith("product=")), "")
    genes.append((s, e, gid.replace("gene-",""), name[:50]))
gdf = pd.DataFrame(genes, columns=["s","e","gene","note"]).sort_values("s").reset_index(drop=True)
gs, ge = gdf.s.values, gdf.e.values
FIX = 0.9
fxpos = inv.pos.values[(inv.delta.values >= FIX) & inv_exon]
gi = np.searchsorted(gs, fxpos, side="right") - 1
hit = {}
for p, i in zip(fxpos, gi):
    if i >= 0 and p <= ge[i]:
        hit[gdf.gene.iloc[i]] = hit.get(gdf.gene.iloc[i], 0) + 1
pg = pd.Series(hit).sort_values(ascending=False)
say(f"\n[cluster @|dAF|>={FIX}] {len(fxpos)} exonic fixed diffs in {len(pg)} genes")
if len(pg):
    top = pg.head(12)
    say(f"top genes carry {int(top.sum())} of {len(fxpos)} ({100*top.sum()/len(fxpos):.0f}%):")
    for gene, n in top.items():
        note = gdf[gdf.gene==gene].note.iloc[0]
        say(f"  {gene}  {n:3d}  {note}")
pg.to_csv(f"{OUT}/inversion_fixed_per_gene.tsv", sep="\t", header=["n_exonic_fixed"])
LOG.close()
