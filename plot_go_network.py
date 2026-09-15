"""Supp: GO enrichment-map network for the Tier-1 sweep genes. Nodes = enriched
biological-process terms (size = # Tier-1 genes, colour = -log10 FDR); edges connect
terms that share Tier-1 genes (Jaccard >= 0.15). Reveals the two connected modules:
neural development and RNA-regulatory/transcriptional control. Run under bioinfo-buddy (networkx)."""
import os
os.environ.setdefault("MPLCONFIGDIR", "/dev/shm/mplcache")
import sys
sys.path.insert(0, "/sietch_colab/data_share/illex/popgen_data/analysis/manuscript")
from figstyle import apply, C, despine
apply()
import numpy as np, pandas as pd, matplotlib.pyplot as plt
import matplotlib.patheffects as pe
import networkx as nx

D = "/sietch_colab/data_share/illex/popgen_data/analysis/steps/14_sweep_seqmodel/results/empirical_scan_fullsfs/hmm_decode"
OUT = "/sietch_colab/data_share/illex/popgen_data/analysis/manuscript/supp_figures/supp_go_network.png"

enr = pd.read_csv(f"{D}/tier1_GO_enrichment.tsv", sep="\t")
enr = enr[(enr.FDR < 0.05) & (enr.category == "biological_process")].copy()
enr = enr.sort_values("FDR").head(15)                       # top terms for legibility
genes = pd.read_csv(f"{D}/tier1_genes.tsv", sep="\t")
genes = genes[genes.annotated == True] if genes.annotated.dtype == bool else genes
gbp = genes.dropna(subset=["GO_biological_process"])

# gene set per enriched term (genes whose GO_BP string contains the term name)
def gene_set(term):
    t = term.lower()
    m = gbp["GO_biological_process"].str.lower().str.contains(t, regex=False, na=False)
    return set(gbp.loc[m, "gene"])
sets = {r.go_id: gene_set(r.term) for _, r in enr.iterrows()}
enr = enr[enr.go_id.isin([k for k, v in sets.items() if len(v) >= 3])]

G = nx.Graph()
for _, r in enr.iterrows():
    G.add_node(r.go_id, term=r.term, n=len(sets[r.go_id]), fdr=r.FDR)
ids = list(G.nodes)
for i in range(len(ids)):
    for j in range(i + 1, len(ids)):
        a, b = sets[ids[i]], sets[ids[j]]
        jac = len(a & b) / len(a | b) if a | b else 0
        if jac >= 0.15:
            G.add_edge(ids[i], ids[j], w=jac)

pos = nx.spring_layout(G, k=2.2, seed=3, weight="w", iterations=400)
sizes = [80 + G.nodes[n]["n"] / max(G.nodes[m]["n"] for m in G) * 900 for n in G]
colors = [-np.log10(G.nodes[n]["fdr"]) for n in G]

fig, ax = plt.subplots(figsize=(9.6, 8.2))
nx.draw_networkx_edges(G, pos, ax=ax, edge_color="#CCCCCC",
                       width=[G[u][v]["w"] * 4 for u, v in G.edges], alpha=0.7)
nodes = nx.draw_networkx_nodes(G, pos, ax=ax, node_size=sizes, node_color=colors,
                               cmap="viridis", edgecolors="white", linewidths=0.8)
import textwrap
for n in G:
    t = G.nodes[n]["term"]
    lab = "\n".join(textwrap.wrap(t, width=15))     # full name, wrapped (no truncation)
    txt = ax.annotate(lab, pos[n], fontsize=6.6, ha="center", va="center", linespacing=0.95,
                      xytext=(0, 0), textcoords="offset points", color="black", zorder=10)
    txt.set_path_effects([pe.withStroke(linewidth=2.8, foreground="white")])
ax.set_axis_off()
cb = fig.colorbar(nodes, ax=ax, fraction=0.04, pad=0.01, shrink=0.7)
cb.set_label(r"$-\log_{10}$ FDR", fontsize=8.5); cb.ax.tick_params(labelsize=7)
ax.set_title("Tier-1 sweep genes: GO enrichment map (shared-gene network)",
             loc="left", fontweight="bold", fontsize=11)
fig.savefig(OUT, dpi=220, bbox_inches="tight")
print("wrote", OUT, "|", G.number_of_nodes(), "nodes", G.number_of_edges(), "edges")
