#!/usr/bin/env python3
"""Supplementary Table S5: functional clusters of the GO terms enriched (FDR < 0.05) in the
background-selection-robust Tier-1 sweep genes.

Input is the masked BGS-robust enrichment table. Every enriched term, in all three ontologies, is
assigned to one cluster by the first keyword rule it matches (rules are applied in the order below),
so the cluster counts sum to the number of enriched terms. Prints the per-cluster summary and the
LaTeX rows, and writes the term-to-cluster assignment next to this script's output table.
"""
import csv
import re
import sys
from pathlib import Path

ENRICH = Path("/sietch_colab/data_share/illex/popgen_data/analysis/steps/14_sweep_seqmodel/results/"
              "empirical_scan_fullsfs/hmm_decode_masked/tier1_BGS_GO_enrichment.tsv")
OUT = Path(__file__).resolve().parent.parent / "data" / "tableS5_go_clusters.tsv"
FDR_MAX = 0.05

# (cluster, regex) in priority order. The first match wins.
RULES = [
    ("Photoreception / sensory response",
     r"light|radiation|opsin|taxis|response to external stimulus"),
    ("Neural development",
     r"neur|nervous|dendrit|axon|synap|learning|memory|cognition|vesicle uncoating|clathrin coat"),
    ("Cell morphogenesis \\& differentiation",
     r"(?<!germ )\bcell (morphogenesis|differentiation|development|projection|part)|cellular developmental"
     r"|cellular anatomical entity morphogenesis|regulation of cell (differentiation|development)"),
    ("RNA / transcription regulation",
     r"RNA|transcription|DNA[- ]binding|DNA binding|nucleobase|nucleic acid|gene expression"
     r"|cis-regulatory|biosynthetic process"),
    ("Organ / anatomical development",
     r"\borgan\b|system development|tissue|embryo|epitheli|anatomical structure|developmental process"
     r"|multicellular organism development|multicellular organismal development|reproduct|oogenesis|germ cell"),
    ("Other regulation / metabolism", r"."),
]

SHORT = {"regulation": "reg.\\", "positive regulation": "pos.\\ reg.\\", "negative regulation": "neg.\\ reg.\\"}


def short(term):
    for k in ("positive regulation", "negative regulation", "regulation"):
        if term.startswith(k + " of "):
            return SHORT[k] + " " + term[len(k) + 4:]
    return term


def sci(x):
    m, e = f"{x:.0e}".split("e")
    return f"${int(m)}\\times10^{{{int(e)}}}$"


def main():
    rows = [r for r in csv.DictReader(open(ENRICH), delimiter="\t") if float(r["FDR"]) < FDR_MAX]
    clusters = {name: [] for name, _ in RULES}
    for r in rows:
        for name, rx in RULES:
            if re.search(rx, r["term"]):
                clusters[name].append(r)
                r["cluster"] = name
                break
    assert sum(len(v) for v in clusters.values()) == len(rows)
    with open(OUT, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["go_id", "term", "category", "fold", "FDR", "cluster"])
        for r in sorted(rows, key=lambda r: (r["cluster"], float(r["FDR"]))):
            w.writerow([r["go_id"], r["term"], r["category"], r["fold"], r["FDR"], r["cluster"]])
    cats = {}
    for r in rows:
        cats[r["category"]] = cats.get(r["category"], 0) + 1
    print(f"{len(rows)} terms", cats, file=sys.stderr)
    order = [n for n, _ in RULES if n != "Other regulation / metabolism"]
    order.sort(key=lambda n: -len(clusters[n]))
    order.append("Other regulation / metabolism")
    for n in order:
        v = clusters[n]
        folds = [float(r["fold"]) for r in v]
        best = min(v, key=lambda r: float(r["FDR"]))
        print(f"    {n} & {len(v)} & {min(folds):.1f}--{max(folds):.1f} & "
              f"{short(best['term'])} ({sci(float(best['FDR']))}) \\\\")


if __name__ == "__main__":
    main()
