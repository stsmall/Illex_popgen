#!/usr/bin/env python
"""Pair Illex CDS with the lifted Sthenoteuthis CDS of the same transcript (one FASTA per gene) for
codon alignment. Keeps one transcript per gene (the lifted one), Illex CDS length a multiple of 3.
Usage: make_pairs.py ILLEX_CDS STHENO_CDS LOC_TSV OUT_DIR"""
import sys, os
ill, sth, loc, out = sys.argv[1:5]; os.makedirs(out, exist_ok=True)
def read(f):
    d, k = {}, None
    for l in open(f):
        if l.startswith(">"): k = l[1:].split()[0]; d[k] = []
        else: d[k].append(l.strip())
    return {k: "".join(v).upper() for k, v in d.items()}
I, S = read(ill), read(sth)
L = {l.split()[0]: l.split()[1:] for l in open(loc)}
seen = set(); n = 0
with open(f"{out}/pairs.tsv", "w") as o:
    o.write("tx\tgene\tchrom\tstart\tend\till_len\tsth_len\n")
    for tx, s in S.items():
        g = tx.rsplit("-mRNA", 1)[0]
        if g in seen or tx not in I or tx not in L or len(I[tx]) % 3 or len(s) < 150: continue
        seen.add(g); n += 1
        open(f"{out}/{tx}.ill.fa", "w").write(f">illex\n{I[tx]}\n")
        open(f"{out}/{tx}.sth.fa", "w").write(f">stheno\n{s}\n")
        o.write(f"{tx}\t{g}\t{L[tx][0]}\t{L[tx][1]}\t{L[tx][2]}\t{len(I[tx])}\t{len(s)}\n")
print(n, "pairs")
