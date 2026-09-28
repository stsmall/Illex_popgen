#!/usr/bin/env python3
"""Full-coverage sex-averaged recombination-map BED for bgspy, from ReLERNN PREDICT.

ReLERNN gives per-15kb-window recombRate in Morgans/bp (~2e-9), but with gaps (windows
dropped by the nSites filter) and chrom names as bytestring reprs (b'1'). bgspy's recmap
loader wants contiguous coverage. We average male+female per window, then emit contiguous
15kb windows 0..L per chrom, filling gaps with that chrom's mean rate.

Output BED: chrom<TAB>start<TAB>end<TAB>rate(M/bp)  -> use `bgspy calcb --conv-factor 1.0`.
Usage: build_recmap.py <male_PREDICT> <female_PREDICT> <seqlens.tsv> <out.bed>
"""
import sys

WIN = 15000
male_f, fem_f, seqlens_f, out_f = sys.argv[1:5]


def parse(path):
    d = {}
    with open(path) as fh:
        for line in fh:
            if line.startswith("chrom") or not line.strip():
                continue
            f = line.rstrip("\n").split("\t")
            chrom = f[0].strip()
            if chrom.startswith("b'") and chrom.endswith("'"):
                chrom = chrom[2:-1]          # b'1' -> 1
            d.setdefault(chrom, {})[int(f[1])] = float(f[4])   # start -> rate
    return d


male, fem = parse(male_f), parse(fem_f)
seqlens = {}
for line in open(seqlens_f):
    c, L = line.split()[:2]
    seqlens[c] = int(L)

rows = []
for c in sorted(set(male) | set(fem), key=lambda c: (len(c), c)):
    L = seqlens.get(c)
    if L is None:
        continue
    avg = {}
    for src in (male.get(c, {}), fem.get(c, {})):
        for s, r in src.items():
            avg.setdefault(s, []).append(r)
    avg = {s: sum(v) / len(v) for s, v in avg.items()}
    if not avg:
        continue
    mean_rate = sum(avg.values()) / len(avg)
    s = 0
    while s < L:
        e = min(s + WIN, L)
        rows.append((c, s, e, avg.get(s, mean_rate)))
        s = e

with open(out_f, "w") as o:
    for c, s, e, r in rows:
        o.write(f"{c}\t{s}\t{e}\t{r:.6e}\n")
print(f"wrote {len(rows)} windows across {len({r[0] for r in rows})} chroms -> {out_f}")
