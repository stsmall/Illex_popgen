#!/usr/bin/env python
"""Build one GONE2 bootstrap-replicate VCF by resampling CHROMOSOMES with replacement.

The 43 autosomal linkage groups present in gone_input.vcf are the natural block for a
recombination/LD-based method: each chromosome is an independent linkage group, so a
chromosome block bootstrap (draw 43 with replacement) is the principled resample.
Drawn chromosomes are relabelled to fresh sequential contig IDs (1..43) so that a
chromosome drawn twice becomes two independent linkage groups, and the header carries a
matching ##contig line (with the source length) for each. Physical positions are kept,
so GONE2's -r constant-rate physical->genetic map is unchanged per block.

Usage: make_boot_vcf.py <rep_index> <out.vcf>
"""
import sys
import random

BD = "/sietch_colab/data_share/illex/popgen_data/analysis/steps/10_gone2/bootstrap"
rep = int(sys.argv[1])
out = sys.argv[2]

# contig id -> length, from the saved header
lengths = {}
hdr_top = []      # header lines except ##contig, in order
for line in open(f"{BD}/header.txt"):
    line = line.rstrip("\n")
    if line.startswith("##contig="):
        # ##contig=<ID=1,length=104810120>
        inner = line[line.index("<") + 1:line.rindex(">")]
        d = dict(kv.split("=", 1) for kv in inner.split(","))
        lengths[d["ID"]] = d["length"]
    else:
        hdr_top.append(line)

# chromosomes that actually have data (body files present)
import os
present = sorted(int(f[:-5]) for f in os.listdir(f"{BD}/chroms") if f.endswith(".body"))
present = [str(c) for c in present]

rng = random.Random(1000 + rep)
draw = [rng.choice(present) for _ in present]     # 43 draws with replacement

with open(out, "w") as o:
    # header: all non-contig lines except the final #CHROM go first, then contigs, then #CHROM
    chrom_line = hdr_top[-1]
    for h in hdr_top[:-1]:
        o.write(h + "\n")
    for new_id, src in enumerate(draw, start=1):
        o.write(f"##contig=<ID={new_id},length={lengths[src]}>\n")
    o.write(chrom_line + "\n")
    for new_id, src in enumerate(draw, start=1):
        nid = str(new_id).encode()
        with open(f"{BD}/chroms/{src}.body", "rb") as bf:
            for bline in bf:
                tab = bline.index(b"\t")
                o.write(nid.decode())
                o.write(bline[tab:].decode())
print(f"rep {rep}: drew {len(draw)} chroms -> {out}", flush=True)
