#!/usr/bin/env python
"""Aggregate GONE2 bootstrap replicate Ne trajectories into a per-generation CI band.

Reads ne/rep*_GONE2_Ne, aligns on Generation, and writes gone2_band.tsv with the
per-generation median and 2.5/97.5 percentiles across replicates. Generations 1-4 are
kept in the file but flagged (the figure drops them as the known GONE2 recent-edge
artifact, matching the point-estimate handling).
"""
import glob
import os
import numpy as np

BD = "/sietch_colab/data_share/illex/popgen_data/analysis/steps/10_gone2/bootstrap"
files = sorted(glob.glob(f"{BD}/ne/rep*_GONE2_Ne"))
assert files, "no replicate Ne files found"

gens = None
mat = []
for f in files:
    g, ne = [], []
    with open(f) as fh:
        next(fh)  # header: Generation\tNe_diploids
        for line in fh:
            a, b = line.split()
            g.append(int(a)); ne.append(float(b))
    g = np.array(g); ne = np.array(ne)
    if gens is None:
        gens = g
    if len(g) != len(gens) or not np.array_equal(g, gens):
        # align on common generations
        common = np.intersect1d(gens, g)
        idx = np.isin(g, common)
        ne = ne[idx]; g = g[idx]
        gens = common
    mat.append(ne)

# re-align all to final gens length
rows = []
for f in files:
    g, ne = [], []
    with open(f) as fh:
        next(fh)
        for line in fh:
            a, b = line.split(); g.append(int(a)); ne.append(float(b))
    d = dict(zip(g, ne))
    rows.append([d[x] for x in gens])
mat = np.array(rows)                      # (nrep, ngen)

lo = np.percentile(mat, 2.5, axis=0)
med = np.percentile(mat, 50, axis=0)
hi = np.percentile(mat, 97.5, axis=0)
point = None
pf = "/sietch_colab/data_share/illex/popgen_data/analysis/steps/10_gone2/gone_input.vcf_GONE2_Ne"
if os.path.exists(pf):
    pd_ = {}
    with open(pf) as fh:
        next(fh)
        for line in fh:
            a, b = line.split(); pd_[int(a)] = float(b)
    point = np.array([pd_.get(x, np.nan) for x in gens])

out = f"{BD}/gone2_band.tsv"
with open(out, "w") as o:
    o.write("Generation\tNe_point\tNe_median\tNe_2.5\tNe_97.5\tn_rep\n")
    for i, x in enumerate(gens):
        p = "" if point is None else f"{point[i]:.6g}"
        o.write(f"{x}\t{p}\t{med[i]:.6g}\t{lo[i]:.6g}\t{hi[i]:.6g}\t{mat.shape[0]}\n")
print(f"wrote {out}  ({mat.shape[0]} replicates, {len(gens)} generations)")
# report band width at a few generations (>=5)
for x in [5, 10, 25, 50, 100]:
    if x in gens:
        i = list(gens).index(x)
        print(f"  gen {x:3d}: median={med[i]:,.0f}  95% CI=[{lo[i]:,.0f}, {hi[i]:,.0f}]  "
              f"width factor={hi[i]/lo[i]:.2f}x")
