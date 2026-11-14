#!/usr/bin/env python
"""Collinear joint-spectrum test of the 40 + 40 sampling composition (review M4), with jsfs_balance.sh.

draw   three division-balanced 40 + 40 sets: 40 individuals allocated over divisions in proportion to all 619
       (largest remainder, capped by the smaller of the 2a/2a and 2+/2+ counts in each division), the same counts drawn
       at random from each karyotype group (seeds 20261114-20261116); the original 04_angsd_chr2 AA40/BB40 lists are
       copied as *_orig. Writes {AA,BB}40_{orig,bal1,bal2,bal3}.txt.
cells  the 12 folded joint-spectrum cells of each AA_BB.<set>.2dSFS with msinv's jsfs_from_2d kernel, copied below (projection
       to 20 chromosomes per class and binning as for the inversion statistics), and the private-rare pair
       jsfs_01 / jsfs_10.
Usage: jsfs_balance.py {draw,cells}   (in analysis/steps/03_karyotype/m4_qc/jsfs_balance)
"""
import shutil
import sys

import numpy as np

ANA = "/sietch_colab/data_share/illex/popgen_data/analysis/steps"
SETS = ("orig", "bal1", "bal2", "bal3")

# Binning kernels copied verbatim from msinv illex/abc/statvec.py (_hyper_binned, _jsfs_fold, JSFS_CELLS) and
# illex/abc/observed.py (jsfs_from_2d) at commit 36a8b77; importing illex needs msinv's compiled core.
PROJ = 20
JSFS_BINS = ((0, 0), (1, 4), (5, 15), (16, 19), (20, 20))


def _jsfs_cells():
    nb = len(JSFS_BINS)
    cells = []
    for a in range(nb):
        for b in range(nb):
            mirror = (nb - 1 - a, nb - 1 - b)
            if (a, b) <= mirror and (a, b) != (0, 0):
                cells.append((a, b))
    return tuple(cells)


JSFS_CELLS = _jsfs_cells()


def _hyper_binned(counts, n):
    from scipy.stats import hypergeom
    kk = np.arange(PROJ + 1)
    u, inv = np.unique(counts, return_inverse=True)
    pmf = hypergeom.pmf(kk[None, :], n, u[:, None], PROJ)
    binned = np.stack([pmf[:, lo:hi + 1].sum(1) for lo, hi in JSFS_BINS], axis=1)
    return binned[inv]


def _jsfs_fold(j):
    nb = len(JSFS_BINS)
    fol = {}
    for a, b in JSFS_CELLS:
        ma, mb = nb - 1 - a, nb - 1 - b
        fol[(a, b)] = j[a, b] + (j[ma, mb] if (ma, mb) != (a, b) else 0.0)
    tot = sum(fol.values())
    return {f"jsfs_{a}{b}": float(v / tot) for (a, b), v in fol.items()}


def jsfs_from_2d(m):
    m = np.asarray(m, float)
    bi = _hyper_binned(np.arange(m.shape[0]), m.shape[0] - 1)
    bs = _hyper_binned(np.arange(m.shape[1]), m.shape[1] - 1)
    return _jsfs_fold(bi.T @ m @ bs)



def draw():
    import pandas as pd
    K = pd.read_csv("/sietch_colab/data_share/illex/popgen_data/illex_manuscript/data/qc/karyotypes_619.tsv", sep="\t",
                    comment="#").set_index("sample").karyotype_final
    div = pd.read_csv("/sietch_colab/data_share/illex/popgen_data/pggpu_illex/popstats/pops.tsv", sep="\t",
                      header=None, index_col=0)[1].reindex(K.index)
    bam = {l.strip().split("/")[-1].replace(".bam", ""): l.strip()
           for l in open(f"{ANA}/04_angsd_chr2/bamlists/bakerMeta_bamlist.txt") if l.strip()}
    aa, bb = K[K == "AA"].index, K[K == "BB"].index
    avail = pd.DataFrame({"AA": div[aa].value_counts(), "BB": div[bb].value_counts()}).fillna(0).astype(int)
    cap = avail.min(axis=1)
    w = div.value_counts().reindex(avail.index).fillna(0)
    alloc = np.minimum(np.floor(40 * w / w.sum()).astype(int), cap)
    while alloc.sum() < 40:
        r = (40 * w / w.sum() - alloc).where((cap - alloc) > 0, -1)
        alloc[r.idxmax()] += 1
    print(alloc.to_dict())
    for rep in (1, 2, 3):
        rng = np.random.default_rng(20261113 + rep)
        for g, idx in (("AA", aa), ("BB", bb)):
            pick = []
            for d, n in alloc.items():
                if n:
                    pick += list(rng.choice(sorted(div[idx][div[idx] == d].index), n, replace=False))
            open(f"{g}40_bal{rep}.txt", "w").write("\n".join(bam[s] for s in pick) + "\n")
    for g in ("AA", "BB"):
        shutil.copy(f"{ANA}/04_angsd_chr2/bamlists/{g}40_bamlist.txt", f"{g}40_orig.txt")


def cells():
    print("set   jsfs_01  jsfs_10  01/10   (original 45-Mb control: 0.352, 0.272)")
    for s in SETS:
        v = np.array(open(f"AA_BB.{s}.2dSFS").read().split(), float)
        assert v.size == 81 * 81
        c = jsfs_from_2d(v.reshape(81, 81))
        print(f"{s:5s} {c['jsfs_01']:.3f}    {c['jsfs_10']:.3f}    {c['jsfs_01'] / c['jsfs_10']:.2f}")


if __name__ == "__main__":
    {"draw": draw, "cells": cells}[sys.argv[1]]()
