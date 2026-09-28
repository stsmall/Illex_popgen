"""Window extraction + gzipped-ms export for the SLiM 9-class sim generator.

Turns a reconstructed tree sequence (n=700 haplotypes over a 2.2 Mb locus with the
sweep at 1.1 Mb) into the 11 diploSHIC/partialSHIC training windows and writes them
as gzipped ms -- the input format both fvec toolchains consume.

Conventions (see PLAN_03 Global Constraints):
  winSize = 1.1 Mb, 11 subwindows of 100 kb, sweep at 1.1 Mb (centre of the 2.2 Mb
  locus). To place the sweep at the centre of subwindow i, slice [offset, offset+win)
  with offset = sweep_pos - (i + 0.5) * sub.
Assumes biallelic (infinite-sites) mutations: genotypes are 0/1.
"""
from __future__ import annotations

import gzip

import numpy as np


def window_offsets(win: int = 1_100_000, sub: int = 100_000,
                   n_sub: int = 11, sweep_pos: int = 1_100_000):
    """Return ``[(subwin_index, offset), ...]`` for the ``n_sub`` analysis windows.

    ``offset = sweep_pos - (i + 0.5) * sub`` places the sweep at the centre of
    subwindow ``i`` of the window ``[offset, offset + win)``.
    """
    return [(i, int(round(sweep_pos - (i + 0.5) * sub))) for i in range(n_sub)]


def genotype_and_positions(ts):
    """``(G [n_site, n_hap] int8, positions [n_site] float)`` for a tree sequence.

    Computed once so a caller (the generator) can window 11 times without
    recomputing the genotype matrix.
    """
    G = ts.genotype_matrix().astype(np.int8)        # (n_site, n_hap)
    pos = np.asarray(ts.tables.sites.position, dtype=float)
    return G, pos


def extract_window(ts, offset: int, win: int):
    """Genotypes + renormalised positions for SNPs in ``[offset, offset + win)``.

    Returns ``(geno [n_hap, n_snp] int8, positions_rel [n_snp] float in (0,1))``,
    with ``positions_rel = (pos - offset) / win`` and sites sorted by position.
    """
    G, pos = genotype_and_positions(ts)
    return window_from(G, pos, offset, win)


def window_from(G, pos, offset: int, win: int):
    """Window a precomputed ``(G [n_site, n_hap], pos [n_site])`` pair.

    Same return contract as :func:`extract_window`; use this in the generator to
    avoid recomputing the genotype matrix for each of the 11 windows.
    """
    pos = np.asarray(pos, dtype=float)
    in_win = (pos >= offset) & (pos < offset + win)
    idx = np.nonzero(in_win)[0]
    n_hap = G.shape[1]
    if idx.size == 0:
        return np.zeros((n_hap, 0), dtype=np.int8), np.zeros(0, dtype=float)
    order = idx[np.argsort(pos[idx])]
    geno = np.asarray(G)[order, :].T.astype(np.int8)   # (n_hap, n_snp)
    positions_rel = (pos[order] - offset) / float(win)
    return geno, positions_rel


def write_ms_gz(records, path: str, n_hap: int, seed: int = 0,
                append: bool = False):
    """Write ms-format replicates to a gzip file.

    ``records``: iterable of ``(geno [n_hap, n_snp] int8, positions_rel [n_snp])``.
    Genotype encoding: ``0`` -> ancestral, ``>0`` -> derived, ``<0`` (e.g. -1) ->
    missing, emitted as ``?`` (post-ascertainment windows carry missing data).

    ``append=False`` (default) writes the ``ms n_hap n_reps`` header then the
    blocks (mode ``wt``). ``append=True`` skips the header and appends blocks
    (mode ``at``) so a campaign can stream many sims into one per-class file;
    diploSHIC/partialSHIC parse replicate blocks by the ``//`` delimiter, so the
    header ``n_reps`` need not be exact when appending.
    """
    records = list(records)
    mode = "at" if append else "wt"
    # transient ms (deleted after fvec); level 1 is ~10x faster than the default 9 on this
    # highly-compressible 0/1 text and the size difference is immaterial.
    with gzip.open(path, mode, compresslevel=1) as fh:
        if not append:
            fh.write(f"ms {n_hap} {len(records)}\n{seed}\n\n")
        for geno, positions_rel in records:
            geno = np.asarray(geno, dtype=np.int8)
            if geno.shape[0] != n_hap:
                raise ValueError(
                    f"record has {geno.shape[0]} haplotypes, expected {n_hap}")
            k = geno.shape[1]
            fh.write("//\n")
            fh.write(f"segsites: {k}\n")
            if k == 0:
                fh.write("\n")            # ms prints a blank line after "segsites: 0"
                continue
            order = np.argsort(positions_rel)
            pos_sorted = np.asarray(positions_rel, dtype=float)[order]
            geno_sorted = geno[:, order]
            fh.write("positions: " + " ".join(np.char.mod("%.6f", pos_sorted)) + "\n")
            # vectorized 0/1/? encoding (pure-Python per-char loops are O(hap*snp) and
            # infeasible at illex density): map <0->'?', 0->'0', >0->'1' via a uint8 buffer.
            codes = np.full(geno_sorted.shape, ord("0"), dtype=np.uint8)
            codes[geno_sorted > 0] = ord("1")
            codes[geno_sorted < 0] = ord("?")
            fh.write(b"\n".join(codes[i].tobytes()
                                for i in range(codes.shape[0])).decode("ascii"))
            fh.write("\n\n")
