"""structured_mask.py -- replicate diploSHIC fvTools structured genotype-masking
(getGenoMaskInfoInWins + maskGenos) in the slim_sim env (numpy + bcftools; no
scikit-allel), so the UNPHASED fvecVcf training path (ms_to_vcf_fvec.py) carries
the SAME real, correlated per-genotype missingness as the empirical scan --
instead of the old MCAR depth-model missingness.

Validated against diploSHIC's chr1 output: 1,408,502 SNPs -> 79 good windows,
calledFrac min 0.4543 / max 0.9686 (identical).

`build` mode reads a real cleaned VCF + accessibility FASTA and writes a packed
bank (.npz). `load_bank` + `apply_structured_mask` are used at fvec time.
"""
import argparse
import subprocess
import sys

import numpy as np

BCF = "/home/ssmall/bin/bcftools"


def read_accessibility_fasta(path, chrom):
    """bool array isAcc[chromLen]: True where the FASTA base != 'N' (accessible)."""
    buf, reading = [], False
    with open(path) as fh:
        for line in fh:
            if line[0] == ">":
                c = line[1:].strip().split()[0]
                if reading:
                    break
                reading = (c == str(chrom))
                continue
            if reading:
                buf.append(line.strip())
    if not buf:
        raise ValueError(f"chrom {chrom!r} not found in {path}")
    return (np.frombuffer("".join(buf).upper().encode(), np.uint8) != ord("N")).copy()


def _stream_missing(vcf, chrom, n_indiv):
    """Yield (pos, missRow bool[n_indiv]) for each site. missing = any allele '.'
    (matches allel/diploSHIC calledGenoFracAtSite: genotype missing if any allele<0)."""
    cmd = f"{BCF} query -r {chrom} -f '%POS[\\t%GT]\\n' {vcf}"
    p = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, text=True,
                         bufsize=1 << 20)
    for line in p.stdout:
        i = line.index("\t")
        gts = line[i + 1:].rstrip("\n").split("\t")
        row = np.fromiter(("." in g for g in gts), bool, n_indiv)
        yield int(line[:i]), row
    p.wait()


def build_bank(vcf, fasta, chrom, n_indiv=350, win_len=1_100_000,
               sub_win_len=100_000, cutoff=0.25, geno_cutoff=0.5):
    """Return (acc_list, miss_list) for GOOD windows, matching getGenoMaskInfoInWins:
      * a SNP is kept if its bp is FASTA-accessible AND calledFrac >= geno_cutoff;
        SNPs failing geno_cutoff flip their bp to inaccessible.
      * a window is kept if it has >=1 kept SNP AND every sub-window's accessible
        fraction >= cutoff.
    acc_list[i]  : bool[win_len]     accessibility for good window i
    miss_list[i] : bool[n_snp,n_indiv] real missing matrix for that window's kept SNPs
    """
    is_acc = read_accessibility_fasta(fasta, chrom)
    L = len(is_acc)
    pos, miss = [], []
    for p_, row in _stream_missing(vcf, chrom, n_indiv):
        pos.append(p_); miss.append(row)
    pos = np.array(pos, np.int64)
    miss = np.array(miss, bool)                       # (n_snp, n_indiv)
    called = 1.0 - miss.mean(axis=1)                  # calledFrac per SNP
    on_acc = is_acc[pos - 1]
    # flip low-call accessible SNP positions to inaccessible
    is_acc[pos[on_acc & (called < geno_cutoff)] - 1] = False
    kept = on_acc & (called >= geno_cutoff)
    kept_pos, kept_idx = pos[kept], np.nonzero(kept)[0]

    acc_list, miss_list = [], []
    last_win_end = L - (L % win_len)
    for w0 in range(0, last_win_end, win_len):
        sel = (kept_pos > w0) & (kept_pos <= w0 + win_len)
        if not np.any(sel):
            continue
        win = is_acc[w0:w0 + win_len]
        if not all(win[s:s + sub_win_len].sum() / sub_win_len >= cutoff
                   for s in range(0, win_len, sub_win_len)):
            continue
        acc_list.append(win.copy())
        miss_list.append(miss[kept_idx[sel]])
    sys.stderr.write("[structured_mask] %s: %d SNPs -> %d good windows "
                     "(avg %.0f kept SNPs/win)\n" %
                     (chrom, len(pos), len(acc_list),
                      np.mean([m.shape[0] for m in miss_list]) if miss_list else 0))
    return acc_list, miss_list


def save_bank(path, acc_list, miss_list, win_len):
    """Pack to a compact .npz (bit-packed accessibility + missing)."""
    acc_packed = np.stack([np.packbits(a) for a in acc_list])
    miss_packed = np.concatenate([np.packbits(m.ravel()) for m in miss_list])
    offsets = np.cumsum([0] + [np.packbits(m.ravel()).size for m in miss_list])
    nsnps = np.array([m.shape[0] for m in miss_list], np.int64)
    np.savez(path, acc_packed=acc_packed, miss_packed=miss_packed,
             offsets=offsets, nsnps=nsnps, win_len=win_len,
             n_indiv=miss_list[0].shape[1] if miss_list else 0)


def load_bank(path):
    """Return (acc_list bool[win_len], miss_list bool[nsnp,n_indiv]) from a saved bank."""
    z = np.load(path)
    win_len = int(z["win_len"]); n_indiv = int(z["n_indiv"])
    acc_list = [np.unpackbits(a)[:win_len].astype(bool) for a in z["acc_packed"]]
    off, nsnps, mp = z["offsets"], z["nsnps"], z["miss_packed"]
    miss_list = []
    for i, ns in enumerate(nsnps):
        chunk = mp[off[i]:off[i + 1]]
        miss_list.append(np.unpackbits(chunk)[:ns * n_indiv].reshape(ns, n_indiv).astype(bool))
    return acc_list, miss_list


def apply_structured_mask(g, pos01, L, acc_slice, missmat):
    """Mask one sim rep with a drawn real window (diploSHIC maskGenos semantics).

    g        : int8 (n_indiv, n_snp) dosage 0/1/2 from parse_ms
    pos01    : float (n_snp,) positions in [0,1)
    acc_slice: bool[L] accessibility of the drawn real window
    missmat  : bool[n_real, n_indiv] real missing pattern of the drawn window
    Returns (g_masked int8 (n_indiv, n_kept), pos01_kept). Accessible sim SNP i takes
    the missing pattern of real SNP (i mod n_real); wrapping matches maskGenos.
    """
    bp = np.clip((pos01 * L).astype(np.int64), 0, L - 1)
    keep = acc_slice[bp]
    g, pos01 = g[:, keep], pos01[keep]
    m = g.shape[1]
    if m == 0:
        return g, pos01
    miss_rows = missmat[np.arange(m) % missmat.shape[0]]      # (m, n_indiv)
    g = np.where(miss_rows.T, np.int8(-1), g)                 # set missing genotypes
    return g, pos01


def main():
    ap = argparse.ArgumentParser(description="Build a structured geno-mask bank (.npz).")
    ap.add_argument("--vcf", required=True); ap.add_argument("--fasta", required=True)
    ap.add_argument("--chrom", required=True); ap.add_argument("--out", required=True)
    ap.add_argument("--n-indiv", type=int, default=350, dest="n_indiv")
    ap.add_argument("--win-len", type=int, default=1_100_000, dest="win_len")
    ap.add_argument("--sub-win-len", type=int, default=100_000, dest="sub_win_len")
    ap.add_argument("--cutoff", type=float, default=0.25)
    ap.add_argument("--geno-cutoff", type=float, default=0.5, dest="geno_cutoff")
    a = ap.parse_args()
    acc, miss = build_bank(a.vcf, a.fasta, a.chrom, a.n_indiv, a.win_len,
                           a.sub_win_len, a.cutoff, a.geno_cutoff)
    save_bank(a.out, acc, miss, a.win_len)
    sys.stderr.write("[structured_mask] wrote bank -> %s\n" % a.out)


if __name__ == "__main__":
    main()
