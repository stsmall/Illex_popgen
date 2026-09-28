#!/usr/bin/env python3
"""ms_to_vcf_fvec.py -- regenerate diploSHIC TRAINING feature vectors via the SCAN path
(fvecVcf on UNPHASED VCFs) so training matches the empirical scan. Root cause of the broken
scan: fvecSim (phased ms) vs fvecVcf (unphased VCF) compute diplo_ZnS etc. ~40x differently
(diploshic-zns-train-scan-mismatch); fix = train on the same you scan on (user 2026-08-08).

Reads one or more ms.gz files (same class-position), converts each replicate to an unphased
diploid dosage with step-13 missing, lays CHUNK reps contiguously on a synthetic 'sim' chrom
(rep j at [j*L,(j+1)*L)), writes ONE VCF + tiled real chr43 mask, runs ONE `fvecVcf diploid`
(amortizing its ~10s startup), and keeps the aligned window per rep (classifiedWinStart=j*L+500001;
100kb-step cross-boundary windows discarded). Fast byte-lookup VCF writer.

Usage: ms_to_vcf_fvec.py --out <fvec> [--n-per-file N] [--chunk 50] [--seed S] [--tmp DIR] --ms f1 [f2 ...]
"""
import argparse, gzip, os, subprocess, sys, tempfile
from itertools import islice
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import ascertain as asc_mod          # noqa: E402
import structured_mask as smask      # noqa: E402
DE = "/home/ssmall/miniforge3/envs/diploshic_env/bin/diploSHIC"
BGZIP = "/home/ssmall/bin/bgzip"; TABIX = "/home/ssmall/bin/tabix"; SAM = "/home/ssmall/bin/samtools"
D = "/sietch_colab/data_share/illex/popgen_data/analysis/steps/14_sweep_seqmodel"
S2P = f"{D}/results/task17_maskval/s2p.tsv"
REALMASK = f"{D}/results/empirical_scan/masks/mask.43.fa"
L = 1_100_000; NSUB = 11; MASK_CHR, MASK_START = "43", 60_000_000
_GTLUT = np.frombuffer(b"0/0\t0/1\t1/1\t./.\t", dtype=np.uint8).reshape(4, 4)


def parse_ms(path, nmax, n_dip=None):
    # binary block-parse: b"".join(rows) -> ONE frombuffer -> reshape. ~10-50x faster than
    # text-mode readline + per-row .encode()/np.stack, which choked at full-SFS density
    # (~15GB decompressed text per neutral file -> hours/file). Same output (diploid dosage, pos).
    with gzip.open(path, "rb") as fh:
        cnt = 0
        for line in fh:
            if not line.startswith(b"//"):
                continue
            if cnt >= nmax:
                break
            k = int(fh.readline().split()[1])
            if k == 0:
                fh.readline(); continue
            pos = np.array(fh.readline().split()[1:], dtype=float)
            rows = []
            while True:
                r = fh.readline()
                if not r or r == b"\n" or r.strip() == b"":
                    break
                rows.append(r.rstrip(b"\r\n"))
            m = len(rows)
            hap = (np.frombuffer(b"".join(rows), np.uint8).reshape(m, k) == ord("1")).astype(np.int8)
            ndip = m // 2 if n_dip is None else min(n_dip, m // 2)
            cnt += 1
            yield hap[:ndip * 2].reshape(ndip, 2, k).sum(1).astype(np.int8), pos


def build_mask(nbatch, out_fa):
    reg = f"{MASK_CHR}:{MASK_START + 1}-{MASK_START + L}"
    r = subprocess.run([SAM, "faidx", REALMASK, reg], capture_output=True, text=True, check=True)
    body = "".join(r.stdout.splitlines()[1:])
    if len(body) < L:
        body += "N" * (L - len(body))
    body *= nbatch
    with open(out_fa, "w") as fh:
        fh.write(">sim\n")
        for i in range(0, len(body), 60):
            fh.write(body[i:i + 60] + "\n")
    subprocess.run([SAM, "faidx", out_fa], check=True)


def build_structured_mask(acc_slices, out_fa):
    """Write per-rep accessibility (drawn real windows) as one FASTA (N=inaccessible)."""
    body = np.concatenate(acc_slices)
    b = np.where(body, np.uint8(ord("A")), np.uint8(ord("N"))).tobytes()
    with open(out_fa, "w") as fh:
        fh.write(">sim\n")
        for i in range(0, len(b), 60):
            fh.write(b[i:i + 60].decode("ascii") + "\n")
    subprocess.run([SAM, "faidx", out_fa], check=True)


def write_rep_rows(fh, g, pos, base, rng, model, mask=None):
    """Append one rep's VCF rows (rep occupies [base, base+L)); bake in missing data.
    mask=(acc_slice, missmat) -> REAL structured missing (diploSHIC maskGenos);
    mask=None -> legacy MCAR depth-model missing."""
    if mask is not None:
        g, pos = smask.apply_structured_mask(g, pos, L, mask[0], mask[1])
    else:
        depth = model["D"][rng.choice(len(model["D"]), size=g.shape, p=model["P"])]
        g = np.where(depth < 2, np.int8(-1), g)
    ac = np.where(g < 0, 0, g).sum(0); an = 2 * (g >= 0).sum(0)
    keep = (ac > 0) & (ac < an)
    mm = model.get("min_maf")
    if mm:                                     # optional MAF cut; default None -> full SFS (unchanged behavior)
        with np.errstate(invalid="ignore", divide="ignore"):
            maf = np.minimum(ac, an - ac) / np.maximum(an, 1)
        keep &= maf > mm
    g, pos = g[:, keep], pos[keep]
    ipos = base + (pos * L).astype(np.int64) + 1
    order = np.argsort(ipos, kind="stable"); ipos = ipos[order]; g = g[:, order]
    uniq = np.concatenate(([True], np.diff(ipos) > 0)); ipos = ipos[uniq]; g = g[:, uniq]
    nsnp = g.shape[1]
    gt = _GTLUT[np.where(g < 0, 3, g).astype(np.intp).T].copy()   # (nsnp, ndip, 4)
    gt[:, -1, 3] = ord("\n")
    gtflat = gt.reshape(nsnp, -1)
    buf = []
    for i in range(nsnp):
        buf.append(f"sim\t{ipos[i]}\t.\tA\tT\t.\tPASS\t.\tGT\t".encode("ascii"))
        buf.append(gtflat[i].tobytes())
        if len(buf) >= 40000:
            fh.write(b"".join(buf)); buf = []
    if buf:
        fh.write(b"".join(buf))


def process_chunk(reps, samples, rng, model, td, bank=None):
    B = len(reps)
    vcf = f"{td}/c.vcf"; maskfa = f"{td}/m.fa"; fvec = f"{td}/c.fvec"
    if bank is not None:
        acc_list, miss_list = bank
        draws = rng.integers(len(acc_list), size=B)       # draw a real window per rep (with replacement)
        build_structured_mask([acc_list[d] for d in draws], maskfa)
    else:
        build_mask(B, maskfa)
    with open(vcf, "wb") as fh:
        fh.write(("##fileformat=VCFv4.2\n##contig=<ID=sim>\n"
                  '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
                  "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + "\t".join(samples) + "\n").encode())
        for j, (g, pos) in enumerate(reps):
            mask = (acc_list[draws[j]], miss_list[draws[j]]) if bank is not None else None
            write_rep_rows(fh, g, pos, j * L, rng, model, mask)
    subprocess.run([BGZIP, "-f", vcf], check=True)
    subprocess.run([TABIX, "-f", "-p", "vcf", vcf + ".gz"], check=True)
    subprocess.run([DE, "fvecVcf", "diploid", vcf + ".gz", "sim", str(B * L), fvec,
                    "--targetPop", "illex", "--sampleToPopFileName", S2P,
                    "--winSize", str(L), "--numSubWins", str(NSUB), "--maskFileName", maskfa,
                    "--unmaskedFracCutoff", "0.25", "--unmaskedGenoFracCutoff", "0.5"],
                   check=True, capture_output=True)
    want = {j * L + 500001: j for j in range(B)}
    header, got = None, {}
    with open(fvec) as fh:
        h = fh.readline().rstrip("\n").split("\t"); header = h
        si = h.index("classifiedWinStart")
        for line in fh:
            f = line.rstrip("\n").split("\t")
            cs = int(f[si])
            if cs in want:
                got[want[cs]] = f
    return header, [got[j] for j in range(B) if j in got]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True); ap.add_argument("--ms", nargs="+", required=True)
    ap.add_argument("--n-per-file", type=int, default=10**9); ap.add_argument("--chunk", type=int, default=50)
    ap.add_argument("--seed", type=int, default=1); ap.add_argument("--tmp", default=None)
    ap.add_argument("--n-dip", type=int, default=None, dest="n_dip",
                    help="subsample each sim rep to this many diploids (first 2*n_dip haplotypes; "
                         "coalescent haps are exchangeable). Default None = all m//2 (=350). Also truncates "
                         "the VCF sample list so column count matches (fvecVcf counts only illex samples present).")
    ap.add_argument("--min-maf", type=float, default=None, dest="min_maf",
                    help="optional MAF cut on emitted sites (default None = full SFS; this script does "
                         "NOT filter by frequency otherwise -- training MAF comes from generate.py's ascertain).")
    ap.add_argument("--geno-mask-bank", default=None, dest="geno_mask_bank",
                    help="path to a structured geno-mask bank (.npz from structured_mask.py). When given, "
                         "bakes REAL correlated per-genotype missingness (diploSHIC maskGenos) into the sim "
                         "VCF instead of the MCAR depth model -> matches the empirical scan's missingness.")
    a = ap.parse_args()
    rng = np.random.default_rng(a.seed); model = asc_mod.load_model()
    model["min_maf"] = a.min_maf              # None -> no MAF filter (overrides the config default, unused here before)
    bank = smask.load_bank(a.geno_mask_bank) if a.geno_mask_bank else None
    samples = [s.split("\t")[0] for s in open(S2P) if s.strip()]
    if a.n_dip is not None:
        samples = samples[:a.n_dip]      # match VCF column count to the subsampled ndip

    def rep_stream():                       # lazy: never materialize all reps at once
        for f in a.ms:
            yield from parse_ms(f, a.n_per_file, a.n_dip)
    sys.stderr.write(f"[ms2vcf] streaming reps from {len(a.ms)} files, chunk={a.chunk}\n")
    rows, header = [], None
    it = rep_stream(); ci = 0
    with tempfile.TemporaryDirectory(dir=a.tmp) as td:
        while True:
            batch = list(islice(it, a.chunk))     # hold only 'chunk' reps in RAM (was: all reps -> 7.6GB
            if not batch:                          # for neutral, thrashing NUMA under concurrency)
                break
            h, r = process_chunk(batch, samples, rng, model, td, bank)
            header = header or h; rows.extend(r)
            del batch
            sys.stderr.write(f"[ms2vcf]   chunk {ci}: {len(rows)} rows total\n"); ci += 1
    with open(a.out, "w") as o:
        o.write("\t".join(header) + "\n")
        for r in rows:
            o.write("\t".join(r) + "\n")
    sys.stderr.write(f"[ms2vcf] wrote {len(rows)} rows -> {a.out}\n")


if __name__ == "__main__":
    main()
