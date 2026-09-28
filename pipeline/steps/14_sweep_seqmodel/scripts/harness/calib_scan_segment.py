#!/usr/bin/env python3
"""calib_scan_segment.py -- PLAN_02 Phase B: tile a calibration segment into the SAME
diploSHIC windows the EMPIRICAL scan uses, and emit per-window class probabilities.

Takes a recapitated+overlaid segment tree sequence (from calib_run_sim.py) and runs it
through the *identical* empirical path: ascertain (step-13 het->hom miscall + realistic
missingness) -> VCF -> `diploSHIC fvecVcf` @ 1.1Mb/100kb -> `diploSHIC predict`. Because
the empirical scan (scan_empirical.sh) also uses fvecVcf with the same flags and the real
chr43 accessibility mask (sliced here to the segment), the resulting per-window prob
vectors live in the identical feature space as the real scan -- so HMM emissions calibrated
on these transfer directly. XPOS (the true sweep centre) is recorded so Phase C can label
each window N/LL/C/LR by distance from it.

Optional `--path sim` reproduces the TRAINING path (tiled ms -> fvecSim + genoMask) on the
SAME ascertained matrix, for a golden fvecSim-vs-fvecVcf parity spot check.

Output: <outdir>/<tag>.preds  (diploSHIC predict table, one row per 100kb window).
Deterministic from --seed.
"""
import argparse, gzip, os, subprocess, sys
import numpy as np
import tskit

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import ascertain as asc_mod          # noqa: E402  (load_model: depth_hist, drop_arr, min_maf)
import ms_export                     # noqa: E402  (window_from, write_ms_gz)

DE = "/home/ssmall/miniforge3/envs/diploshic_env/bin/diploSHIC"
BGZIP = "/home/ssmall/bin/bgzip"
TABIX = "/home/ssmall/bin/tabix"
SAMTOOLS = "/home/ssmall/bin/samtools"
D = "/sietch_colab/data_share/illex/popgen_data/analysis/steps/14_sweep_seqmodel"
S2P_DEF = f"{D}/results/task17_maskval/s2p.tsv"
MODEL_DEF = f"{D}/results/fullrun/illexModel"
REALMASK_DEF = f"{D}/results/empirical_scan/masks/mask.43.fa"      # real chr43 accessibility
# training genoMask (chr1 real missing pattern) + mask FASTA, used only by --path sim
TRAIN_VCFMASK = f"{D}/results/task17_maskval/cleaned.1.vcf.gz"
TRAIN_MASKFA = f"{HERE}/talapas/mask/mask.fa"


def segment_ascertain(ts, seed, bake_missing=True, maf_filter=False):
    """Full step-13 ascertainment over a WHOLE segment (one depth draw per genotype, so a
    SNP has a single depth across all tiling windows). Returns diploid dosage g[ndip,nsnp]
    in {0,1,2} with missing = -1, and absolute positions (bp within the segment).
    het->hom miscall on called (depth>=2) genos; depth<2 -> missing iff bake_missing.
    maf_filter: apply the training MAF>min_maf cut over COMPLETE genotypes (matches ascertain.py
    exactly -- the model was trained on MAF-filtered sims, so the empirical scan must match too).
    Keeps polymorphic sites (ac>0 & ac<an over called)."""
    rng = np.random.default_rng(seed)
    m = asc_mod.load_model()
    D_, P_, drop_arr, min_maf = m["D"], m["P"], m["drop_arr"], m["min_maf"]
    # sim_mutations' default (JC69) yields allele INDICES 0..k (recurrent muts); diploSHIC is
    # biallelic ancestral/derived, so fold any derived index -> 1 (the ms_export convention).
    G = (ts.genotype_matrix() > 0).astype(np.int64)           # (nsite, nhap) 0/1
    pos = np.asarray(ts.tables.sites.position, dtype=float)   # bp within [0,L)
    nhap = G.shape[1]; ndip = nhap // 2
    g = G[:, :ndip * 2].reshape(G.shape[0], ndip, 2).sum(2).T   # (ndip, nsite) dosage 0/1/2
    depth = D_[rng.choice(len(D_), size=g.shape, p=P_)]
    called = depth >= 2
    drop = drop_arr[np.clip(depth.astype(int), 0, len(drop_arr) - 1)]
    het = g == 1
    het2hom = het & called & (rng.random(g.shape) < drop)
    g = np.where(het2hom, 2 * rng.integers(0, 2, g.shape), g)   # complete 0/1/2
    if maf_filter:                                            # MAF over COMPLETE (as ascertain.py)
        ac0 = g.sum(0); an0 = 2 * ndip
        kmaf = np.minimum(ac0, an0 - ac0) / an0 > min_maf
        g, pos, called = g[:, kmaf], pos[kmaf], called[:, kmaf]
    if bake_missing:
        g = np.where(~called, -1, g)                          # depth<2 -> missing
    ac = np.where(g < 0, 0, g).sum(0)
    an = 2 * (g >= 0).sum(0)
    keep = (ac > 0) & (ac < an)                               # polymorphic among called
    g, pos = g[:, keep], pos[keep]
    ipos = (pos + 1).astype(int)                              # 1-based VCF POS
    _, first = np.unique(ipos, return_index=True)             # dedupe integer POS (tabix-safe)
    first.sort()
    return g[:, first].astype(np.int8), pos[first]


def write_vcf(g, pos, start, chrom_name, samples, out_vcf):
    """Write a bgzipped+tabixed VCF from diploid dosage g[ndip,nsnp] ({0,1,2,-1}) at
    POS = pos+1 (1-based within the synthetic segment). GT: 0->0/0 1->0/1 2->1/1 -1->./."""
    ndip, nsnp = g.shape
    assert len(samples) == ndip, f"{len(samples)} samples != {ndip} diploids"
    GT = np.empty(g.shape, dtype=object)
    GT[g == 0] = "0/0"; GT[g == 1] = "0/1"; GT[g == 2] = "1/1"; GT[g < 0] = "./."
    plain = out_vcf[:-3] if out_vcf.endswith(".gz") else out_vcf
    ipos = (pos + 1).astype(int)
    with open(plain, "w") as fh:
        fh.write("##fileformat=VCFv4.2\n")
        fh.write(f"##contig=<ID={chrom_name}>\n")
        fh.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
        fh.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t"
                 + "\t".join(samples) + "\n")
        for j in range(nsnp):
            fh.write(f"{chrom_name}\t{ipos[j]}\t.\tA\tT\t.\tPASS\t.\tGT\t"
                     + "\t".join(GT[:, j]) + "\n")
    subprocess.run([BGZIP, "-f", plain], check=True)
    subprocess.run([TABIX, "-f", "-p", "vcf", out_vcf], check=True)


def slice_real_mask(realmask_fa, chrom, start, L, chrom_name, out_fa):
    """Slice the real per-chrom accessibility mask to [start, start+L) and rename to
    chrom_name, so fvecVcf masks the sim windows with the region's TRUE accessibility."""
    reg = f"{chrom}:{start + 1}-{start + L}"
    r = subprocess.run([SAMTOOLS, "faidx", realmask_fa, reg],
                       capture_output=True, text=True, check=True)
    body = "".join(r.stdout.splitlines()[1:])
    if len(body) < L:                                         # pad past chrom end with N
        body = body + "N" * (L - len(body))
    with open(out_fa, "w") as fh:
        fh.write(f">{chrom_name}\n")
        for i in range(0, len(body), 60):
            fh.write(body[i:i + 60] + "\n")
    subprocess.run([SAMTOOLS, "faidx", out_fa], check=True)


def run_vcf_path(g, pos, a, samples, chrom_name):
    vcf = os.path.join(a.outdir, a.tag + ".seg.vcf.gz")
    maskfa = os.path.join(a.outdir, a.tag + ".segmask.fa")
    fvec = os.path.join(a.outdir, a.tag + ".vcf.fvec")
    preds = os.path.join(a.outdir, a.tag + ".vcf.preds")
    write_vcf(g, pos, a.start, chrom_name, samples, vcf)
    slice_real_mask(a.realmask, a.chrom, a.start, a.L, chrom_name, maskfa)
    subprocess.run([DE, "fvecVcf", "diploid", vcf, chrom_name, str(a.L), fvec,
                    "--targetPop", "illex", "--sampleToPopFileName", a.s2p,
                    "--winSize", str(a.win), "--numSubWins", str(a.numSubWins),
                    "--maskFileName", maskfa,
                    "--unmaskedFracCutoff", "0.25", "--unmaskedGenoFracCutoff", "0.5"],
                   check=True)
    subprocess.run([DE, "predict", a.model + ".json", a.model + ".weights.h5",
                    fvec, preds, "--numSubWins", str(a.numSubWins)], check=True)
    return preds


def run_sim_path(g, pos, a):
    """Training path on the SAME ascertained matrix: tiled ms -> fvecSim + genoMask."""
    # g[ndip,nsnp] -> hap matrix [nsite,nhap] for ms_export.window_from
    ndip, nsnp = g.shape
    a_hap = np.where(g < 0, -1, np.where(g == 2, 1, 0)).astype(np.int8)   # hap A hom-alt
    b_hap = np.where(g < 0, -1, np.where(g >= 1, 1, 0)).astype(np.int8)   # hap B carries alt
    Ghap = np.empty((nsnp, ndip * 2), dtype=np.int8)
    Ghap[:, 0::2] = a_hap.T; Ghap[:, 1::2] = b_hap.T
    nhap = ndip * 2
    step = a.win // (a.numSubWins)     # 100kb
    offsets = list(range(0, a.L - a.win + 1, step))
    records = [ms_export.window_from(Ghap, pos, off, a.win) for off in offsets]
    ms = os.path.join(a.outdir, a.tag + ".sim.msOut.gz")
    ms_export.write_ms_gz(records, ms, nhap, seed=a.seed)
    fvec = os.path.join(a.outdir, a.tag + ".sim.fvec")
    preds = os.path.join(a.outdir, a.tag + ".sim.preds")
    subprocess.run([DE, "fvecSim", "diploid", ms, fvec,
                    "--totalPhysLen", str(a.win), "--numSubWins", str(a.numSubWins),
                    "--maskFileName", TRAIN_MASKFA, "--chrArmsForMasking", "1",
                    "--vcfForMaskFileName", TRAIN_VCFMASK, "--popForMask", "illex",
                    "--sampleToPopFileName", a.s2p,
                    "--unmaskedGenoFracCutoff", "0.5", "--unmaskedFracCutoff", "0.25"],
                   check=True)
    subprocess.run([DE, "predict", a.model + ".json", a.model + ".weights.h5",
                    fvec, preds, "--numSubWins", str(a.numSubWins)], check=True)
    return preds


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--trees", required=True)
    ap.add_argument("--chrom", required=True)
    ap.add_argument("--start", type=int, required=True)
    ap.add_argument("--L", type=int, required=True)
    ap.add_argument("--xpos", type=int, default=None, help="true sweep centre bp (default L/2)")
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--path", choices=["vcf", "sim", "both"], default="vcf")
    ap.add_argument("--no-missing", action="store_true", help="skip missingness (parity)")
    ap.add_argument("--maf-filter", action="store_true",
                    help="apply training MAF>0.01 cut (match the MAF-trained model)")
    ap.add_argument("--win", type=int, default=1_100_000)
    ap.add_argument("--numSubWins", type=int, default=11)
    ap.add_argument("--s2p", default=S2P_DEF)
    ap.add_argument("--model", default=MODEL_DEF)
    ap.add_argument("--realmask", default=REALMASK_DEF)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--tag", default=None)
    a = ap.parse_args()
    os.makedirs(a.outdir, exist_ok=True)
    if a.tag is None:
        a.tag = os.path.basename(a.trees).replace(".recap.trees", "")
    if a.xpos is None:
        a.xpos = a.L // 2
    chrom_name = "seg"
    samples = [s.split("\t")[0] for s in open(a.s2p) if s.strip()]

    ts = tskit.load(a.trees)
    g, pos = segment_ascertain(ts, a.seed, bake_missing=not a.no_missing,
                               maf_filter=a.maf_filter)
    print(f"[scanB] {a.tag}: {ts.num_samples} haps, seg L={a.L} xpos={a.xpos} "
          f"ascertained {g.shape[1]} SNPs (missing={'no' if a.no_missing else 'yes'}, "
          f"maf={'>0.01' if a.maf_filter else 'full'})", flush=True)

    outs = {}
    if a.path in ("vcf", "both"):
        outs["vcf"] = run_vcf_path(g, pos, a, samples, chrom_name)
    if a.path in ("sim", "both"):
        outs["sim"] = run_sim_path(g, pos, a)
    for k, p in outs.items():
        n = sum(1 for _ in open(p)) - 1
        print(f"[scanB] {k} preds -> {p} ({n} windows)", flush=True)


if __name__ == "__main__":
    main()
