#!/usr/bin/env python3
"""fit_params.py -- PLAN_02 Phase C+D: from the scanned calibration segments (manifest.tsv),
LABEL each window N/LL/C/LR by genetic distance from the known sweep centre, then FIT the
recomb-aware decoder's parameters and the neutral null threshold. Writes calibration/params.json
+ null_maxima_{hmm,hsmm}.npy with full provenance. No fabricated values: every parameter is
derived from the labeled sims or a documented genome constant.

Fitting recipe (all data-driven):
  d_C, d_L  : genetic distances where the pooled diploSHIC signal crosses over
              (pC=pL -> d_C ; pL=pN -> d_L), from hard+soft segments binned by |gdist|.
  emission  : transform chosen (prior_corrected vs raw) by labeled-state argmax accuracy;
              priors fixed to the balanced 5-class training prior (0.2,0.4,0.4).
  HMM lambdas: lambda_c=d_C ; lambda_ll=lambda_lr=(d_L-d_C) ; lambda_n=large (neutral dwell).
  HMM topo  : narrow_entry/abort_* from empirical labeled-run transition frequencies.
  beta      : genome length-weighted mean recomb rate (physical fallback for zero-rho).
  HSMM dur  : lognormal(median,sigma) of the genetic SPAN of C/LL/LR runs across segments.
  null      : neutral segments decoded -> per-segment max region score -> block-bootstrap the
              genome-wide max ; threshold = (1-alpha) quantile.
"""
import argparse, json, math, os, sys, time
import numpy as np, pandas as pd

D = "/sietch_colab/data_share/illex/popgen_data/analysis/steps/14_sweep_seqmodel"
sys.path.insert(0, f"{D}/scripts")
from seqmodel.std_diploshic import standardize
from seqmodel.emissions import collapse, emission_logprob
from seqmodel.calib_label import seg_mapdf, label_windows
from seqmodel.mapcoord import split_blocks, deltas, genetic_positions
from seqmodel.decode import viterbi_hmm, viterbi_hsmm
from seqmodel.regions import extract_regions
from seqmodel.calibrate_null import genome_max, threshold as null_threshold
from seqmodel import STATE_IDX

SCHEMA = f"{D}/config/probability_schema.yaml"
MEAN_RATE = 2.51859e-09          # genome length-weighted mean ReLERNN rate (M/bp)
PRIORS = (0.2, 0.4, 0.4)         # balanced 5-class training prior, collapsed (N,L,C)
EPSF = 1e-9


def load_seg(row, genome_map, L):
    """Standardized preds + seg-coord map + signed genetic distance (cM) per window."""
    preds = standardize(row["preds"], SCHEMA)          # chrom='seg', local coords
    smap = seg_mapdf(genome_map, "43", int(row["start"]), L, seg_name="seg")
    xpos = int(row.get("xpos", L // 2))
    lab = label_windows(preds, smap, xpos, d_C=1.0, d_L=1.0)   # d huge -> only need gdist_cM
    lab = lab[np.isfinite(lab["gdist_cM"])].reset_index(drop=True)
    lab["pC"] = lab["p_hard"] + lab["p_soft"]
    lab["pL"] = lab["p_linked_hard"] + lab["p_linked_soft"]
    lab["pN"] = lab["p_neutral"]
    return preds, smap, lab, xpos


def crossover(x, ya, yb):
    """smallest x>0 where mean(ya) drops below mean(yb), via binned interpolation."""
    order = np.argsort(x)
    x, ya, yb = x[order], ya[order], yb[order]
    diff = ya - yb
    below = np.where(diff < 0)[0]
    if len(below) == 0:
        return float(x.max())
    i = below[0]
    if i == 0:
        return float(x[0])
    # linear interp of diff crossing 0 between i-1 and i
    x0, x1, d0, d1 = x[i - 1], x[i], diff[i - 1], diff[i]
    if d1 == d0:
        return float(x1)
    return float(x0 + (0 - d0) * (x1 - x0) / (d1 - d0))


def fit_labels(sweeps):
    """Pool hard+soft; return d_C, d_L in Morgans + the binned profile for provenance."""
    g = np.concatenate([abs(s["lab"]["gdist_cM"].to_numpy()) / 100.0 for s in sweeps])  # Morgans
    pC = np.concatenate([s["lab"]["pC"].to_numpy() for s in sweeps])
    pL = np.concatenate([s["lab"]["pL"].to_numpy() for s in sweeps])
    pN = np.concatenate([s["lab"]["pN"].to_numpy() for s in sweeps])
    # bin by |gdist| (Morgans)
    nb = 40
    edges = np.linspace(0, np.quantile(g, 0.99), nb + 1)
    mid = 0.5 * (edges[:-1] + edges[1:])
    bpC = np.array([pC[(g >= edges[i]) & (g < edges[i + 1])].mean() if ((g >= edges[i]) & (g < edges[i + 1])).any() else np.nan for i in range(nb)])
    bpL = np.array([pL[(g >= edges[i]) & (g < edges[i + 1])].mean() if ((g >= edges[i]) & (g < edges[i + 1])).any() else np.nan for i in range(nb)])
    bpN = np.array([pN[(g >= edges[i]) & (g < edges[i + 1])].mean() if ((g >= edges[i]) & (g < edges[i + 1])).any() else np.nan for i in range(nb)])
    ok = np.isfinite(bpC) & np.isfinite(bpL) & np.isfinite(bpN)
    mid, bpC, bpL, bpN = mid[ok], bpC[ok], bpL[ok], bpN[ok]
    d_C = crossover(mid, bpC, bpL)             # center dominates within d_C
    d_L = crossover(mid, bpL, bpN)             # linked dominates within d_L
    d_L = max(d_L, d_C * 1.5)
    profile = [dict(gdist_M=float(m), pC=float(c), pL=float(l), pN=float(n))
               for m, c, l, n in zip(mid, bpC, bpL, bpN)]
    return d_C, d_L, profile


def relabel(sweeps, d_C, d_L):
    for s in sweeps:
        lab = s["lab"].copy()
        dM = np.abs(lab["gdist_cM"].to_numpy()) / 100.0
        signed = lab["gdist_cM"].to_numpy()
        reg = np.where(dM <= d_C, "C", np.where(dM <= d_L, np.where(signed < 0, "LL", "LR"), "N"))
        lab["region"] = reg
        s["lab2"] = lab.reset_index(drop=True)


def emission_accuracy(sweeps, transform):
    """argmax(emission state) vs true label, collapsing LL/LR->L; over labeled windows."""
    order = {"N": 0, "L": 1, "C": 2}
    correct = tot = 0
    for s in sweeps:
        lab = s["lab2"]
        col = collapse(lab)
        e = emission_logprob(col, transform=transform, priors=PRIORS)  # (n,4) N,LL,C,LR
        # collapse emission to N,L,C (LL,LR share L)
        e3 = np.column_stack([e[:, 0], e[:, 1], e[:, 2]])
        pred = e3.argmax(1)
        true = lab["region"].map({"N": 0, "LL": 1, "LR": 1, "C": 2}).to_numpy()
        m = np.isin(true, [0, 1, 2])
        correct += (pred[m] == true[m]).sum(); tot += m.sum()
    return correct / max(tot, 1)


def run_lengths_genetic(sweeps, state_set):
    """genetic spans (Morgans) of contiguous runs of labeled states in `state_set`."""
    spans = []
    for s in sweeps:
        lab = s["lab2"]
        gpos = genetic_positions(lab, s["smap"], "seg")
        reg = lab["region"].to_numpy()
        i = 0
        while i < len(reg):
            if reg[i] in state_set:
                j = i
                while j + 1 < len(reg) and reg[j + 1] in state_set:
                    j += 1
                span = float(gpos[j] - gpos[i]) if j > i else float((lab["end"].iloc[i] - lab["start"].iloc[i]) * MEAN_RATE)
                spans.append(max(span, EPSF))
                i = j + 1
            else:
                i += 1
    return np.array(spans)


def lognorm_fit(spans, floor_med):
    if len(spans) < 2:
        return float(floor_med), 0.6
    ls = np.log(np.clip(spans, EPSF, None))
    return float(np.exp(ls.mean())), float(max(ls.std(ddof=1), 0.2))


def transition_freqs(sweeps):
    """empirical labeled-run transition counts to estimate narrow_entry/abort_left/abort_center."""
    from collections import Counter
    c = Counter()
    for s in sweeps:
        reg = list(s["lab2"]["region"])
        runs = []
        i = 0
        while i < len(reg):
            j = i
            while j + 1 < len(reg) and reg[j + 1] == reg[i]:
                j += 1
            runs.append(reg[i]); i = j + 1
        for a, b in zip(runs[:-1], runs[1:]):
            c[(a, b)] += 1
    def frac(numer, denom_keys):
        den = sum(c[k] for k in denom_keys)
        return (c[numer] / den) if den else 0.0
    # narrow_entry: N->C among N-exits into sweep (N->C vs N->LL)
    narrow_entry = frac(("N", "C"), [("N", "C"), ("N", "LL")])
    # abort_left: LL->N among LL exits (LL->N vs LL->C)
    abort_left = frac(("LL", "N"), [("LL", "N"), ("LL", "C")])
    # abort_center: C->N among C exits (C->N vs C->LR)
    abort_center = frac(("C", "N"), [("C", "N"), ("C", "LR")])
    return dict(counts={f"{a}->{b}": n for (a, b), n in c.items()},
                narrow_entry=narrow_entry, abort_left=abort_left, abort_center=abort_center)


def build_params(d_C, d_L, transform, hsmm_dur, topo, max_dur_windows, mean_win_gap):
    lam_c = max(d_C, 1e-4)
    lam_flank = max(d_L - d_C, 5e-4)
    hmm = dict(lambda_n=max(10 * d_L, 0.05), lambda_ll=lam_flank, lambda_c=lam_c,
               lambda_lr=lam_flank, beta=MEAN_RATE,
               skip_weight=0.01,
               narrow_entry=float(np.clip(topo["narrow_entry"], 0.01, 0.5)),
               abort_left=float(np.clip(topo["abort_left"], 0.02, 0.6)),
               abort_center=float(np.clip(topo["abort_center"], 0.02, 0.6)))
    # HSMM penalties from canonical vs off-path frequency (log-odds, floored)
    ne = hmm["narrow_entry"]; al = hmm["abort_left"]; ac = hmm["abort_center"]
    hsmm = dict(ll_med=hsmm_dur["LL"][0], ll_sig=hsmm_dur["LL"][1],
                c_med=hsmm_dur["C"][0], c_sig=hsmm_dur["C"][1],
                lr_med=hsmm_dur["LR"][0], lr_sig=hsmm_dur["LR"][1],
                max_dur_windows=int(max_dur_windows),
                trans_penalty=float(math.log(max(1 - ne, 0.1)) - 0.0),
                narrow_penalty=float(math.log(max(ne, 0.01))),
                abort_penalty=float(math.log(max(min(al, ac), 0.01))))
    emission = dict(transform=transform, priors=list(PRIORS), temperature=1.0)
    decode = dict(max_phys_gap=300000, fallback_rate=MEAN_RATE)
    return dict(emission=emission, hmm=hmm, hsmm=hsmm, decode=decode,
                labels=dict(d_C_Morgans=d_C, d_L_Morgans=d_L))


def decode_neutral_maxima(neutrals, params, genome_map, L, model):
    """decode each neutral segment, return per-segment max region score."""
    maxima = []
    for row in neutrals:
        preds = standardize(row["preds"], SCHEMA)
        smap = seg_mapdf(genome_map, "43", int(row["start"]), L, seg_name="seg")
        col = collapse(preds)
        emis = emission_logprob(col, transform=params["emission"]["transform"],
                                priors=tuple(params["emission"]["priors"]))
        blocks = split_blocks(preds, smap, "seg", params["decode"]["max_phys_gap"])
        path = ["GAP"] * len(preds)
        for (a, b) in blocks:
            wb = preds.iloc[a:b + 1]
            if model == "hmm":
                dg, dx = deltas(wb, smap, "seg") if b > a else (np.array([]), np.array([]))
                p, _ = viterbi_hmm(emis[a:b + 1], dg, dx, params["hmm"])
            else:
                gp = genetic_positions(wb, smap, "seg")
                p, _ = viterbi_hsmm(emis[a:b + 1], gp - gp[0], params["hsmm"])
            path[a:b + 1] = p
        if model == "hmm":
            gf = genetic_positions(preds, smap, "seg"); dgf = np.diff(gf)
            dxf = np.diff(((preds["start"].to_numpy(np.int64) + preds["end"].to_numpy(np.int64)) // 2).astype(float))
            R = extract_regions(preds, path, emis, col, params=params["hmm"], model="hmm", dg=dgf, dx=dxf)
        else:
            gf = genetic_positions(preds, smap, "seg")
            R = extract_regions(preds, path, emis, col, params=params["hsmm"], model="hsmm", gpos=gf)
        maxima.append(genome_max(R))
    return np.array(maxima, float)


def bootstrap_genome_max(seg_maxima, n_blocks_per_genome, B, rng):
    out = np.empty(B)
    for b in range(B):
        draw = rng.choice(seg_maxima, size=n_blocks_per_genome, replace=True)
        out[b] = draw.max()
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--manifest", default=f"{D}/calibration/sims/manifest.tsv")
    ap.add_argument("--mapdf", default=f"{D}/config/mapdf.tsv")
    ap.add_argument("--L", type=int, default=3_000_000)
    ap.add_argument("--windows", default=f"{D}/results/empirical_scan_fullsfs/outlier_scan_45/windows.tsv")
    ap.add_argument("--alpha", type=float, default=0.05)
    ap.add_argument("--B", type=int, default=2000)
    ap.add_argument("--out", default=f"{D}/calibration")
    a = ap.parse_args()
    genome_map = pd.read_csv(a.mapdf, sep="\t"); genome_map["chrom"] = genome_map["chrom"].astype(str)
    man = pd.read_csv(a.manifest, sep="\t")
    man = man[man["status"] == "OK"].reset_index(drop=True)
    sweeps, neutrals = [], []
    for _, row in man.iterrows():
        if row["klass"] in ("hard", "soft"):
            preds, smap, lab, xpos = load_seg(row, genome_map, a.L)
            sweeps.append(dict(row=row, preds=preds, smap=smap, lab=lab, xpos=xpos, klass=row["klass"]))
        elif row["klass"] == "neutral":
            neutrals.append(row)
    print(f"[fit] {len(sweeps)} sweep segs (hard={sum(s['klass']=='hard' for s in sweeps)} "
          f"soft={sum(s['klass']=='soft' for s in sweeps)}), {len(neutrals)} neutral segs")

    # ---- d_C, d_L ----
    d_C, d_L, profile = fit_labels(sweeps)
    relabel(sweeps, d_C, d_L)
    print(f"[fit] d_C={d_C:.5f} M ({d_C*100:.3f} cM), d_L={d_L:.5f} M ({d_L*100:.3f} cM)")

    # ---- emission transform ----
    acc_pc = emission_accuracy(sweeps, "prior_corrected")
    acc_raw = emission_accuracy(sweeps, "raw")
    transform = "prior_corrected" if acc_pc >= acc_raw else "raw"
    print(f"[fit] emission argmax acc: prior_corrected={acc_pc:.3f} raw={acc_raw:.3f} -> {transform}")

    # ---- HSMM durations + max_dur_windows ----
    dur = {st: lognorm_fit(run_lengths_genetic(sweeps, {st}),
                           floor_med=(d_C if st == "C" else (d_L - d_C)))
           for st in ("C", "LL", "LR")}
    # max_dur_windows: max labeled run length in windows across C/LL/LR
    maxrun = 1
    for s in sweeps:
        reg = list(s["lab2"]["region"]); i = 0
        while i < len(reg):
            j = i
            while j + 1 < len(reg) and reg[j + 1] == reg[i]:
                j += 1
            if reg[i] in ("C", "LL", "LR"):
                maxrun = max(maxrun, j - i + 1)
            i = j + 1
    max_dur_windows = max(maxrun + 2, 6)
    print(f"[fit] durations (median M, sigma): " + ", ".join(f"{k}={v[0]:.4g}/{v[1]:.2g}" for k, v in dur.items())
          + f"  max_dur_windows={max_dur_windows}")

    # ---- topology freqs ----
    topo = transition_freqs(sweeps)
    print(f"[fit] topo: narrow_entry={topo['narrow_entry']:.3f} abort_left={topo['abort_left']:.3f} "
          f"abort_center={topo['abort_center']:.3f}")

    params = build_params(d_C, d_L, transform, dur, topo, max_dur_windows, mean_win_gap=None)

    # ---- null threshold (block-bootstrap genome max) ----
    win = pd.read_csv(a.windows, sep="\t"); win["chrom"] = win["chrom"].astype(str)
    win_auto = win[~win["chrom"].isin(["42"])]        # autosomal windows (exclude sex chr42)
    G = len(win_auto)
    # windows per neutral segment (median observed)
    seg_wins = []
    for row in neutrals:
        try:
            seg_wins.append(sum(1 for _ in open(row["preds"])) - 1)
        except Exception:
            pass
    w_per_seg = int(np.median(seg_wins)) if seg_wins else 30
    n_blocks = max(int(round(G / max(w_per_seg, 1))), 1)
    rng = np.random.default_rng(7)
    nullinfo = {}
    for model in ("hmm", "hsmm"):
        seg_max = decode_neutral_maxima(neutrals, params, genome_map, a.L, model)
        boot = bootstrap_genome_max(seg_max, n_blocks, a.B, rng)
        thr = null_threshold(boot, a.alpha)
        np.save(os.path.join(a.out, f"null_maxima_{model}.npy"), boot)
        nullinfo[f"threshold_{model}"] = float(thr)
        nullinfo[f"seg_max_{model}_nonzero"] = int((seg_max > 0).sum())
        nullinfo[f"seg_max_{model}_max"] = float(seg_max.max()) if len(seg_max) else 0.0
        print(f"[fit] null {model}: {len(seg_max)} neutral segs, {int((seg_max>0).sum())} with a region, "
              f"seg_max in [{seg_max.min():.2f},{seg_max.max():.2f}]; genome-max thr(alpha={a.alpha})={thr:.3f}")
    nullinfo.update(alpha=a.alpha, B=a.B, n_neutral_segments=len(neutrals),
                    windows_per_segment=w_per_seg, genome_autosome_windows=int(G),
                    n_blocks_per_genome=n_blocks, segment_L=a.L, model_primary="hmm")
    params["null"] = nullinfo

    params["provenance"] = dict(
        created=time.strftime("%Y-%m-%dT%H:%M:%S"),
        plan="PLAN_02 Phase C+D (pilot)", manifest=a.manifest,
        n_hard=int(sum(s['klass'] == 'hard' for s in sweeps)),
        n_soft=int(sum(s['klass'] == 'soft' for s in sweeps)),
        n_neutral=len(neutrals), segment_L=a.L, chrom="43",
        model_vcf=f"{D}/results/vcfretrain_fullsfs/illexModel_vcf",
        priors="balanced 5-class training (0.2,0.4,0.4)",
        mean_recomb_rate_Mbp=MEAN_RATE,
        sweep_sel_prior="10**U(-4,-2)", note="pilot-scale calibration; see report",
        label_profile=profile, topo_counts=topo["counts"])

    outp = os.path.join(a.out, "params.json")
    json.dump(params, open(outp, "w"), indent=2)
    print(f"[fit] wrote {outp}")


if __name__ == "__main__":
    main()
