#!/usr/bin/env python
"""Recompute stats (incl. the corrected NARROW-CORE haplotype homozygosity) on each
grid cell's SAVED sweep_final tree -- no SLiM re-run -- rewriting summary + CSV, then
rebuild grid_summary.{json,csv}. Uses the updated run_slim_sweep.emit_stats so the
output structure matches a fresh driver run exactly."""
import csv
import json
import os
import types

import tskit

import run_slim_sweep as R

ROOT = "/sietch_colab/data_share/illex/popgen_data/analysis/steps/14_sweep_seqmodel/results/slim_validation_soft"
CELLS = [(0.001, 0.02, 42), (0.001, 0.10, 42), (0.01, 0.02, 42), (0.01, 0.10, 42),
         (0.05, 0.02, 42), (0.05, 0.10, 42), (0.01, 0.10, 43)]


def cellname(sel, f0, seed):
    return f"sel{sel:g}_f{f0:g}_s{seed}"


def mk_args(sel, f0, seed):
    return types.SimpleNamespace(
        seed=seed, Q=100.0, N0=6808096.0, NREF=547928.0, TGROW=769519.0,
        MU=3e-9, R=2.1e-9, L=1100000, SEL=sel, TAU_GEN=20000.0,
        sweep_type="soft", F0=f0, FIX_THRESH=0.95, n_dip=350, n_windows=21,
        hap_core_bp=2000)


def main():
    rows = []
    for sel, f0, seed in CELLS:
        name = cellname(sel, f0, seed)
        outdir = os.path.join(ROOT, name)
        spath = os.path.join(outdir, f"summary_seed{seed}.json")
        tpath = os.path.join(outdir, f"sweep_final_seed{seed}.trees")
        if not (os.path.exists(spath) and os.path.exists(tpath)):
            print(f"[recompute] MISSING {name} (summary or tree); skipping")
            continue
        old = json.load(open(spath))
        sweep_info = old["sweep"]                     # keep SLiM's SWEEP_RESULT tokens
        ts = tskit.load(tpath)
        args = mk_args(sel, f0, seed)
        stats_csv = os.path.join(outdir, f"windowed_stats_seed{seed}.csv")
        summary = R.emit_stats(args, ts, sweep_info, tpath, stats_csv, spath)
        rows.append({
            "cell": name, "SEL": sel, "F0": f0, "seed": seed,
            "status": sweep_info.get("status"), "freq": sweep_info.get("freq"),
            "f0_realized": sweep_info.get("f0_realized"),
            "allele_age_ticks": sweep_info.get("allele_age_ticks"),
            "t_est": sweep_info.get("t_est"), "nrestart": sweep_info.get("nrestart"),
            "n_sites": summary["n_sites"],
            "flank_pi": round(summary["flank_pi"], 6), "flank_tajd": round(summary["flank_tajd"], 4),
            "center_pi": round(summary["center_pi"], 6),
            "ratio_center_flank": round(summary["ratio_center_flank"], 4),
            "ratio_core_flank": round(summary["ratio_core_flank"], 4),
            "core_center_hapH1": round(summary["core_center_hapH1"], 4),
            "core_flank_hapH1": round(summary["core_flank_hapH1"], 4),
            "ratio_core_hapH1": round(summary["ratio_core_hapH1_center_flank"], 3),
            "core_center_hapH12": round(summary["core_center_hapH12"], 4),
            "core_center_hap_ndistinct": summary["core_center_hap_ndistinct"],
            "core_flank_hap_ndistinct": round(summary["core_flank_hap_ndistinct"], 1),
        })
    with open(os.path.join(ROOT, "grid_summary.json"), "w") as fh:
        json.dump(rows, fh, indent=2)
    if rows:
        with open(os.path.join(ROOT, "grid_summary.csv"), "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
            w.writeheader()
            w.writerows(rows)
    print("\n[recompute] rebuilt grid_summary.{json,csv}\n")
    print(json.dumps(rows, indent=2))


if __name__ == "__main__":
    main()
