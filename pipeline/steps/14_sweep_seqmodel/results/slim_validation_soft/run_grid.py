#!/usr/bin/env python
"""Orchestrate the SOFT-sweep validation grid over (SEL, F0) at n=700, Q=100.
Runs run_slim_sweep.py per cell into a per-cell outdir (filenames are keyed only by
seed, so same-seed cells MUST use separate dirs), <=4 concurrent, thread-capped.
Aggregates per-cell summary JSONs into grid_summary.{json,csv}."""
import csv
import json
import os
import subprocess
import sys
from concurrent.futures import ProcessPoolExecutor, as_completed

PY = "/home/ssmall/miniforge3/envs/slim_sim/bin/python"
DRIVER = "/sietch_colab/data_share/illex/popgen_data/analysis/steps/14_sweep_seqmodel/scripts/harness/slim/run_slim_sweep.py"
ROOT = "/sietch_colab/data_share/illex/popgen_data/analysis/steps/14_sweep_seqmodel/results/slim_validation_soft"

# (SEL, F0, seed): 3x2 grid at seed 42 + one intermediate cell repeated at seed 43
CELLS = []
for sel in (0.001, 0.01, 0.05):
    for f0 in (0.02, 0.10):
        CELLS.append((sel, f0, 42))
CELLS.append((0.01, 0.10, 43))   # extra seed for one (intermediate) cell


def cellname(sel, f0, seed):
    return f"sel{sel:g}_f{f0:g}_s{seed}"


def run_cell(cell):
    sel, f0, seed = cell
    name = cellname(sel, f0, seed)
    outdir = os.path.join(ROOT, name)
    os.makedirs(outdir, exist_ok=True)
    env = dict(os.environ)
    env.update(OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1", MKL_NUM_THREADS="1",
               NUMEXPR_NUM_THREADS="1")
    cmd = [PY, DRIVER, "--sweep-type", "soft", "--seed", str(seed),
           "--SEL", repr(sel), "--F0", repr(f0), "--outdir", outdir]
    logpath = os.path.join(outdir, "run.log")
    with open(logpath, "w") as log:
        log.write("CMD: " + " ".join(cmd) + "\n\n")
        log.flush()
        rc = subprocess.run(cmd, stdout=log, stderr=subprocess.STDOUT, env=env).returncode
    return name, rc, logpath


def main():
    print(f"[grid] {len(CELLS)} cells; 4 concurrent", flush=True)
    results = {}
    with ProcessPoolExecutor(max_workers=4) as ex:
        futs = {ex.submit(run_cell, c): c for c in CELLS}
        for fut in as_completed(futs):
            name, rc, logpath = fut.result()
            results[name] = rc
            print(f"[grid] DONE {name} rc={rc} (log {logpath})", flush=True)

    # aggregate summaries
    rows = []
    for sel, f0, seed in CELLS:
        name = cellname(sel, f0, seed)
        spath = os.path.join(ROOT, name, f"summary_seed{seed}.json")
        if not os.path.exists(spath):
            rows.append({"cell": name, "SEL": sel, "F0": f0, "seed": seed,
                         "status": "NO_SUMMARY", "rc": results.get(name)})
            continue
        s = json.load(open(spath))
        sw = s.get("sweep", {})
        rows.append({
            "cell": name, "SEL": sel, "F0": f0, "seed": seed,
            "status": sw.get("status"), "freq": sw.get("freq"),
            "f0_realized": sw.get("f0_realized"), "allele_age_ticks": sw.get("allele_age_ticks"),
            "t_est": sw.get("t_est"), "nrestart": sw.get("nrestart"),
            "n_sites": s.get("n_sites"), "n_haplotypes": s.get("n_haplotypes"),
            "flank_pi": s.get("flank_pi"), "flank_tajd": s.get("flank_tajd"),
            "center_pi": s.get("center_pi"),
            "ratio_center_flank": s.get("ratio_center_flank"),
            "ratio_core_flank": s.get("ratio_core_flank"),
            "center_hapH1": s.get("center_hapH1"), "flank_hapH1": s.get("flank_hapH1"),
            "center_hapH12": s.get("center_hapH12"), "flank_hapH12": s.get("flank_hapH12"),
            "center_hap_ndistinct": s.get("center_hap_ndistinct"),
            "flank_hap_ndistinct": s.get("flank_hap_ndistinct"),
            "ratio_hapH1_center_flank": s.get("ratio_hapH1_center_flank"),
            "ratio_hapH12_center_flank": s.get("ratio_hapH12_center_flank"),
        })
    with open(os.path.join(ROOT, "grid_summary.json"), "w") as fh:
        json.dump(rows, fh, indent=2)
    cols = list(rows[0].keys())
    with open(os.path.join(ROOT, "grid_summary.csv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols)
        w.writeheader()
        w.writerows(rows)
    print("[grid] wrote grid_summary.json / grid_summary.csv", flush=True)
    print(json.dumps(rows, indent=2), flush=True)


if __name__ == "__main__":
    main()
