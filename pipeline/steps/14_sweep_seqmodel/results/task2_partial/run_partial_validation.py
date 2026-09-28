#!/usr/bin/env python
"""Task 2 validation: PARTIAL-sweep mode in the SLiM recipes.

Runs, into per-cell subdirs of results/task2_partial/:
  (A) REGRESSION cells (MODE=completed) that must reproduce the committed
      slim_validation / slim_validation_soft numbers (completed path unchanged);
  (B) PARTIAL cells (hard + soft) across SEL spanning the prior, which must
      terminate with present beneficial freq in [F_LO,F_HI], single m2 origin,
      n=700 reconstruct, flank pi~0.009 / TajD~-2.

Each cell -> run_slim_sweep.py into its own dir (filenames keyed only by seed, so
same-seed cells MUST use separate dirs). <=4 concurrent, thread-capped. After each
run the RAW SLiM trees (pre-overlay) are loaded to count m2 origins (must be 1).
Aggregates per-cell summaries into partial_validation_summary.{json,csv}.
"""
import csv
import json
import os
import subprocess
from concurrent.futures import ProcessPoolExecutor, as_completed

PY = "/home/ssmall/miniforge3/envs/slim_sim/bin/python"
DRIVER = "/sietch_colab/data_share/illex/popgen_data/analysis/steps/14_sweep_seqmodel/scripts/harness/slim/run_slim_sweep.py"
ROOT = "/sietch_colab/data_share/illex/popgen_data/analysis/steps/14_sweep_seqmodel/results/task2_partial"

# (tag, sweep_type, SEL, F0, seed, MODE)
CELLS = [
    # (A) regression: completed path must be unchanged
    ("reg_hard_completed", "hard", 0.01, None, 42, "completed"),
    ("reg_soft_completed", "soft", 0.01, 0.05, 42, "completed"),
    # (B) partial hard: SEL spanning the prior (+ a second seed for the middle)
    ("hard_partial_sel0.001", "hard", 0.001, None, 42, "partial"),
    ("hard_partial_sel0.01",  "hard", 0.01,  None, 42, "partial"),
    ("hard_partial_sel0.01_s43", "hard", 0.01, None, 43, "partial"),
    ("hard_partial_sel0.05",  "hard", 0.05,  None, 42, "partial"),
    # (B) partial soft: SEL spanning the prior at F0=0.05 (+ a second seed)
    ("soft_partial_sel0.001", "soft", 0.001, 0.05, 42, "partial"),
    ("soft_partial_sel0.01",  "soft", 0.01,  0.05, 42, "partial"),
    ("soft_partial_sel0.01_s43", "soft", 0.01, 0.05, 43, "partial"),
    ("soft_partial_sel0.05",  "soft", 0.05,  0.05, 42, "partial"),
]


def count_m2_origins(raw_trees):
    """Number of forward (m2) mutation origins in the RAW SLiM trees (pre-overlay).
    Forward mutation rate is 0, so every mutation in the raw ts is an m2 origin;
    a single-origin sweep => exactly 1."""
    try:
        import tskit
        ts = tskit.load(raw_trees)
        return int(ts.num_mutations)
    except Exception as e:  # pragma: no cover
        return f"ERR:{e}"


def run_cell(cell):
    tag, stype, sel, f0, seed, mode = cell
    outdir = os.path.join(ROOT, tag)
    os.makedirs(outdir, exist_ok=True)
    env = dict(os.environ)
    env.update(OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1", MKL_NUM_THREADS="1",
               NUMEXPR_NUM_THREADS="1")
    cmd = [PY, DRIVER, "--sweep-type", stype, "--seed", str(seed),
           "--SEL", repr(sel), "--MODE", mode, "--outdir", outdir, "--keep-slim-trees"]
    if f0 is not None:
        cmd += ["--F0", repr(f0)]
    logpath = os.path.join(outdir, "run.log")
    with open(logpath, "w") as log:
        log.write("CMD: " + " ".join(cmd) + "\n\n")
        log.flush()
        rc = subprocess.run(cmd, stdout=log, stderr=subprocess.STDOUT, env=env).returncode
    # count m2 origins in the raw slim trees
    raw = os.path.join(outdir, f"slim_raw_seed{seed}.trees")
    m2 = count_m2_origins(raw) if os.path.exists(raw) else "NO_RAW"
    return tag, rc, logpath, m2


def main():
    print(f"[val] {len(CELLS)} cells; 4 concurrent", flush=True)
    m2counts = {}
    rcs = {}
    with ProcessPoolExecutor(max_workers=4) as ex:
        futs = {ex.submit(run_cell, c): c for c in CELLS}
        for fut in as_completed(futs):
            tag, rc, logpath, m2 = fut.result()
            rcs[tag] = rc
            m2counts[tag] = m2
            print(f"[val] DONE {tag} rc={rc} m2_origins={m2} (log {logpath})", flush=True)

    rows = []
    for tag, stype, sel, f0, seed, mode in CELLS:
        spath = os.path.join(ROOT, tag, f"summary_seed{seed}.json")
        base = {"cell": tag, "sweep_type": stype, "SEL": sel, "F0": f0, "seed": seed,
                "MODE": mode, "rc": rcs.get(tag), "m2_origins": m2counts.get(tag)}
        if not os.path.exists(spath):
            base["status"] = "NO_SUMMARY"
            rows.append(base)
            continue
        s = json.load(open(spath))
        sw = s.get("sweep", {})
        base.update({
            "status": sw.get("status"),
            "present_freq": s.get("present_freq") or sw.get("freq"),
            "onset_tick": s.get("onset_tick") or sw.get("onset_tick"),
            "f0_realized": sw.get("f0_realized"),
            "t_est": sw.get("t_est"),
            "nrestart": sw.get("nrestart"),
            "n_haplotypes": s.get("n_haplotypes"),
            "n_sites": s.get("n_sites"),
            "flank_pi": s.get("flank_pi"),
            "flank_tajd": s.get("flank_tajd"),
            "center_pi": s.get("center_pi"),
            "center_core_pi": s.get("center_core_pi"),
            "ratio_center_flank": s.get("ratio_center_flank"),
            "ratio_core_flank": s.get("ratio_core_flank"),
            "core_center_hapH1": s.get("core_center_hapH1"),
            "core_flank_hapH1": s.get("core_flank_hapH1"),
            "ratio_core_hapH1_center_flank": s.get("ratio_core_hapH1_center_flank"),
        })
        rows.append(base)

    with open(os.path.join(ROOT, "partial_validation_summary.json"), "w") as fh:
        json.dump(rows, fh, indent=2)
    cols = sorted({k for r in rows for k in r})
    # stable, readable column order
    lead = ["cell", "sweep_type", "MODE", "SEL", "F0", "seed", "status",
            "present_freq", "onset_tick", "f0_realized", "t_est", "nrestart",
            "m2_origins", "n_haplotypes", "n_sites", "flank_pi", "flank_tajd",
            "center_pi", "center_core_pi", "ratio_center_flank", "ratio_core_flank",
            "core_center_hapH1", "core_flank_hapH1", "ratio_core_hapH1_center_flank", "rc"]
    cols = lead + [c for c in cols if c not in lead]
    with open(os.path.join(ROOT, "partial_validation_summary.csv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols)
        w.writeheader()
        w.writerows(rows)
    print("[val] wrote partial_validation_summary.json / .csv", flush=True)
    print(json.dumps(rows, indent=2), flush=True)


if __name__ == "__main__":
    main()
