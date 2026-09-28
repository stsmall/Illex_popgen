#!/usr/bin/env bash
# Wait for the running gen+scan campaign (run_pilot_campaign.py) to finish, then
# fit the HMM/HSMM params and decode the empirical genome -> Markov-corrected regions.
# Idempotent: safe to re-run; writes FINISH_DONE / FINISH_FAIL marker. nohup-friendly.
set -uo pipefail
D=/sietch_colab/data_share/illex/popgen_data/analysis/steps/14_sweep_seqmodel
CAL=$D/calibration
PY=/home/ssmall/miniforge3/envs/slim_sim/bin/python   # same env the driver uses
LOG=$CAL/finish.log
OUT=$D/results/empirical_scan_fullsfs/hmm_decode
WIN=$D/results/empirical_scan_fullsfs/outlier_scan_45/windows.tsv
MAP=$D/config/mapdf.tsv
L=2500000                       # must match the campaign --L
DRIVER_PID=${1:-1841881}
exec >>"$LOG" 2>&1
echo "[finish $(date +%FT%T)] waiting on driver PID $DRIVER_PID (scan phase)"

# 1) wait for the driver to exit (scan phase complete)
while kill -0 "$DRIVER_PID" 2>/dev/null; do sleep 120; done
tot=$(tail -n +2 "$CAL/sims/manifest.tsv" 2>/dev/null | wc -l)
nok=$(grep -c $'\tOK\t' "$CAL/sims/manifest.tsv" 2>/dev/null || echo 0)
echo "[finish $(date +%FT%T)] driver exited; manifest rows=$tot scanned-OK(tab)=$nok"
if [ "$tot" -eq 0 ]; then echo "[finish] FAIL: empty manifest"; touch "$CAL/FINISH_FAIL"; exit 1; fi

# 2) fit params (writes params.json + null_maxima_*.npy INTO $CAL)
echo "[finish $(date +%FT%T)] fitting params (L=$L)..."
"$PY" "$CAL/fit_params.py" --manifest "$CAL/sims/manifest.tsv" --mapdf "$MAP" \
    --L "$L" --windows "$WIN" --out "$CAL"
rc=$?
if [ $rc -ne 0 ] || [ ! -s "$CAL/params.json" ]; then
    echo "[finish] FAIL: fit_params rc=$rc"; touch "$CAL/FINISH_FAIL"; exit 1
fi
echo "[finish] params.json written"

# 3) decode empirical (HMM + HSMM) -> regions.tsv
mkdir -p "$OUT"
echo "[finish $(date +%FT%T)] decoding empirical genome (both models)..."
"$PY" "$CAL/decode_empirical.py" --windows "$WIN" --mapdf "$MAP" \
    --params "$CAL/params.json" --model both --out "$OUT"
rc=$?
if [ $rc -ne 0 ] || [ ! -s "$OUT/regions.tsv" ]; then
    echo "[finish] FAIL: decode rc=$rc"; touch "$CAL/FINISH_FAIL"; exit 1
fi
echo "[finish $(date +%FT%T)] DONE."
wc -l "$OUT/regions.tsv"; head -1 "$OUT/regions.tsv"
touch "$CAL/FINISH_DONE"
