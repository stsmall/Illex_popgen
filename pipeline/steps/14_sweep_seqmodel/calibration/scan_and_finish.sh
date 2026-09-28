#!/usr/bin/env bash
# Scan existing calibration trees at HIGH concurrency (compute-bound, parallelizes well),
# then fit HMM/HSMM params + decode the empirical genome -> Markov-corrected regions.tsv.
# Writes FINISH_DONE / FINISH_FAIL. nohup-friendly.
set -uo pipefail
D=/sietch_colab/data_share/illex/popgen_data/analysis/steps/14_sweep_seqmodel
CAL=$D/calibration
PY=/home/ssmall/miniforge3/envs/slim_sim/bin/python
LOG=$CAL/finish.log
OUT=$D/results/empirical_scan_fullsfs/hmm_decode
WIN=$D/results/empirical_scan_fullsfs/outlier_scan_45/windows.tsv
MAP=$D/config/mapdf.tsv
L=2500000
CONC=${1:-17}
exec >>"$LOG" 2>&1
rm -f "$CAL/FINISH_DONE" "$CAL/FINISH_FAIL"
echo "[scanfin $(date +%FT%T)] scanning existing trees conc=$CONC L=$L"

# 1) scan all existing trees -> manifest
"$PY" "$CAL/scan_existing.py" --segs "$CAL/sims/segs" --L "$L" --conc "$CONC" \
    --manifest "$CAL/sims/manifest.tsv"
rc=$?
nok=$(tail -n +2 "$CAL/sims/manifest.tsv" 2>/dev/null | awk -F'\t' '$5=="OK"' | wc -l)
echo "[scanfin $(date +%FT%T)] scan done rc=$rc OK=$nok"
if [ "$nok" -lt 6 ]; then echo "[scanfin] FAIL: too few OK segments ($nok)"; touch "$CAL/FINISH_FAIL"; exit 1; fi

# 2) fit params
echo "[scanfin $(date +%FT%T)] fitting params..."
"$PY" "$CAL/fit_params.py" --manifest "$CAL/sims/manifest.tsv" --mapdf "$MAP" \
    --L "$L" --windows "$WIN" --out "$CAL"
rc=$?
if [ $rc -ne 0 ] || [ ! -s "$CAL/params.json" ]; then
    echo "[scanfin] FAIL: fit_params rc=$rc"; touch "$CAL/FINISH_FAIL"; exit 1
fi

# 3) decode empirical (HMM + HSMM)
mkdir -p "$OUT"
echo "[scanfin $(date +%FT%T)] decoding empirical genome..."
"$PY" "$CAL/decode_empirical.py" --windows "$WIN" --mapdf "$MAP" \
    --params "$CAL/params.json" --model both --out "$OUT"
rc=$?
if [ $rc -ne 0 ] || [ ! -s "$OUT/regions.tsv" ]; then
    echo "[scanfin] FAIL: decode rc=$rc"; touch "$CAL/FINISH_FAIL"; exit 1
fi
echo "[scanfin $(date +%FT%T)] DONE."
wc -l "$OUT/regions.tsv"
touch "$CAL/FINISH_DONE"
