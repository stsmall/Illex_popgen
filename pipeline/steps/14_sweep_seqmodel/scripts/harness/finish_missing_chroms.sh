#!/usr/bin/env bash
# Detached finisher: wait for the chr2 + chr42 full-SFS fvecs, then run predict_fullsfs.sh
# (predicts only the new fvecs -- skips the 43 existing preds -- and rebuilds genome.preds to 45 chroms).
# Real detached process (setsid), polls on disk; not a harness background-waiter.
set -uo pipefail
D=/sietch_colab/data_share/illex/popgen_data/analysis/steps/14_sweep_seqmodel
OUT=$D/results/empirical_scan_fullsfs; FVD=$OUT/fvec
LOG=$OUT/finish_missing.log; : > "$LOG"
say(){ echo "$(date -u +%FT%TZ) $*" >> "$LOG"; }
say "FINISHER_START waiting for chr2.fvec + chr42.fvec"
# wait up to 6h for both fvecs; bail if the fvec jobs died without producing them
for i in $(seq 1 720); do
  n2=$([ -s "$FVD/chr2.fvec" ] && echo 1 || echo 0)
  n42=$([ -s "$FVD/chr42.fvec" ] && echo 1 || echo 0)
  if [ "$n2" = 1 ] && [ "$n42" = 1 ]; then say "BOTH_FVECS_PRESENT"; break; fi
  # if the fvec driver is gone and fvecs still missing, keep waiting a bit (subset is long) but log
  [ $((i % 12)) -eq 0 ] && say "waiting... chr2.fvec=$n2 chr42.fvec=$n42 (tick $i)"
  sleep 30
done
if [ ! -s "$FVD/chr2.fvec" ] || [ ! -s "$FVD/chr42.fvec" ]; then
  say "TIMEOUT_OR_MISSING chr2=$([ -s "$FVD/chr2.fvec" ]&&echo Y||echo N) chr42=$([ -s "$FVD/chr42.fvec" ]&&echo Y||echo N); NOT running predict"
  exit 1
fi
say "RUN_PREDICT (rebuilds genome.preds)"
bash "$D/scripts/harness/predict_fullsfs.sh" >> "$LOG" 2>&1
np=$(tail -n +2 "$OUT/genome.preds" 2>/dev/null | cut -f1 | sort -u | wc -l)
say "PREDICT_DONE genome.preds now spans $np chroms; chr2 win=$(tail -n +2 "$OUT/genome.preds"|awk -F'\t' '$1=="2"'|wc -l) chr42 win=$(tail -n +2 "$OUT/genome.preds"|awk -F'\t' '$1=="42"'|wc -l)"
say "FINISHER_DONE"
