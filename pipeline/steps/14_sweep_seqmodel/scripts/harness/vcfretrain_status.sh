#!/usr/bin/env bash
# Status for the full-SFS fvecVcf retrain. Usage: vcfretrain_status.sh [tag]
TAG=${1:-vcfretrain_fullsfs}
D=/sietch_colab/data_share/illex/popgen_data/analysis/steps/14_sweep_seqmodel
FR=$D/results/$TAG
NJOB=$(wc -l < "$FR/jobs.txt" 2>/dev/null || echo '?')
DONE=$(ls "$FR/parts"/*.fvec 2>/dev/null | grep -vc '\.tmp$')
FAIL=$(grep -l -e '^FAIL' "$FR/parts"/*.log 2>/dev/null | wc -l)
# running worker procs (comm-filter, not pgrep -f)
RUN=$(ps -eo comm,args | awk '/ms_to_vcf_fvec\.py/ && !/awk/' | wc -l)
echo "=== $TAG  $(date -u +%FT%TZ) ==="
echo "parts: $DONE/$NJOB done   failing-logs: $FAIL   worker-procs: $RUN"
tail -4 "$FR/progress" 2>/dev/null | sed 's/^/  /'
# if PARTS_DONE reached, show train phase
grep -qE 'RETRAIN_ALL_DONE|TRAIN_FAIL|FVEC_INCOMPLETE|MAKETS_FAIL' "$FR/progress" 2>/dev/null \
  && echo ">>> TERMINAL STATE REACHED"
