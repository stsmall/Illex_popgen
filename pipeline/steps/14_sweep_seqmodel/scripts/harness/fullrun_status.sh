#!/usr/bin/env bash
# Snapshot of the box diploSHIC chain (run_all.sh). Used by the monitor cron.
FR=/sietch_colab/data_share/illex/popgen_data/analysis/steps/14_sweep_seqmodel/results/fullrun
P=$FR/run_all.progress
echo "=== progress (last 10) ==="; tail -10 "$P" 2>/dev/null
alive=$(ps -eo args 2>/dev/null | grep -c '[r]un_all.sh')
nchunk=$(ls "$FR"/rawfvec/*.fvec 2>/dev/null | wc -l)
nlock=$(ls -d "$FR"/rawfvec/*.lock 2>/dev/null | wc -l)
if   grep -q RUNALL_DONE "$P" 2>/dev/null; then st=DONE
elif grep -qE 'FAIL|abort|INCOMPLETE' "$P" 2>/dev/null; then st=FAIL
elif grep -q TRAIN_START "$P" 2>/dev/null; then st=TRAINING
elif grep -q MAKETS_START "$P" 2>/dev/null; then st=MAKETS
elif grep -q FVEC_CHUNKS_DONE "$P" 2>/dev/null; then st=MERGING
elif grep -q FVEC_START "$P" 2>/dev/null; then st="FVEC(${nchunk}/46 chunks, ${nlock} active)"
else st=STARTING; fi
echo "STATE=$st  run_all_alive=$alive"
free -g | awk 'NR==2{print "mem: used="$3"G avail="$7"G"}'
# surface any failing chunk logs
fails=$(grep -lE 'Error|Traceback|Killed' "$FR"/fveclogs/*.log 2>/dev/null | wc -l)
echo "chunk logs with errors: $fails"
# active-chunk rep progress (min/max of the true 'starting rep N of M' line)
if [ -d "$FR/fveclogs" ]; then
  for l in "$FR"/fveclogs/*.log; do
    b=$(basename "$l" .log); [ -s "$FR/rawfvec/$b.fvec" ] && continue
    grep 'starting rep' "$l" 2>/dev/null | tail -1 | grep -oE 'rep [0-9]+ of [0-9]+'
  done 2>/dev/null | awk '{d=$2;t=$4} {sd+=d; st+=t; n++; if(d/t<mn||n==1)mn=d/t; if(d/t>mx)mx=d/t}
    END{if(n)printf "active chunks: %d, rep progress %.0f%%-%.0f%% (mean %.0f%%)\n",n,100*mn,100*mx,100*sd/st}'
fi
# train accuracy if trained
[ -f "$FR/train.log" ] && grep -iE 'accuracy|val_acc|confusion|test' "$FR/train.log" 2>/dev/null | tail -4 | sed 's/^/  train: /'
