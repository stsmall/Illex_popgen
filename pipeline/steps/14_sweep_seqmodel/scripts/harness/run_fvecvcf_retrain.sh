#!/usr/bin/env bash
# Regenerate diploSHIC TRAINING features via the SCAN path (fvecVcf on unphased VCFs) so training
# matches the empirical scan (fix diploshic-zns-train-scan-mismatch: fvecSim phased vs fvecVcf
# unphased; fvecVcf diploid IGNORES phase, so training MUST use fvecVcf). Splits each class-position
# into sub-jobs (<= MAXREPS reps) run in parallel to amortize fvecVcf's ~10s startup and beat the
# single-sequential-job bottleneck; concats parts; strips the 4 fvecVcf meta cols (train wants
# pure-numeric fvecSim format); makeTrainingSets + train -> illexModel_vcf.
# Run detached:  setsid bash run_fvecvcf_retrain.sh <tag> >/dev/null 2>&1 </dev/null &
set -uo pipefail
TAG=${1:-vcfretrain_full}
D=/sietch_colab/data_share/illex/popgen_data/analysis/steps/14_sweep_seqmodel
H=$D/scripts/harness
PY=/home/ssmall/miniforge3/envs/slim_sim/bin/python
DE=/home/ssmall/miniforge3/envs/diploshic_env/bin
BANK=${BANK:-$D/results/task17_maskval/geno_mask_bank.chr1.npz}   # structured (correlated) real missingness overlay
NDIP=${NDIP:-}   # optional: subsample each sim to this many diploids (chrZ-males retrain uses 171). Empty = all m//2 (=350).
# NOTE: if NDIP is set, BANK MUST be a matching n_indiv bank (e.g. the chrZ 171-sample bank) or fvecVcf/mask shapes mismatch.
FR=$D/results/$TAG; PARTS=$FR/parts; FV=$FR/fvec; FVN=$FR/fvec_nometa; TS=$FR/trainingSets; TMP=$FR/tmp
mkdir -p "$PARTS" "$FV" "$FVN" "$TS" "$TMP"
PROG=$FR/progress; : > "$PROG"
say(){ echo "$(date -u +%FT%TZ) $*" >> "$PROG"; }
# full-scale defaults (override via env)
NEUT_TASKS=${NEUT_TASKS:-$(seq 0 18)}          # 110/task -> ~2090 neutral
FOCAL_TASKS=${FOCAL_TASKS:-$(seq 0 199)}       # 10/task -> ~2000 focal hard_5/soft_5
LINKED_TASKS=${LINKED_TASKS:-$(seq 0 19)}      # 10/task -> ~200 per linked position (x10 -> 2000 pooled)
TASKS_PER_JOB=${TASKS_PER_JOB:-4}              # focal/linked: 4 tasks (~40 reps) per sub-job
NEUT_TASKS_PER_JOB=${NEUT_TASKS_PER_JOB:-1}    # neutral: 1 task (~110 reps) per sub-job
CHUNK=${CHUNK:-10}; CONC=${CONC:-6}          # low CONC: fvec stats are memory-bandwidth-bound;
# CONC=10 gave ~15x per-worker slowdown (net negative) on this 2-socket EPYC. Small CHUNK keeps the
# per-worker GenotypeArray small (less cache/bandwidth thrash). NUMA: pin workers across both sockets.
# node0 CPUs=0-63,128-191  node1 CPUs=64-127,192-255 (taskset; memory follows first-touch, no numactl).
SOCK0_CPUS=${SOCK0_CPUS:-0-63,128-191}; SOCK1_CPUS=${SOCK1_CPUS:-64-127,192-255}
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 TF_CPP_MIN_LOG_LEVEL=3
export PY H TMP CHUNK PARTS BANK SOCK0_CPUS SOCK1_CPUS NDIP
say "RETRAIN_START tag=$TAG conc=$CONC bank=$BANK"

# ---- manifest of sub-jobs: partname  position  seed  ms-files... ----
MAN=$FR/jobs.txt; : > "$MAN"; seed=1
emit_split(){  # position  "tasks"  tasks_per_job  dir  fname
  local pos=$1 tasks="$2" tpj=$3 dir=$4 fname=$5
  local arr=($tasks) i=0 p=0
  while [ $i -lt ${#arr[@]} ]; do
    local files=() j t
    for ((j=0; j<tpj && i<${#arr[@]}; j++, i++)); do t=${arr[$i]}; files+=("$H/$dir/task$t/$fname"); done
    echo "${pos}.part${p} $pos $((seed++)) ${files[*]}" >> "$MAN"; p=$((p+1))
  done
}
emit_split neutral "$NEUT_TASKS" "$NEUT_TASKS_PER_JOB" neutral neutral.msOut.gz
for cls in hard soft; do
  for i in 0 1 2 3 4 5 6 7 8 9 10; do
    if [ "$i" = 5 ]; then tks="$FOCAL_TASKS"; else tks="$LINKED_TASKS"; fi
    emit_split "${cls}_$i" "$tks" "$TASKS_PER_JOB" "$cls" "${cls}_$i.msOut.gz"
  done
done
say "manifest: $(wc -l < "$MAN") sub-jobs"

do_job(){
  local part=$1 pos=$2 seed=$3; shift 3; local ms=("$@")
  local out=$PARTS/$part.fvec
  [ -s "$out" ] && { echo "skip $part"; return; }
  # NUMA: alternate sockets by seed parity so concurrent workers split across both memory controllers
  local cpus; [ $((seed % 2)) -eq 0 ] && cpus="$SOCK0_CPUS" || cpus="$SOCK1_CPUS"
  # Cap focal/linked to REPS_FOCAL reps/file: some hard/soft ms files have ~160 reps (not the design's
  # 10) -> those jobs did 300+ reps = 16x wasted work (makeTrainingSets subsamples to 2000/class anyway).
  # Neutral is uncapped (wants all ~110 reps/file for the 2090 target).
  local npf=""; [ "$pos" != "neutral" ] && npf="--n-per-file ${REPS_FOCAL:-10}"
  local ndarg=""; [ -n "${NDIP:-}" ] && ndarg="--n-dip $NDIP"
  taskset -c "$cpus" "$PY" "$H/ms_to_vcf_fvec.py" --out "$out.tmp" --ms "${ms[@]}" --chunk "$CHUNK" --seed "$seed" \
      --tmp "$TMP" --geno-mask-bank "$BANK" $npf $ndarg > "$PARTS/$part.log" 2>&1 \
    && { mv "$out.tmp" "$out"; echo "done $part"; } || echo "FAIL $part"
}
export -f do_job; export PARTS REPS_FOCAL
say "FVEC_START ($(wc -l < "$MAN") jobs, conc=$CONC)"
xargs -P "$CONC" -L 1 bash -c 'do_job "$@"' _ < "$MAN"
say "PARTS_DONE ($(ls "$PARTS"/*.fvec 2>/dev/null | wc -l)/$(wc -l < "$MAN"))"

# ---- concat parts -> per-position fvec (header once) ----
for pos in neutral $(for c in hard soft; do for i in 0 1 2 3 4 5 6 7 8 9 10; do echo "${c}_$i"; done; done); do
  pf=("$PARTS/$pos.part"*.fvec)
  [ -e "${pf[0]}" ] || { say "CONCAT_MISS $pos"; continue; }
  { head -1 "${pf[0]}"; for f in "${pf[@]}"; do tail -n +2 "$f"; done; } > "$FV/$pos.fvec"
done
nf=$(ls "$FV"/*.fvec 2>/dev/null | wc -l)
say "FVEC_DONE ($nf/23; neutral=$(($(wc -l <"$FV/neutral.fvec" 2>/dev/null)-1)) hard_5=$(($(wc -l <"$FV/hard_5.fvec" 2>/dev/null)-1)))"
[ "$nf" -eq 23 ] || { say "FVEC_INCOMPLETE abort"; exit 1; }

# ---- strip 4 meta cols (fvecVcf) -> pure-numeric (fvecSim format train wants) ----
for f in "$FV"/*.fvec; do cut -f5- "$f" > "$FVN/$(basename "$f")"; done
export OMP_NUM_THREADS=8 OPENBLAS_NUM_THREADS=8 CUDA_VISIBLE_DEVICES=""
rm -rf "$TS"; mkdir -p "$TS"
say "MAKETS_START"
"$DE/diploSHIC" makeTrainingSets "$FVN/neutral.fvec" "$FVN/soft" "$FVN/hard" 5 0,1,2,3,4,6,7,8,9,10 "$TS/" >> "$FR/train.log" 2>&1 \
  && say "MAKETS_DONE $(for c in "$TS"/*.fvec; do printf '%s=%d ' "$(basename "$c" .fvec)" "$(($(wc -l <"$c")-1))"; done)" || { say "MAKETS_FAIL"; exit 1; }
say "TRAIN_START"
"$DE/diploSHIC" train "$TS/" "$TS/" "$FR/illexModel_vcf" --epochs 100 --numSubWins 11 --confusionFile "$FR/confusion_vcf.png" >> "$FR/train.log" 2>&1
if [ -s "$FR/illexModel_vcf.json" ]; then say "TRAIN_DONE acc=$(tr '\r' '\n' < "$FR/train.log" | grep -oE 'diploSHIC accuracy: [0-9.]+' | tail -1)"; say "RETRAIN_ALL_DONE"; else say "TRAIN_FAIL"; fi
