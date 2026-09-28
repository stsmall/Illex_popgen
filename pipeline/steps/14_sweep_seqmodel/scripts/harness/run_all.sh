#!/usr/bin/env bash
# Box diploSHIC chain, RAM-safe. fvec is 53s/rep + a fixed ~183s load/invocation + ~4.9GB
# RAM/proc, so: group raw task-files into ~46 chunks of ~250-330 blocks (amortizes the load
# to ~1%, keeps -P CONC fed), combine_ms.sh each (accurate header -> no EOF), fvec, merge to
# 23 per-class-position fvecs, makeTrainingSets (balances all 5 classes to 2000), train.
# Run detached:  setsid bash run_all.sh >/dev/null 2>&1 </dev/null &
set -uo pipefail
D=/sietch_colab/data_share/illex/popgen_data/analysis/steps/14_sweep_seqmodel
H=$D/scripts/harness
FR=$D/results/fullrun; CH=$FR/chunkms; RAWFV=$FR/rawfvec; FV=$FR/fvec; LOGD=$FR/fveclogs
mkdir -p "$CH" "$RAWFV" "$FV" "$LOGD"
DE=/home/ssmall/miniforge3/envs/diploshic_env/bin
COMBINE=$H/combine_ms.sh; NHAP=700; CONC=14
VCF=$D/results/task17_maskval/cleaned.1.vcf.gz
S2P=$D/results/task17_maskval/s2p.tsv
MASK=$H/talapas/mask/mask.fa
PROG=$FR/run_all.progress
export COMPRESS="pigz -1 -p 2"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
export H CH RAWFV LOGD DE COMBINE NHAP VCF S2P MASK COMPRESS
say(){ echo "$(date -u +%FT%TZ) $*" >> "$PROG"; }
say "RUNALL_START pid=$$ conc=$CONC"

# ---- chunk manifest: name  cls  tfirst  tlast  chunkid ----
MAN=$FR/chunk_manifest.txt; : > "$MAN"
addchunks(){ local name=$1 cls=$2 nt=$3 per=$4 c=0 t=0 last
  while [ $t -lt "$nt" ]; do last=$((t+per-1)); [ $last -ge "$nt" ] && last=$((nt-1))
    echo "$name $cls $t $last $c" >> "$MAN"; c=$((c+1)); t=$((last+1)); done; }
addchunks neutral neutral 28 3            # 110 blk/task -> ~330/chunk, 10 chunks
addchunks hard_5 hard 200 25              # 10 blk/task  -> 250/chunk, 8 chunks (focal, full 2000)
addchunks soft_5 soft 200 25
for i in 0 1 2 3 4 6 7 8 9 10; do         # linked: 30 tasks = 1 chunk of 300 blk each
  addchunks "hard_$i" hard 30 30; addchunks "soft_$i" soft 30 30; done
say "chunks: $(wc -l < "$MAN")"

# ---- worker: combine a task group -> chunk ms -> fvec -> chunk fvec ----
do_chunk(){ local name=$1 cls=$2 tf=$3 tl=$4 cid=$5
  local cms=$CH/${name}.c${cid}.msOut.gz out=$RAWFV/${name}.c${cid}.fvec lock=$RAWFV/${name}.c${cid}.lock
  [ -s "$out" ] && { echo "skip $name.c$cid"; return; }
  mkdir "$lock" 2>/dev/null || return
  local files=() t; for ((t=tf; t<=tl; t++)); do files+=("$H/$cls/task$t/${name}.msOut.gz"); done
  bash "$COMBINE" "$cms" "$NHAP" "${files[@]}" >/dev/null 2>&1
  if "$DE/diploSHIC" fvecSim diploid "$cms" "$out.tmp" --totalPhysLen 1100000 --numSubWins 11 \
       --maskFileName "$MASK" --chrArmsForMasking 1 --vcfForMaskFileName "$VCF" \
       --popForMask illex --sampleToPopFileName "$S2P" \
       --unmaskedGenoFracCutoff 0.5 --unmaskedFracCutoff 0.25 > "$LOGD/${name}.c${cid}.log" 2>&1; then
    mv "$out.tmp" "$out"; rm -f "$cms"; echo "done $name.c$cid ($(($(wc -l <"$out")-1)) rows)"
  else rm -f "$out.tmp" "$cms"; echo "FAIL $name.c$cid"; fi
  rmdir "$lock" 2>/dev/null; }
export -f do_chunk

say "FVEC_START"
rm -rf "$RAWFV"/*.lock 2>/dev/null
xargs -P "$CONC" -L 1 bash -c 'do_chunk "$@"' _ < "$MAN"
rm -rf "$RAWFV"/*.lock 2>/dev/null
say "FVEC_CHUNKS_DONE ($(ls "$RAWFV"/*.fvec 2>/dev/null | wc -l) chunk fvecs)"

# ---- merge chunk fvecs -> per class-position fvec (header once + all data rows) ----
names="neutral hard_5 soft_5"; for i in 0 1 2 3 4 6 7 8 9 10; do names="$names hard_$i soft_$i"; done
for name in $names; do
  cf=("$RAWFV/${name}.c"*.fvec)
  [ -e "${cf[0]}" ] || { say "MERGE_MISS $name"; continue; }
  { head -1 "${cf[0]}"; for f in "${cf[@]}"; do tail -n +2 "$f"; done; } > "$FV/$name.fvec"
done
nf=$(ls "$FV"/*.fvec 2>/dev/null | wc -l)
say "FVEC_DONE ($nf/23; neutral=$(($(wc -l <"$FV/neutral.fvec")-1)) hard_5=$(($(wc -l <"$FV/hard_5.fvec")-1)) hard_0=$(($(wc -l <"$FV/hard_0.fvec")-1)))"
[ "$nf" -eq 23 ] || { say "FVEC_INCOMPLETE abort"; exit 1; }

# ---- makeTrainingSets + train (CPU-pinned to spare the GPUs) ----
export OMP_NUM_THREADS=8 OPENBLAS_NUM_THREADS=8 TF_CPP_MIN_LOG_LEVEL=3 CUDA_VISIBLE_DEVICES=""
TS=$FR/trainingSets; rm -rf "$TS"; mkdir -p "$TS"
say "MAKETS_START"
"$DE/diploSHIC" makeTrainingSets "$FV/neutral.fvec" "$FV/soft" "$FV/hard" 5 0,1,2,3,4,6,7,8,9,10 "$TS/" >> "$FR/train.log" 2>&1
say "MAKETS_DONE $(for c in "$TS"/*.fvec; do printf '%s=%d ' "$(basename "$c" .fvec)" "$(($(wc -l <"$c")-1))"; done)"
say "TRAIN_START"
"$DE/diploSHIC" train "$TS/" "$TS/" "$FR/illexModel" --epochs 100 --numSubWins 11 --confusionFile "$FR/confusion.png" >> "$FR/train.log" 2>&1
if [ -s "$FR/illexModel.json" ] && [ -s "$FR/illexModel.weights.h5" ]; then say "TRAIN_DONE model=$FR/illexModel.json"; say "RUNALL_DONE"
else say "TRAIN_FAIL (see $FR/train.log)"; fi
