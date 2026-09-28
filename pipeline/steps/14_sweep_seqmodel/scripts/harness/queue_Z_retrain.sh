#!/usr/bin/env bash
# ORCHESTRATOR: chrZ-males sweep-scan at n=N diploid, done right (self-contained).
# All 633 baker samples are illecebrosus; the "350 illex set" was only a subsample sized to the
# n=350 sims. The real male pool is 330 confident ZZ males (12_sex/males.txt) -> use ALL of them.
# diploSHIC features are n-dependent, so chrZ needs its OWN model trained at n=N.
# Chains, unattended:
#   1. build the N-male chrZ VCF subset (all.Z.vcf.gz -> SUBLIST, biallelic SNPs, GT only)
#   2. build an N-sample chrZ geno-mask bank (real correlated missingness for the overlay)
#   3. RETRAIN at NDIP=N with that bank  -> results/<TAG>/illexModel_vcf
#   4. build the N-male empirical chrZ fvec + PREDICT with the Z model
#      -> results/<OUTDIR>/preds/chrZ.preds  (SEPARATE from the autosomal genome.preds:
#         different model/n/sex chromosome -- do NOT merge).
# Config (env, defaults = the 330-male run):
#   N=330  SUBLIST=<12_sex/males.txt>  CONC=4  CHUNK=10
# Shared box: CONC low, OMP=1, taskset-pinned (from run_fvecvcf_retrain.sh).
# Detached:  N=330 setsid bash queue_Z_retrain.sh >/dev/null 2>&1 </dev/null &
set -uo pipefail
D=/sietch_colab/data_share/illex/popgen_data/analysis/steps/14_sweep_seqmodel
B=/sietch_colab/data_share/illex/popgen_data
H=$D/scripts/harness
PY=/home/ssmall/miniforge3/envs/slim_sim/bin/python
DE=/home/ssmall/miniforge3/envs/diploshic_env/bin
BCF=/home/ssmall/bin/bcftools
N=${N:-330}
SUBLIST=${SUBLIST:-$B/analysis/steps/12_sex/males.txt}
TAG=${TAG:-vcfretrain_Z${N}}
OUTDIR=${OUTDIR:-empirical_scan_Z${N}}
ZVCF=$D/results/task17_maskval/Zmales${N}.Z.vcf.gz
MASKZ=$D/results/empirical_scan/masks/mask.Z.fa
BANK=$D/results/task17_maskval/geno_mask_bank.chrZ.${N}.npz
S2Z=$D/results/task17_maskval/s2p_Z${N}.tsv     # sample<TAB>illex for the N males (for fvecVcf targetPop)
MODEL=$D/results/$TAG/illexModel_vcf
ZOUT=$D/results/$OUTDIR; mkdir -p "$ZOUT"/{fvec,preds,logs}
LOG=$D/results/$TAG.queue.log; mkdir -p "$D/results/$TAG"; : > "$LOG"
say(){ echo "$(date -u +%FT%TZ) $*" >> "$LOG"; }
FAI=/sietch_colab/data_share/illex/assembly/current_assembly_version/Illex_F24.primary.clean.fa.fai
say "QUEUE_START tag=$TAG N=$N sublist=$SUBLIST"

# sample list: SUBLIST samples that are actually in the Z callset
"$BCF" query -l "$B/seq_data/baker_2025/vcfs/all.Z.vcf.gz" | sort > "$ZOUT/callset_samples.txt"
sort "$SUBLIST" | comm -12 - "$ZOUT/callset_samples.txt" > "$ZOUT/Zmales.list"
nreq=$(wc -l < "$ZOUT/Zmales.list")
say "SAMPLE_LIST n=$nreq (of N=$N requested)"
[ "$nreq" -ge 1 ] || { say "NO_SAMPLES abort"; exit 1; }
awk '{print $1"\tillex"}' "$ZOUT/Zmales.list" > "$S2Z"

# ---- 1. build the N-male chrZ VCF ----
if [ ! -f "$ZVCF.done" ]; then
  say "BUILD_ZVCF"
  "$BCF" view -S "$ZOUT/Zmales.list" "$B/seq_data/baker_2025/vcfs/all.Z.vcf.gz" -Ou 2>>"$ZOUT/logs/subset.log" \
    | "$BCF" view -m2 -M2 -v snps -Ou 2>>"$ZOUT/logs/subset.log" \
    | "$BCF" annotate -x '^FORMAT/GT,INFO' -Oz -o "$ZVCF" 2>>"$ZOUT/logs/subset.log" \
    && "$BCF" index -t "$ZVCF" && echo done > "$ZVCF.done" || { say "ZVCF_FAIL abort"; exit 1; }
fi
nz=$("$BCF" query -l "$ZVCF" | wc -l); say "ZVCF_READY samples=$nz sites=$("$BCF" index -n "$ZVCF")"

# ---- 2. build the N-sample chrZ geno-mask bank ----
if [ ! -s "$BANK" ]; then
  say "BUILD_BANK n_indiv=$nz"
  "$PY" "$H/structured_mask.py" --vcf "$ZVCF" --fasta "$MASKZ" --chrom Z --out "$BANK" --n-indiv "$nz" >> "$LOG" 2>&1 \
    || { say "BANK_FAIL abort"; exit 1; }
fi
say "BANK_READY size=$(ls -la "$BANK" | awk '{print $5}')"

# ---- 3. retrain at NDIP=nz with the Z bank ----
say "RETRAIN_LAUNCH (NDIP=$nz BANK=chrZ.$N CONC=${CONC:-4})"
NDIP="$nz" BANK="$BANK" CONC="${CONC:-4}" CHUNK="${CHUNK:-10}" bash "$H/run_fvecvcf_retrain.sh" "$TAG" >> "$LOG" 2>&1
[ -s "$MODEL.json" ] || { say "RETRAIN_FAIL (no $MODEL.json) abort"; exit 1; }
say "RETRAIN_DONE model=$MODEL"

# ---- 4. chrZ-male empirical fvec (n=nz) + predict with the Z model ----
c=Z; L=$(awk '$1=="Z"{print $2}' "$FAI")
if [ ! -s "$ZOUT/fvec/chr$c.fvec" ]; then
  say "ZSCAN_FVEC"
  "$DE/diploSHIC" fvecVcf diploid "$ZVCF" "$c" "$L" "$ZOUT/fvec/chr$c.fvec" \
    --targetPop illex --sampleToPopFileName "$S2Z" \
    --winSize 1100000 --numSubWins 11 --maskFileName "$MASKZ" \
    --unmaskedFracCutoff 0.25 --unmaskedGenoFracCutoff 0.5 > "$ZOUT/logs/chr$c.fvec.log" 2>&1
fi
[ -s "$ZOUT/fvec/chr$c.fvec" ] || { say "ZSCAN_FVEC_FAIL abort"; exit 1; }
say "ZSCAN_FVEC_DONE ($(($(wc -l <"$ZOUT/fvec/chr$c.fvec")-1)) win)"
export OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2 CUDA_VISIBLE_DEVICES="" TF_CPP_MIN_LOG_LEVEL=3
say "ZSCAN_PREDICT"
"$DE/diploSHIC" predict "$MODEL.json" "$MODEL.weights.h5" "$ZOUT/fvec/chr$c.fvec" "$ZOUT/preds/chr$c.preds" --numSubWins 11 > "$ZOUT/logs/chr$c.predict.log" 2>&1
if [ -s "$ZOUT/preds/chr$c.preds" ]; then
  say "ZSCAN_DONE preds=$(($(wc -l <"$ZOUT/preds/chr$c.preds")-1)) win classDist=[$(tail -n +2 "$ZOUT/preds/chr$c.preds"|cut -f5|sort|uniq -c|tr '\n' ';')]"
else say "ZSCAN_PREDICT_FAIL"; fi
say "QUEUE_DONE"
