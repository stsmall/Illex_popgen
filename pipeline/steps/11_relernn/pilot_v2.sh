#!/usr/bin/env bash
# Pilot v2 chr45 male: drop zero-call individuals (fast via bcftools stats PSC), then ReLERNN
# SIMULATE(no --forceDiploid) -> TRAIN -> PREDICT -> BSCORRECT. Sequential/blocking (no wait-loop bug).
set -uo pipefail
A=/sietch_colab/data_share/illex/popgen_data/analysis
BB=/home/ssmall/miniforge3/envs/bioinfo-buddy/bin
B=/home/ssmall/miniforge3/envs/grenepipe/envs/relernn_env/bin
P=$A/steps/11_relernn/pilot_chr45_male
export PYTHONNOUSERSITE=1 CUDA_VISIBLE_DEVICES=0
cd $P
# 1) per-sample non-missing counts -> keep samples with >=1 called genotype (PSC: nRefHom+nNonRefHom+nHets)
$BB/bcftools stats -s - chr45_male.vcf.gz 2>/dev/null | awk '/^PSC/{nm=$4+$5+$6; if(nm>0) print $3}' > keep.txt
echo "[$(date)] keep $(wc -l < keep.txt) / 330 males (dropped zero-call)"
$BB/bcftools view -S keep.txt chr45_male.vcf.gz -Oz -o chr45_male.kept.vcf.gz 2>/dev/null
$BB/bgzip -dk -f chr45_male.kept.vcf.gz
# 2) SIMULATE (blocking)
rm -rf proj; echo "[$(date)] SIMULATE"
$B/ReLERNN_SIMULATE -v chr45_male.kept.vcf -g genome.bed -m inacc.bed -d proj \
  -u 3e-9 -l 1 -t 12 --unphased --maskThresh 0.9 --nTrain 2000 --nVali 200 --nTest 200 > sim2.log 2>&1
ns=$(ls proj/train 2>/dev/null | wc -l)
echo "[$(date)] SIMULATE done: train files=$ns"
if [ "$ns" -eq 0 ]; then echo "SIMULATE FAILED"; grep -iE 'Error|ZeroDiv|thetaW|nSamps|both haploid' sim2.log | grep -vi tensorflow | tail -3; exit 1; fi
# 3) TRAIN/PREDICT/BSCORRECT
echo "[$(date)] TRAIN"; $B/ReLERNN_TRAIN -d proj -t 12 --nEpochs 50 --nValSteps 10 > train2.log 2>&1
echo "[$(date)] PREDICT"; $B/ReLERNN_PREDICT -v chr45_male.kept.vcf -d proj -u 3e-9 -l 1 > predict2.log 2>&1
echo "[$(date)] BSCORRECT"; $B/ReLERNN_BSCORRECT -d proj -t 12 --nSlice 20 --nReps 20 > bscorrect2.log 2>&1
echo "[$(date)] PILOT v2 DONE"; find proj -name '*PREDICT*' 2>/dev/null
