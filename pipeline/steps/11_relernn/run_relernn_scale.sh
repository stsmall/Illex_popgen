#!/usr/bin/env bash
# Scaled ReLERNN: per-sex autosomes (per-window landscape) + Z by sex (male diploid / female haploid).
# Uses demographic history (-n stairway). Drops zero-call-per-chrom individuals. Retry on GPU death.
set -uo pipefail
A=/sietch_colab/data_share/illex/popgen_data/analysis
BB=/home/ssmall/miniforge3/envs/bioinfo-buddy/bin
B=/home/ssmall/miniforge3/envs/grenepipe/envs/relernn_env/bin
ENV=/home/ssmall/miniforge3/envs/grenepipe/envs/relernn_env
RD=$A/steps/11_relernn; SD=$A/steps/12_sex
DEM=$RD/demHist.stairway.txt
export PYTHONNOUSERSITE=1 TF_USE_LEGACY_KERAS=1 XLA_FLAGS=--xla_gpu_cuda_data_dir=$ENV
AUT="1 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40 41 43 44 45"

# --- one full ReLERNN pipeline with retry (SIMULATE-n -> TRAIN -> PREDICT -> BSCORRECT) ---
run_pipe(){ local tag=$1 vcf=$2 gbed=$3 mask=$4 gpu=$5; local d=$RD/run_$tag; mkdir -p $d
  local L=$RD/logs; mkdir -p $L
  for attempt in 1 2 3; do
    rm -rf $d/proj
    echo "[$(date)] [$tag] attempt $attempt SIMULATE (gpu$gpu)"
    CUDA_VISIBLE_DEVICES=$gpu $B/ReLERNN_SIMULATE -v $vcf -g $gbed -m $mask -d $d/proj \
      -n $DEM -u 3e-9 -l 1 -t 12 --unphased --maskThresh 0.9 --nTrain 12000 --nVali 1000 --nTest 1000 > $L/$tag.sim.log 2>&1
    [ "$(ls $d/proj/train 2>/dev/null|wc -l)" -eq 0 ] && { echo "[$tag] SIM fail"; grep -iE 'Error|ZeroDiv|both haploid' $L/$tag.sim.log|tail -2; continue; }
    echo "[$(date)] [$tag] TRAIN"
    CUDA_VISIBLE_DEVICES=$gpu $B/ReLERNN_TRAIN -d $d/proj -t 12 --nEpochs 100 --nValSteps 20 > $L/$tag.train.log 2>&1 || { echo "[$tag] TRAIN died, retry"; continue; }
    echo "[$(date)] [$tag] PREDICT"
    CUDA_VISIBLE_DEVICES=$gpu $B/ReLERNN_PREDICT -v $vcf -d $d/proj --unphased > $L/$tag.pred.log 2>&1 || { echo "[$tag] PRED died, retry"; continue; }
    echo "[$(date)] [$tag] BSCORRECT"
    CUDA_VISIBLE_DEVICES=$gpu $B/ReLERNN_BSCORRECT -d $d/proj -t 12 --nSlice 20 --nReps 20 > $L/$tag.bs.log 2>&1 || { echo "[$tag] BS died, retry"; continue; }
    local out=$(find $d/proj -name '*BSCORRECTED.txt'|head -1)
    [ -n "$out" ] && { echo "[$(date)] [$tag] DONE -> $out"; awk 'NR==1{for(i=1;i<=NF;i++)if($i=="recombRate")r=i} NR>1{s+=$r;n++} END{if(n)printf "[%s] mean r=%.3e/bp (%.4f cM/Mb) n=%d\n","'$tag'",s/n,s/n*1e8,n}' "$out"; return 0; }
  done
  echo "[$(date)] [$tag] FAILED after 3 attempts"; return 1
}

# --- prep: drop zero-call-per-chrom union for a per-sex autosomal concat VCF ---
prep_auto(){ local sex=$1; local dir=$RD/$sex
  : > $dir/zerocall_union.txt
  for c in $AUT; do
    $BB/bcftools stats -s - $dir/parts/$c.vcf.gz 2>/dev/null | awk '/^PSC/{if($4+$5+$6==0)print $3}' >> $dir/zerocall_union.txt
  done
  sort -u $dir/zerocall_union.txt > $dir/drop.txt
  $BB/bcftools query -l $dir/$sex.vcf.gz | grep -vxF -f $dir/drop.txt > $dir/keep.txt 2>/dev/null || $BB/bcftools query -l $dir/$sex.vcf.gz > $dir/keep.txt
  $BB/bcftools view -S $dir/keep.txt $dir/$sex.vcf.gz -Oz -o $dir/$sex.kept.vcf.gz 2>/dev/null; $BB/bgzip -dk -f $dir/$sex.kept.vcf.gz
  echo "[$sex auto] dropped $(wc -l<$dir/drop.txt) zero-call; kept $(wc -l<$dir/keep.txt) indiv, $($BB/bcftools index -n $dir/$sex.kept.vcf.gz 2>/dev/null) SNPs"
}
# --- prep chrZ per-sex VCF (biallelic MAF>0.01 F_MISSING<0.5, drop zero-call) ---
prep_Z(){ local sex=$1; local dir=$RD/Z_$sex; mkdir -p $dir
  $BB/bcftools view -S $SD/${sex}s.txt --force-samples -m2 -M2 -v snps $A/steps/00_callset/filtered/Z/variants_filt.vcf.gz -Ou 2>/dev/null \
   | $BB/bcftools view -e 'F_MISSING>0.5' -Ou 2>/dev/null | $BB/bcftools +fill-tags -Ou -- -t MAF 2>/dev/null \
   | $BB/bcftools view -e 'INFO/MAF<0.01' -Oz -o $dir/Z_$sex.vcf.gz 2>/dev/null
  $BB/bcftools stats -s - $dir/Z_$sex.vcf.gz 2>/dev/null | awk '/^PSC/{if($4+$5+$6==0)print $3}' > $dir/drop.txt
  $BB/bcftools query -l $dir/Z_$sex.vcf.gz | grep -vxF -f $dir/drop.txt > $dir/keep.txt 2>/dev/null || $BB/bcftools query -l $dir/Z_$sex.vcf.gz > $dir/keep.txt
  $BB/bcftools view -S $dir/keep.txt $dir/Z_$sex.vcf.gz -Oz -o $dir/Z_$sex.kept.vcf.gz 2>/dev/null; $BB/bgzip -dk -f $dir/Z_$sex.kept.vcf.gz
  local Lz=$(awk '$1=="Z"{print $2}' /sietch_colab/data_share/illex/andy_links/Illex_F24.primary.clean.fa.fai)
  printf 'Z\t0\t%s\n' "$Lz" > $dir/genome.bed; awk '$1=="Z"' $RD/inaccessible.bed > $dir/inacc.bed
  echo "[Z $sex] dropped $(wc -l<$dir/drop.txt); kept $(wc -l<$dir/keep.txt), $($BB/bcftools index -n $dir/Z_$sex.kept.vcf.gz 2>/dev/null) SNPs"
}

# wait for autosomal per-sex concat VCFs (from build_persex_vcf.sh)
while [ ! -s $RD/male/male.vcf.gz ] || [ ! -s $RD/female/female.vcf.gz ]; do sleep 30; done
echo "[$(date)] per-sex autosomal VCFs ready; prepping"
prep_auto male; prep_auto female; prep_Z male; prep_Z female

# launch 4 pipelines: autosomes on gpu0/1 (parallel), Z on gpu2 (sequential)
run_pipe male_auto   $RD/male/male.kept.vcf     $RD/genome.bed        $RD/inaccessible.bed 0 &
run_pipe female_auto $RD/female/female.kept.vcf $RD/genome.bed        $RD/inaccessible.bed 1 &
( run_pipe Z_male $RD/Z_male/Z_male.kept.vcf $RD/Z_male/genome.bed $RD/Z_male/inacc.bed 2
  run_pipe Z_female $RD/Z_female/Z_female.kept.vcf $RD/Z_female/genome.bed $RD/Z_female/inacc.bed 2 ) &
wait
echo "[$(date)] ALL ReLERNN SCALE RUNS DONE"
