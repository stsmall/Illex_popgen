#!/usr/bin/env bash
# Run mkado vcf (MK test) for illecebrosus ingroup vs each outgroup.
# Outputs: per-gene TSV + volcano, and asymptotic aggregate + alpha(x) plot.
set -euo pipefail

ROOT=/sietch_colab/data_share/illex/popgen_data/mkado_illex
IN=$ROOT/inputs
REF=/sietch_colab/data_share/illex/andy_links/Illex_F24.primary.clean.fa
GFF=$ROOT/Illex_F24.gene_lnc_pseudo.func.fix.sq3.FINAL.v2.fixID.gff3
VCF=$IN/baker_illecebrosus.cds.vcf.gz
ENV=/home/ssmall/miniforge3/envs/mkado-vcf
WORKERS=16

MK() { mamba run --no-banner -p "$ENV" mkado "$@"; }

run_outgroup() {
  local name=$1 ogvcf=$2
  local out=$ROOT/results/$name
  mkdir -p "$out"
  echo "[$(date)] === $name : per-gene + volcano ==="
  MK vcf --vcf "$VCF" --ref "$REF" --gff "$GFF" --outgroup-vcf "$ogvcf" \
      --per-gene -f tsv -w $WORKERS \
      --volcano "$out/$name.volcano.pdf" \
      --output "$out/$name.per_gene.tsv"

  echo "[$(date)] === $name : asymptotic aggregate ==="
  MK vcf --vcf "$VCF" --ref "$REF" --gff "$GFF" --outgroup-vcf "$ogvcf" \
      -a -f tsv -w $WORKERS \
      --plot-asymptotic "$out/$name.asymptotic.pdf" \
      --output "$out/$name.asymptotic.tsv"
  echo "[$(date)] === $name DONE ==="
}

echo "[$(date)] mkado $(MK --version 2>&1)"
run_outgroup coindetii  "$IN/coindetii.vcf.gz"
run_outgroup argentinus "$IN/argentinus.cds.vcf.gz"
echo "[$(date)] === ALL MKADO RUNS DONE ==="
