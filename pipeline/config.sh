#!/usr/bin/env bash
# Shared config for the illex popgen pipeline. Source at top of every run.sh.
set -euo pipefail
export POPGEN=/sietch_colab/data_share/illex/popgen_data
export ANA=$POPGEN/analysis
export REF=/sietch_colab/data_share/illex/andy_links/Illex_F24.primary.clean.fa
export GFF=$POPGEN/degenotate_illex/Illex_F24.gene_lnc_pseudo.func.fix.sq3.FINAL.v2.fixID.gff3
export ACCESS=$POPGEN/degenotate_illex/accessible_sites.bed
export DEGEN=$POPGEN/degenotate_illex/degeneracy-all-sites.bed
export GENES_BED=$POPGEN/degenotate_illex/Illex.genes.bed
export INTERGENIC_BED=$POPGEN/degenotate_illex/Illex.intergenic-access.bed
export INGROUP_TSV=$POPGEN/docs/samples_table.lp.tsv
export BAKER_VCFDIR=$POPGEN/seq_data/baker_2025/vcfs
export SMALL_VCFDIR=$POPGEN/seq_data/small_2026/vcfs
export CHR_ALL="Z 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40 41 42 43 44 45"
# Per-chromosome analyses use CHR_ALL (keep everything incl. chr Z).
# CHR_COMBINED = chrom set for genome-wide-COMBINED estimates only (global PCA /
# admixture / genome-wide beagle): drop the SEX chrom (chr42); KEEP chr Z (it is a
# real chromosome, not sex-linked). The chr2 inversion is excluded via region bed.
export CHR_COMBINED=$(echo $CHR_ALL | tr ' ' '\n' | grep -vx 42 | tr '\n' ' ')
# Back-compat alias (was mistakenly excluding Z); now == CHR_COMBINED.
export CHR_AUTO="$CHR_COMBINED"
export INV_CHR=2
export INV_START=60000000
export INV_END=80000000
export SEX_CHR=42
# conda envs (abspaths; -n is shadowed by grenepipe base)
export ENV_BCF=/home/ssmall/miniforge3/envs/mkado-vcf
export ENV_PGGPU=/home/ssmall/miniforge3/envs/pg_gpu
export ENV_PCANGSD=/home/ssmall/miniforge3/envs/pcangsd_lcpipe
export ENV_POPGEN=/home/ssmall/miniforge3/envs/popgen
BCF()   { mamba run --no-banner -p "$ENV_BCF" "$@"; }
export -f BCF 2>/dev/null || true   # bash-only; harmless no-op under zsh
