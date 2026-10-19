#!/usr/bin/env bash
# Genome-wide per-SNP geographic FST on EVERY biallelic SNP with pooled MAF > 0.05.
# Hudson estimator per division pair, ratio of averages per SNP (Bhatia et al. 2013).
# Per chrom: plink2 MAF filter on the 632 individuals with a division (3N, n=1, excluded)
# -> per-division allele counts (--loop-cats) -> hudson_maf05_agg.py; also writes
# per-sample missingness for the neutral null. 4 chroms x 8 threads = 32 cores.
set -uo pipefail
HERE=$(cd "$(dirname "$0")" && pwd)
OUT=${OUT:-/sietch_colab/data_share/illex/popgen_data/analysis/steps/17_ensemble_scan/fst_maf05}
VCFDIR=${VCFDIR:-/sietch_colab/data_share/illex/popgen_data/pggpu_illex/popstats/snps_stitch}
PLINK=${PLINK:-/home/ssmall/bin/plink2}
PY=${PY:-/home/ssmall/miniforge3/envs/mkado-vcf/bin/python}
DIV=$HERE/../config/divisions_632.tsv
mkdir -p "$OUT/work" "$OUT/perchrom" "$OUT/logs"
awk 'NR>1{print $1}' "$DIV" > "$OUT/keep_632.txt"
export OUT VCFDIR PLINK PY DIV HERE
do_one(){
  local c=$1 w=$OUT/work/chr$1
  [ -s "$OUT/perchrom/hudson_$c.tsv" ] && { echo "SKIP $c"; return; }
  "$PLINK" --vcf "$VCFDIR/all.$c.snps.vcf.gz" --chr-set 45 --max-alleles 2 --set-all-var-ids @:# \
      --keep "$OUT/keep_632.txt" --maf 0.05 --make-pgen --threads 8 --out "$w" > /dev/null 2>&1 || { echo "PGEN_FAIL $c"; return; }
  "$PLINK" --pfile "$w" --pheno "$DIV" --loop-cats POP --freq counts --threads 8 --out "$w" > /dev/null 2>&1 || { echo "FREQ_FAIL $c"; return; }
  "$PLINK" --pfile "$w" --missing sample-only --threads 8 --out "$w" > /dev/null 2>&1
  "$PY" "$HERE/hudson_maf05_agg.py" "$c" "$w" "$OUT/perchrom" || { echo "AGG_FAIL $c"; return; }
  /bin/cp -f "$w.smiss" "$OUT/perchrom/miss_$c.smiss"
  rm -f "$w".pgen "$w".pvar "$w".psam "$w".D*.acount
  echo "DONE $c"
}
export -f do_one
printf '%s\n' $(seq 1 45) Z | xargs -P 4 -I{} bash -c 'do_one "$@"' _ {} 2>&1 | tee -a "$OUT/logs/run_gw_maf05.log"
touch "$OUT/ALL_DONE"
