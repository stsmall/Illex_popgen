#!/usr/bin/env bash
# Permutation null for per-SNP geographic Hudson FST (MAF > 0.05): division labels shuffled
# among the 632 individuals NPERM times (perm_00 = real labels). Every permutation keeps the
# real genotypes, genotype-calling behaviour, per-SNP missingness, frequency spectrum and
# linkage, and removes only geography. Output per chrom: perm/NN/<chrom>.fst.npy aligned SNP
# by SNP, plus <chrom>.pos.npy and <chrom>.maf.npy. 4 chroms x 8 threads = 32 cores.
set -uo pipefail
HERE=$(cd "$(dirname "$0")" && pwd)
OUT=${OUT:-/sietch_colab/data_share/illex/popgen_data/analysis/steps/17_ensemble_scan/fst_maf05/perm}
VCFDIR=${VCFDIR:-/sietch_colab/data_share/illex/popgen_data/pggpu_illex/popstats/snps_stitch}
PLINK=${PLINK:-/home/ssmall/bin/plink2}
PY=${PY:-/home/ssmall/miniforge3/envs/mkado-vcf/bin/python}
NPERM=${NPERM:-20}
mkdir -p "$OUT/labels" "$OUT/work" "$OUT/logs"
"$PY" "$HERE/make_permutations.py" "$OUT/labels" "$NPERM"
awk 'NR>1{print $1}' "$OUT/labels/perm_00.tsv" > "$OUT/keep_632.txt"
export OUT VCFDIR PLINK PY NPERM HERE
do_one(){
  local c=$1 w=$OUT/work/chr$1
  [ -f "$OUT/done_$c" ] && { echo "SKIP $c"; return; }
  "$PLINK" --vcf "$VCFDIR/all.$c.snps.vcf.gz" --chr-set 45 --max-alleles 2 --set-all-var-ids @:# \
      --keep "$OUT/keep_632.txt" --maf 0.05 --make-pgen --threads 8 --out "$w" > /dev/null 2>&1 || { echo "PGEN_FAIL $c"; return; }
  for k in $(seq -f "%02g" 0 "$NPERM"); do
    mkdir -p "$OUT/$k"
    [ -s "$OUT/$k/$c.fst.npy" ] && continue
    "$PLINK" --pfile "$w" --pheno "$OUT/labels/perm_$k.tsv" --loop-cats POP --freq counts --threads 8 \
        --out "$w.p$k" > /dev/null 2>&1 || { echo "FREQ_FAIL $c $k"; return; }
    "$PY" "$HERE/hudson_maf05_agg.py" "$c" "$w.p$k" "$OUT/$k" npy || { echo "AGG_FAIL $c $k"; return; }
    rm -f "$w.p$k".D*.acount "$w.p$k".log
  done
  rm -f "$w".pgen "$w".pvar "$w".psam "$w".log
  touch "$OUT/done_$c"; echo "DONE $c"
}
export -f do_one
printf '%s\n' $(seq 1 45) Z | xargs -P 4 -I{} bash -c 'do_one "$@"' _ {} 2>&1 | grep -v Warning >> "$OUT/logs/run_perm_null.log"
touch "$OUT/ALL_DONE"
