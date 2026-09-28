#!/usr/bin/env bash
# Genome-wide per-SNP geographic FST scan on the FULL per-chrom callsets.
# Replicates the chr2 full-callset method: plink2 --fst POP method=wc report-variants,
# categorical D-prefixed division pheno, --thin 0.02 seed 20240915, --max-alleles 2.
# Per chrom: run plink2 (45 pair .fst.var) -> aggregate to per-SNP mean -> delete raw.
# Concurrency 4 x 8 threads = 32 cores (<=40 budget on shared box).
set -uo pipefail
SP=/tmp/claude-1003/-sietch-colab-data-share-illex-popgen-data-mkado-illex/6d7afefb-6654-47e0-96f2-bc56ef699fbb/scratchpad
GW=$SP/gwfst
VCFDIR=/sietch_colab/data_share/illex/popgen_data/pggpu_illex/popstats/snps_stitch
PLINK=/home/ssmall/bin/plink2
PY=/home/ssmall/miniforge3/envs/mkado-vcf/bin/python
POP=$GW/pop2.txt
PROG=$GW/progress.log
: > "$PROG"
say(){ echo "$(date -u +%FT%TZ) $*" >> "$PROG"; }
export GW VCFDIR PLINK PY POP PROG
say "GW_FST_START"

do_one(){
  local c=$1
  local wd=$GW/chr$c
  mkdir -p "$wd"
  if [ -s "$GW/perchrom_$c.tsv.done" ]; then echo "$(date -u +%FT%TZ) SKIP $c" >>"$PROG"; return; fi
  echo "$(date -u +%FT%TZ) PLINK_START $c" >> "$PROG"
  "$PLINK" --vcf "$VCFDIR/all.$c.snps.vcf.gz" \
      --chr-set 45 --max-alleles 2 --set-all-var-ids @:# \
      --thin 0.02 --seed 20240915 \
      --pheno "$POP" \
      --fst POP method=wc report-variants \
      --threads 8 --out "$wd/fst_$c" >> "$wd/plink.log" 2>&1
  local rc=$?
  if [ $rc -ne 0 ]; then echo "$(date -u +%FT%TZ) PLINK_FAIL $c rc=$rc" >> "$PROG"; return; fi
  echo "$(date -u +%FT%TZ) AGG $c" >> "$PROG"
  "$PY" "$GW/agg_chrom.py" "$c" "$wd" >> "$wd/agg.log" 2>&1
  if [ -s "$wd/perchrom_$c.tsv" ]; then
    cp "$wd/perchrom_$c.tsv" "$GW/perchrom_$c.tsv"
    touch "$GW/perchrom_$c.tsv.done"
    rm -f "$wd"/fst_$c.*.fst.var "$wd"/fst_$c-temporary.*      # free disk
    local n=$(($(wc -l < "$GW/perchrom_$c.tsv")-1))
    echo "$(date -u +%FT%TZ) CHROM_DONE $c nsnp=$n" >> "$PROG"
  else
    echo "$(date -u +%FT%TZ) AGG_FAIL $c" >> "$PROG"
  fi
}
export -f do_one

printf '%s\n' $(seq 1 45) Z | xargs -P 4 -I{} bash -c 'do_one "$@"' _ {}
say "GW_FST_DONE"
touch "$GW/ALL_DONE"
