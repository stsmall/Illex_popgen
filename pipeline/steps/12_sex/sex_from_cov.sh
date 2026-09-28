#!/usr/bin/env bash
# Assign sex from sex-chromosome coverage: per-BAM samtools idxstats -> normalized read density
# on chrZ and chr42 relative to autosomes (chr1,3-41,43-45; excl Z,42). Bimodal ratio -> sex.
set -uo pipefail
A=/sietch_colab/data_share/illex/popgen_data/analysis
BB=/home/ssmall/miniforge3/envs/bioinfo-buddy/bin
SD=$A/steps/12_sex; mkdir -p $SD
BL=$A/steps/04_angsd_chr2/bamlists/bakerMeta_bamlist.txt
export BB
one(){ local bam=$1
  local s=$(basename "$bam" | sed 's/\.[^.]*$//')
  $BB/samtools idxstats "$bam" 2>/dev/null | awk -v S="$s" '
    $1=="Z"{zr=$3; zl=$2; next}
    $1=="42"{wr=$3; wl=$2; next}
    $1 ~ /^([0-9]+)$/ && $1!=42 {ar+=$3; al+=$2}
    END{
      if(al>0 && zl>0 && wl>0){
        ad=ar/al; zd=zr/zl; wd=wr/wl;
        printf "%s\t%.4f\t%.4f\t%.0f\n", S, (ad>0?zd/ad:0), (ad>0?wd/ad:0), ar
      }
    }'
}
export -f one
echo -e "sample\tZ_ratio\tchr42_ratio\tautosome_reads" > $SD/sex_cov.tsv
cat $BL | xargs -P24 -I{} bash -c 'one "$@"' _ {} >> $SD/sex_cov.tsv
n=$(($(wc -l < $SD/sex_cov.tsv)-1))
echo "computed ratios for $n individuals -> sex_cov.tsv"
echo "=== Z_ratio distribution (histogram-ish) ==="
tail -n +2 $SD/sex_cov.tsv | awk '{printf "%.1f\n",$2}' | sort -n | uniq -c
echo "=== chr42_ratio distribution ==="
tail -n +2 $SD/sex_cov.tsv | awk '{printf "%.1f\n",$3}' | sort -n | uniq -c
