#!/usr/bin/env bash
# Callability-filtered re-rank of the SF2 (CLR) + RAiSD (mu) scans: annotate each grid
# window with accessible fraction + local SNP density, drop uncallable/sparse windows
# (the chr37-type artifacts), re-rank the survivors.
set -uo pipefail
B=/sietch_colab/data_share/illex/popgen_data
R=$B/analysis/steps/17_ensemble_scan/results
ACC=$B/degenotate_illex/accessible_sites.bed
FAI=/sietch_colab/data_share/illex/assembly/current_assembly_version/Illex_F24.primary.clean.fa.fai
BT=/home/ssmall/miniforge3/envs/annotation-expression-buddy-barrnap_trnascan/bin/bedtools
CH="1 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40 41 43 44 45"
ACC_MIN=0.5          # window must be >=50% accessible
DENS_MIN=5.0         # window must have >=5 SNPs/kb
OUT=$R/filtered; mkdir -p "$OUT"
export LC_ALL=C

awk 'BEGIN{OFS="\t"}{print $1,$2}' "$FAI" | sort -k1,1 > "$OUT/genome.txt"
sort -k1,1 -k2,2n "$ACC" > "$OUT/acc.sorted.bed"

# genome SNP bed from the per-chrom FreqFiles
for c in $CH; do awk -v c="$c" 'BEGIN{OFS="\t"} NR>1{print c,($1-1),$1}' "$R/freq/$c.sf2.freq"; done \
  | sort -k1,1 -k2,2n > "$OUT/snps.sorted.bed"
echo "SNP bed: $(wc -l < "$OUT/snps.sorted.bed") SNPs"

rerank() {   # method  halfwin
  m=$1; w=$2
  for c in $CH; do
    L=$(awk -v c="$c" '$1==c{print $2}' "$FAI")
    if [ "$m" = sf2 ]; then
      tail -n +2 "$R/sf2/$c.clr" | awk -v c="$c" -v L="$L" -v w="$w" 'BEGIN{OFS="\t"}{p=int($1);s=p-w;if(s<0)s=0;e=p+w;if(e>L)e=L;print c,s,e,$2,p}'
    else
      grep -v '^//' "$R/raisd/RAiSD_Report.$c" | awk -v c="$c" -v L="$L" -v w="$w" 'BEGIN{OFS="\t"} NF==7{p=$1;s=p-w;if(s<0)s=0;e=p+w;if(e>L)e=L;print c,s,e,$7,p}'
    fi
  done | sort -k1,1 -k2,2n > "$OUT/${m}_win.bed"
  "$BT" coverage -a "$OUT/${m}_win.bed" -b "$OUT/acc.sorted.bed" -sorted -g "$OUT/genome.txt" > "$OUT/${m}_acc.txt" 2>/dev/null
  "$BT" coverage -a "$OUT/${m}_win.bed" -b "$OUT/snps.sorted.bed" -sorted -counts -g "$OUT/genome.txt" > "$OUT/${m}_snp.txt" 2>/dev/null
  # acc.txt: 5 win cols + [nOv coveredbp winlen frac] -> frac=$9 ; snp.txt: 5 win cols + count -> count=$6
  paste "$OUT/${m}_acc.txt" "$OUT/${m}_snp.txt" | awk -v amin="$ACC_MIN" -v dmin="$DENS_MIN" 'BEGIN{OFS="\t"}
    {chrom=$1;stat=$4;pos=$5;frac=$9;snpc=$15;wlen=$3-$2; dens=snpc/(wlen/1000.0);
     tot++; if(frac>=amin && dens>=dmin){print stat,chrom,pos,frac,dens; kept++}}
    END{printf "  %s: %d/%d windows pass (acc>=%.0f%%, dens>=%.0f/kb)\n",m,kept,tot,100*amin,dmin > "/dev/stderr"}' m="$m" \
    | sort -gr > "$OUT/${m}_callable.tsv"
}
rerank sf2 25000
rerank raisd 12500
echo "=== TOP callable SF2 (LR  chr  pos  acc%  SNP/kb) ==="; head -8 "$OUT/sf2_callable.tsv" | awk -F'\t' '{printf "  %-8.1f chr%s:%-10s acc=%.0f%% %.0f/kb\n",$1,$2,$3,100*$4,$5}'
echo "=== TOP callable RAiSD (mu  chr  pos  acc%  SNP/kb) ==="; head -8 "$OUT/raisd_callable.tsv" | awk -F'\t' '{printf "  %-8.0f chr%s:%-10s acc=%.0f%% %.0f/kb\n",$1,$2,$3,100*$4,$5}'
