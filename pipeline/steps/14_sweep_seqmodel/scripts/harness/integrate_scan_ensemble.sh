#!/usr/bin/env bash
# Integrate the diploSHIC empirical scan (100kb classified windows) with the SF2 (CLR) + RAiSD (mu)
# + Bvalcalc (B) ensemble, plus callability + CDS density. One row per diploSHIC window:
#   chrom start end predClass p_neut p_lsoft p_lhard p_soft p_hard maxSF2_LR maxRAiSD_mu meanB callable_frac cds_frac
# Ranks concordant sweep candidates (diploSHIC hard/soft AND callable AND external support) + chr43:72-73Mb.
set -uo pipefail
B=/sietch_colab/data_share/illex/popgen_data
O=$B/analysis/steps/14_sweep_seqmodel/results/empirical_scan
E=$B/analysis/steps/17_ensemble_scan/results
BT=/home/ssmall/miniforge3/envs/annotation-expression-buddy-barrnap_trnascan/bin/bedtools
FAI=/sietch_colab/data_share/illex/assembly/current_assembly_version/Illex_F24.primary.clean.fa.fai
ACC=$B/degenotate_illex/accessible_sites.bed
CDS=$B/analysis/steps/16_bgspy/inputs/annot_cds_auto.bed
INT=$O/integration; mkdir -p "$INT"; cd "$INT"
export LC_ALL=C
CH="1 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40 41 43 44 45"

awk 'BEGIN{OFS="\t"}{print $1,$2}' "$FAI" | sort -k1,1 > genome.txt
# diploSHIC 100kb windows + class/probs (cols 4-9: predClass p_neut p_lsoft p_lhard p_soft p_hard)
for c in $CH; do tail -n +2 "$O/preds/chr$c.preds" 2>/dev/null; done \
  | awk 'BEGIN{OFS="\t"}{print $1,$2-1,$3,$5,$6,$7,$8,$9,$10}' | sort -k1,1 -k2,2n > win.bed
# ensemble tracks as point/interval BEDs
for c in $CH; do tail -n +2 "$E/sf2/$c.clr" 2>/dev/null | awk -v c="$c" 'BEGIN{OFS="\t"}{p=int($1); if(p<1)p=1; print c,p-1,p,$2}'; done | sort -k1,1 -k2,2n > sf2.bed
for c in $CH; do grep -v '^//' "$E/raisd/RAiSD_Report.$c" 2>/dev/null | awk -v c="$c" 'BEGIN{OFS="\t"} NF==7{print c,$1-1,$1,$7}'; done | sort -k1,1 -k2,2n > raisd.bed
grep -v '^#' "$E/illex_Bmap.csv" | awk -F',' 'BEGIN{OFS="\t"} $3!="nan" && $3!=""{print $1,$2-1,$2+999,$3}' | sort -k1,1 -k2,2n > bval.bed
sort -k1,1 -k2,2n "$ACC" > acc.sorted.bed
sort -k1,1 -k2,2n "$CDS" > cds.sorted.bed

# map ensemble stats onto windows (order preserved = win.bed order)
"$BT" map -a win.bed -b sf2.bed   -c 4 -o max  -g genome.txt -null NA > m1
"$BT" map -a m1      -b raisd.bed  -c 4 -o max  -g genome.txt -null NA > m2
"$BT" map -a m2      -b bval.bed   -c 4 -o mean -g genome.txt -null NA > m3
"$BT" coverage -a win.bed -b acc.sorted.bed -sorted -g genome.txt | awk '{print $NF}' > acc.frac
"$BT" coverage -a win.bed -b cds.sorted.bed -sorted -g genome.txt | awk '{print $NF}' > cds.frac
paste m3 acc.frac cds.frac > body.tsv
{ printf 'chrom\tstart\tend\tpredClass\tp_neut\tp_lsoft\tp_lhard\tp_soft\tp_hard\tmaxSF2_LR\tmaxRAiSD_mu\tmeanB\tcallable_frac\tcds_frac\n'; cat body.tsv; } > integrated.tsv
echo "wrote $INT/integrated.tsv ($(($(wc -l < integrated.tsv)-1)) windows)"

# ---- reports ----
echo "=== predClass distribution (all windows) ==="
awk -F'\t' 'NR>1{c[$4]++} END{for(k in c) printf "  %-12s %d\n",k,c[k]}' integrated.tsv | sort -k2,2nr
echo "=== TOP sweep candidates: diploSHIC hard/soft, callable>=0.5, ranked by sweep prob (with ensemble support) ==="
awk -F'\t' 'NR>1 && ($4=="hard"||$4=="soft") && $13>=0.5 {sp=($4=="hard"?$9:$8);
  printf "  chr%s:%d-%d %-5s p=%.3f  SF2=%s RAiSD=%s B=%s acc=%.2f cds=%.2f\n",$1,$2,$3,$4,sp,$10,$11,$12,$13,$14}' integrated.tsv \
  | sort -t= -k2,2gr | head -25
echo "=== chr43:72-73Mb (the SF2 lead region) ==="
awk -F'\t' 'NR>1 && $1=="43" && $3>72000000 && $2<73000000 {
  printf "  chr43:%d-%d %-10s p_hard=%.3f p_soft=%.3f SF2=%s RAiSD=%s B=%s acc=%.2f\n",$2,$3,$4,$9,$8,$10,$11,$12,$13}' integrated.tsv
echo "INTEGRATION_DONE"
