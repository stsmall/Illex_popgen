#!/usr/bin/env bash
# FST between arrangements — subsample 40/arrangement + ACCESSIBLE SNP sites only,
# region 2:60-80Mb. Shrinks the 2D-SFS to 81x81 over ~1.85M SNPs -> fast/stable.
set -uo pipefail
ANA=/sietch_colab/data_share/illex/popgen_data/analysis
REF=/sietch_colab/data_share/illex/andy_links/Illex_F24.primary.clean.fa
ACC=/sietch_colab/data_share/illex/popgen_data/degenotate_illex/accessible_sites.bed
BBIN=/home/ssmall/miniforge3/envs/bioinfo-buddy/bin
ANGSD=/home/ssmall/programs/angsd/angsd
REALSFS=/home/ssmall/programs/angsd/misc/realSFS
D=$ANA/steps/04_angsd_chr2; BL=$D/bamlists; SAF=$D/saf_sub; F=$D/fst_sub
VCF=$ANA/steps/00_callset/filtered/2/variants_filt.vcf.gz
REG=2:60000000-80000000; TH=8
mkdir -p "$SAF" "$F" "$D/logs"

# 1. subsample 40 per arrangement (baker panmictic -> any 40 equivalent; deterministic)
for g in AA AB BB; do head -40 "$BL/${g}_bamlist.txt" > "$BL/${g}40_bamlist.txt"; done
echo "subsampled: AA40=$(wc -l <$BL/AA40_bamlist.txt) AB40=$(wc -l <$BL/AB40_bamlist.txt) BB40=$(wc -l <$BL/BB40_bamlist.txt)"

# 2. ACCESSIBLE biallelic SNP sites in 60-80Mb -> ANGSD sites + index
SITES=$D/inv_snp_accessible.sites
if [[ ! -s "${SITES}.idx" ]]; then
  "$BBIN/bcftools" view -r "$REG" -m2 -M2 -v snps "$VCF" 2>/dev/null \
    | "$BBIN/bcftools" query -f '%CHROM\t%POS\n' \
    | awk 'BEGIN{OFS="\t"}{print $1,$2-1,$2}' \
    | "$BBIN/bedtools" intersect -a - -b "$ACC" -u \
    | awk 'BEGIN{OFS="\t"}{print $1,$3}' > "$SITES"
  echo "accessible SNP sites in $REG: $(wc -l < "$SITES")"
  "$ANGSD" sites index "$SITES"
fi

# 3. SAF per subsampled arrangement (region, accessible SNP sites, folded)
saf_one() {
  local g=$1 out="$SAF/${g}40.inv"
  [[ -s "$out.saf.idx" ]] && return 0
  "$ANGSD" -b "$BL/${g}40_bamlist.txt" -sites "$SITES" -r "$REG" -anc "$REF" -ref "$REF" -out "$out" \
    -doSaf 1 -doCounts 1 -GL 1 -P "$TH" -minInd 1 -setMinDepthInd 1 -minQ 20 -minMapQ 20 \
    -remove_bads 1 -only_proper_pairs 1 > "$D/logs/saf_sub.${g}.log" 2>&1
}
for g in AA AB BB; do echo "[$(date)] SAF ${g}40"; saf_one "$g"; done

# 4. 2D-SFS + FST per pair (small matrix -> fast)
for pair in "AA BB" "AA AB" "AB BB"; do
  set -- $pair; p1=$1; p2=$2; o="$F/${p1}_${p2}"
  "$REALSFS" "$SAF/${p1}40.inv.saf.idx" "$SAF/${p2}40.inv.saf.idx" -fold 1 -P "$TH" > "$o.2dSFS" 2>"$D/logs/fstsub.${p1}_${p2}.log"
  "$REALSFS" fst index "$SAF/${p1}40.inv.saf.idx" "$SAF/${p2}40.inv.saf.idx" -sfs "$o.2dSFS" -fstout "$o" -fold 1 2>>"$D/logs/fstsub.${p1}_${p2}.log"
  "$REALSFS" fst stats "$o.fst.idx" > "$o.fst.txt" 2>>"$D/logs/fstsub.${p1}_${p2}.log"
  echo "[$(date)] FST ${p1} vs ${p2} (40/arr, accessible SNPs, 60-80Mb): $(cat "$o.fst.txt")"
done
echo "[$(date)] DONE"
