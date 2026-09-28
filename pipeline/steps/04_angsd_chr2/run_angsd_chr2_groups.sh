#!/usr/bin/env bash
# chr2 ANGSD by-group diversity (baker-633): metapop (baker pooled) + karyotype AA/AB/BB.
# Folded SFS (-fold 1) -> pi/theta_w/Tajima's D per group; realSFS 2D -> FST for karyotype pairs.
# Sites = chr2 accessible (from accessible_sites.bed). SHARED MACHINE: 4 groups x -P8 ~= 32 cores.
set -uo pipefail
source /sietch_colab/data_share/illex/popgen_data/analysis/config.sh
ANGSD=/home/ssmall/programs/angsd/angsd
REALSFS=/home/ssmall/programs/angsd/misc/realSFS
THETASTAT=/home/ssmall/programs/angsd/misc/thetaStat
SAMBIN=/home/ssmall/miniforge3/envs/bioinfo-buddy/bin
D=$ANA/steps/04_angsd_chr2
BL=$D/bamlists; SAF=$D/saf; TH=$D/thetas; FST=$D/fst; SITES=$D/chr2.accessible.sites
mkdir -p "$SAF" "$TH" "$FST" "$D/logs"
CHR=2; THREADS=8

# --- build chr2 accessible ANGSD sites (per-position) + index (once) ---
if [[ ! -s "${SITES}.idx" ]]; then
  echo "[$(date)] building chr2 accessible sites file"
  awk '$1=="2"{for(p=$2+1;p<=$3;p++) print "2\t"p}' "$ACCESS" > "$SITES"
  echo "  $(wc -l < "$SITES") sites; indexing"
  "$ANGSD" sites index "$SITES"
fi

run_saf() {  # $1 = group name
  local g=$1 bam="$BL/${1}_bamlist.txt" out="$SAF/${1}.${CHR}"
  [[ -s "$out.saf.idx" ]] && { echo "[$g] SAF exists"; return 0; }
  echo "[$(date)] [$g] SAF ($(wc -l < "$bam") bams)"
  "$ANGSD" -b "$bam" -sites "$SITES" -r "$CHR" -anc "$REF" -ref "$REF" -out "$out" \
    -doSaf 1 -doCounts 1 -GL 1 -P "$THREADS" \
    -minInd 1 -setMinDepthInd 1 -minQ 20 -minMapQ 20 -remove_bads 1 -only_proper_pairs 1 \
    > "$D/logs/saf.${g}.log" 2>&1
  "$REALSFS" "$out.saf.idx" -fold 1 -P "$THREADS" -tole 1e-8 -maxIter 200 > "$TH/${g}.${CHR}.sfs" 2>>"$D/logs/saf.${g}.log"
  "$REALSFS" saf2theta "$out.saf.idx" -sfs "$TH/${g}.${CHR}.sfs" -outname "$TH/${g}.${CHR}" -fold 1 2>>"$D/logs/saf.${g}.log"
  "$THETASTAT" do_stat "$TH/${g}.${CHR}.thetas.idx" -outnames "$TH/${g}.${CHR}.avg" 2>>"$D/logs/saf.${g}.log"
  "$THETASTAT" do_stat "$TH/${g}.${CHR}.thetas.idx" -win 50000 -step 10000 -outnames "$TH/${g}.${CHR}.win50k" 2>>"$D/logs/saf.${g}.log"
  echo "[$(date)] [$g] DONE"
}
export -f run_saf; export ANGSD REALSFS THETASTAT SITES SAF TH BL REF D CHR THREADS

echo "[$(date)] === SAF/theta for 4 groups (parallel) ==="
printf '%s\n' bakerMeta AA AB BB | xargs -P 4 -I{} bash -c 'run_saf "$@"' _ {}

echo "[$(date)] === FST for karyotype pairs ==="
fst_pair() {  # $1 $2
  local p1=$1 p2=$2 o="$FST/${1}_${2}.${CHR}"
  "$REALSFS" "$SAF/${p1}.${CHR}.saf.idx" "$SAF/${p2}.${CHR}.saf.idx" -fold 1 -P "$THREADS" > "$o.2dSFS" 2>"$D/logs/fst.${p1}_${p2}.log"
  "$REALSFS" fst index "$SAF/${p1}.${CHR}.saf.idx" "$SAF/${p2}.${CHR}.saf.idx" -sfs "$o.2dSFS" -fstout "$o" -fold 1 2>>"$D/logs/fst.${p1}_${p2}.log"
  "$REALSFS" fst stats "$o.fst.idx" > "$o.fst.txt" 2>>"$D/logs/fst.${p1}_${p2}.log"
  echo "[$(date)] FST ${p1} vs ${p2}: $(cat "$o.fst.txt")"
}
export -f fst_pair
for pair in "AA BB" "AA AB" "AB BB"; do fst_pair $pair; done

echo "[$(date)] === chr2-average theta summary (pi/thetaW/TajD per group) ==="
for g in bakerMeta AA AB BB; do
  f="$TH/${g}.${CHR}.avg.pestPG"
  [[ -s "$f" ]] && awk -v g="$g" 'NR==2{print g": "$0}' "$f"
done | tee "$D/chr2_group_theta_summary.txt"
echo "[$(date)] ALL DONE"
