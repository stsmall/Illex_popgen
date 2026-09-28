# Illex_F24 annotation v3 (2026-09-28)

`Illex_F24.gene_lnc_pseudo.func.fix.sq3.FINAL.v3.gff3` = v2 unchanged, plus:

1. **14 `protein_match` features on chromosome 2** (`invrecov_01`..`invrecov_14`, source `DIAMOND_blastx`): protein-coding loci
   recovered by DIAMOND blastx of unannotated, non-repeat sequence in the chr2 inversion against seven cephalopod proteomes,
   function transferred from SwissProt. Interval-level evidence only (no exon model), so gene counts are unchanged from v2.
   Source table: msinv `results/illex/missed_genes_inversion.csv` (class `gene_candidate`).
2. **Note attribute on `LOC_0000SQ312948` and `LOC_0000SQ312965`** (chr16:35.3 Mb): blastx of the transcript against nr matches a
   cephalopod-specific uncharacterised protein family (Euprymna XP_079992720.1, 40% identity, E=2e-53).

`Illex_F24.v3_additions.gff3` contains only the added and modified records, and `Illex_F24.gene_lnc_pseudo.func.fix.sq3.FINAL.v3.gff3.gz` is the full annotation.
