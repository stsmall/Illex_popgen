# Illex manuscript: rules for editing

This repository holds the manuscript (`main.tex`), supplement (`supplement.tex`), bibliography, figures, and the
scripts behind every result. These rules apply to anyone, human or AI, editing the text.

## Writing style

- Write full sentences that explain the data. State what was found and what it means.
- Do not use semicolons. Split the clause into its own sentence, or join it with "and", "but", or "because".
- Do not use colons in prose. Colons are allowed only in genomic coordinates (chr2:60.54 Mb), ratios, and table cells.
- Do not use em-dashes (`---`). Use a comma, parentheses, or a new sentence. The en-dash `--` is correct in numeric
  ranges (60.54--79.50 Mb) and stays.
- Do not use `\emph{}` or italics for emphasis. Italics are for gene names, species names, and variables only.
  Bold is for caption titles and panel letters only.
- Keep parentheses for citations, figure references, and short numeric detail, not for extra clauses.
- Present the final, correct analysis. Do not describe fixes, reruns, or corrections ("after masking we removed 14
  regions"). Report the corrected result as the result.
- Keep the main-text Methods brief. Put parameters, thresholds, and command options in the Supplementary Methods.

## Facts must match the code

- Every parameter, threshold, sample count, and software option in the text must match the script that produced it.
  The scripts are in `pipeline/` (analysis) and `scripts/` (figures and tables). The chromosome-2 inversion analysis is
  in https://github.com/stsmall/msinv under `illex/` and `results/illex/`.
- Before changing a number or method description, find it in the code and cite the file in the commit message.
  If it cannot be found, do not state it as fact. Ask the authors.
- Do not describe a result from memory or from an earlier draft. Re-derive it or quote the output file.
- Every `\cite` key must exist in `references.bib`, and every bib entry must be checked against CrossRef or the
  publisher, with a DOI.

## Key facts that are easy to get wrong

- All 633 individuals were sequenced by Baker et al. 2025 (ENA PRJEB101095). Baker analysed 540 of them.
- The callset is GATK HaplotypeCaller joint genotyping (gVCF, all sites retained) run through grenepipe.
- Diversity statistics come from ANGSD run directly on the BAMs (`-GL 1`, the SAMtools model), not from the VCF.
- Windowed per-arrangement statistics, Hudson FST, LD, and PCA use pg_gpu on called genotypes. Per-SNP geographic FST
  uses plink2 (Weir--Cockerham).
- The accessibility mask is the sequence outside Earl Grey repeats (1.43 Gb). For the sweep scan and gene-set analyses
  it is extended by the 22 anomalous-diversity blocks (Supplementary Table S9).
- The demography was fitted before the 22 blocks were masked. Do not refit it.
- *Cpes* (chr1:85.5 Mb) and *GALT* (chr25:65.0 Mb) are existing annotated genes. Fine-scale localisation reassigned
  sweep targets to them. They are not newly discovered genes.
- The annotation is `Illex_F24` version 3 (`data/annotation/`). The 14 loci recovered in the inversion are
  `protein_match` features and are not counted among the inversion's 89 gene models.

## Supplement numbering

The supplement is a separate document, so the main text cites it with hard-coded labels (Supplementary Fig.~S12).
After adding or removing a supplementary figure or table, check `supplement.aux` and update every hard-coded label in
`main.tex`. Never `\cref` a main-text label from `supplement.tex`. Write "main-text Fig.~3" instead.

## Building and checking

```
make            # builds main.pdf and supplement.pdf
```

Before committing, confirm that both documents build with no undefined references, no LaTeX errors, and no overfull
boxes wider than 30 pt, and that the text contains no semicolons, prose colons, em-dashes, or `\emph`.
