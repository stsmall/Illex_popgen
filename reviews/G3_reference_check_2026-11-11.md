# Reference check for the G3 mock review (2026-11-11)

Companion to [`G3_review_2026-11-11.md`](G3_review_2026-11-11.md). Checked against `main.pdf` and `supplement.pdf` at
commit `8a87bd4`. Each cited work has a link (DOI, arXiv or repository). For each citation the report gives the
manuscript location (main text page:line, supplement page and paragraph), the claim, a verbatim passage from the cited
work with its location, and a rating (SUPPORTS / PARTIAL / DOES NOT SUPPORT / CANNOT CHECK). Where only the abstract
was reachable, the entry says so. Tick the box in the action list when an item is fixed.

## Action list

Status 2026-11-11 (manuscript commit after `72fc198`): ticked items are fixed. Open: `small2026` accession or preprint (the "recently made available" wording is removed), the "(needs reference)" placeholder, the msprime 1.4.1 vs 1.4.2 conflict and the SQANTI3, LiftOn and Pfam-A versions (not verifiable on the analysis machine), and the CANNOT CHECK items. The StringTie item needs no change, because the text never says StringTie2 (only the cited paper's title does).

- [x] `msinv`: bib, Data availability and Table S1 cite v0.3.5 / commit 35412bf. Update to the 0.4.1 release and commit that produced the ABC-SMC and NPE results.
- [x] `matuszewski2015` (main 12:505–506) DOES NOT SUPPORT "a moving optimum can hold an allele at an intermediate frequency". Replace the citation or the claim.
- [x] `baker2025` (main 17:744–746, 2:74–75): Baker retained 616 after QC then filtered to 540, never states 633 and never attributes missingness to the reference. Its SV caveat is general. Rephrase both.
- [x] `adrion2020relernn` (main 8:348–349, 18:799): ReLERNN does not use LD summaries. Say "from patterns of linkage in genotype data".
- [x] `kelleher2018` (main 18:804–805) is the tree-sequence-recording paper, not msprime. Cite `kelleher2016` or `baumdicker2022` for msprime.
- [x] `albertin2015` (main 11:462): Octopus genome. Drop it for the Euprymna protein family (keep `albertin2022`).
- [x] `mathger2010`, `ramirez2015` (main 3:103–106): cuttlefish and octopus under "In squids". Change to "In cephalopods".
- [x] `stefansson2005` (main 12:494): does not state a recent increase of 17q21 H2. Rephrase ("under ongoing positive selection").
- [x] `navarro2000`, `peischl2013` (main 6:219–220): Peischl is a simulation-algorithm paper. Cite an analytic source for diversity scaling with arrangement frequency (e.g. Navarro et al. 2000 equations, or Guerrero et al. 2012).
- [x] `waples1998`, `hauser2008`, `gagnaire2015` (main 2:39–41): none says marine datasets were limited to mtDNA or few loci, and Hauser & Carvalho stress extensive marine structure. Rephrase or recite.
- [ ] `small2026` is "in preparation" yet cited for a genome "recently made available", the 4.14 Gb assembly, k-mer heterozygosity, long reads and Z lncRNAs. Add an assembly accession or preprint DOI.
- [ ] Main 2:38–39 "(needs reference)" placeholder.
- [x] `degenotate` URL 404: use https://github.com/harvardinformatics/degenotate.
- [x] `poplin2018gatk`: bioRxiv first posted 2017.
- [x] `chao2024lifton`: add Genome Research 35(2):311–325 (issue 2025).
- [x] `talts2018`: add DOI 10.48550/arXiv.1804.06788. `lemos2023`: add PMLR 202:19256–19273.
- [x] `korunes2019` (supp p. 3): paper reports 188–390 bp tracts, not 300 bp. Say "300 bp, within the 188–390 bp estimated by Korunes and Noor".
- [x] `miyata1987` (supp p. 7, Table S12 caption): X/Y/autosome derivation. Say "following Miyata et al. (1987)" for the Z version.
- [x] `kelleher2016` cited for msinv in Table S1: use `kelleher2018` (tskit) for tree-sequence output, which is already on the tskit row.
- [ ] Version conflicts: msprime 1.4.1 (supp p. 4) vs 1.4.2 (Table S1). StringTie2 (supp p. 3) vs StringTie 3.0.0 (Table S1). Add versions for SQANTI3, LiftOn and the Pfam-A release.
- [x] Table S1 caption draft text "Citations are pointers to confirm at final formatting". The comment in software_table.tex says it is included by main.tex but supplement.tex inputs it.
- [x] `baumdicker2022`: author "Kretzschumar" copies a publisher typo (Kretzschmar).
- [x] Uncited bib entries: `keightley2018estsfs`, `schrider2017`. Remove or cite.
- [ ] CANNOT CHECK, confirm by hand: `arkhipkin2015`, `odor1998`, `schweinsberg2003`, `small2026`, `nei1986` (abstract only), `ohta1971`.

---

# Part 1. Main-text citations


Scope: all 104 keys cited in main.tex (\citep, \citet, \citealt). Page:line refer to the built main.pdf (printed page number : left-margin line number). Metadata was checked against Crossref for every DOI (script output agreed on title, authors, year, volume and pages, apart from the items listed below). Full text came from Europe PMC/PMC (OA subset or efetch) or the publisher page. Where only the abstract could be reached, the entry says so.

## Main: Summary

Ratings (104 keys): **SUPPORTS 87 | PARTIAL 12 | DOES NOT SUPPORT 1 | CANNOT CHECK 4**

### Problems

**Metadata / bibliography**
1. `degenotate`: the howpublished URL `https://github.com/gwct/degenotate` returns HTTP 404. The repository is `https://github.com/harvardinformatics/degenotate`. Version 1.3 is real (released 2024-02-05).
2. `msinv`: the bib and Data-availability text both give "v0.3.5, commit 35412bf". That commit exists (2026-10-29, "illex: parametric bootstrap"), sits 13 commits after the v0.3.5 tag, and its pyproject still says 0.3.5. msinv `main` is now 0.4.1, which adds fixes for crossover suppression and coalescence in heterokaryotypes. The latest Illex NPE/ABC work is at commit 58a695e (2026-11-10). Check that the cited version/commit is the one that produced the reported ABC-SMC/NPE results; it probably needs updating.
3. `small2026` is `@unpublished`, "Manuscript in preparation", and its title is about Z dosage compensation. The text cites it for: a genome "recently made available" (p3:111-112), assembly size and chromosome count (p17:747), 1.06% k-mer heterozygosity (p8:333-334), long reads from a second individual (p9:380), and Z lncRNAs/Zfest (p11:475-476). None of this can be checked. Saying the genome was "recently made available" while citing a manuscript in preparation is inconsistent. Give an assembly accession or a preprint DOI if one exists.
4. `poplin2018gatk`: bioRxiv 10.1101/201178 was first posted 2017-11-14 (Crossref year 2017). The bib's 2018 is a later version. Minor.
5. `haller2019slim`: the key name is stale. The entry is correctly SLiM 5 (MBE 43:msaf313, 2026) and renders as "Haller et al., 2026". Cosmetic only.
6. No DOI by design (acceptable): `durkan2019` (NeurIPS; arXiv:1906.04032 verified), `odor1998` (FAO Tech. Paper 376; URL resolves, HTTP 200), `coindetii2025` (GenBank GCA_977009265.1 verified: Illex coindetii, xcIllCoin1.1, chromosome level, Genoscope, released 2025-11-08, PRJEB89560, all matching the bib), `msinv`, `degenotate`, `small2026`.
7. Every DOI resolves to the intended work. No wrong DOIs, years, volumes or pages beyond item 4. Crossref title differences for bhatia2013, degiorgio2016, hart2020entap and hermisson2005 are markup or subtitle only.

**Missing citation in text**
8. p2:38-39: "There are examples of large-Ne species with no geographic structure in terrestrial systems (needs reference)". This placeholder has no citation.

**Unsupported or weakly supported claims**
9. `matuszewski2015` (p12:505-506), **DOES NOT SUPPORT**: the paper is cited for "a moving selective optimum can hold an allele at an intermediate frequency". It models the distribution of adaptive substitutions (fixations) from standing variation under a moving optimum and says nothing about maintaining alleles at intermediate frequency.
10. `baker2025`, PARTIAL: (a) p17:744-746 says Baker analysed 540 of the 633 "after removing individuals with more than 25% missing genotypes, a consequence of mapping to the distantly related Thysanoteuthis reference". Baker retained 616 after unspecified QC, then applied the >25%-missing filter to reach 540. Baker never states 633 and never attributes the missingness to the reference. (b) p2:74-75 says Baker "explicitly acknowledged reference divergence as a constraint on their detection" of SVs. Baker gives a general caution about the distant reference and says remapping could reveal signals in non-conserved regions, but does not tie this explicitly to SV detection. Everything else checks out: 540 specimens, >40 sites, ~51 My divergence, 46% alignment, no large SVs, PRJEB101095, 2022-2023, 11 NAFO divisions.
11. `kelleher2018` (p18:804-805) is cited as the msprime reference for recapitation. It is the tree-sequence-recording (pedigree recording) paper. The msprime citation should be `kelleher2016` or `baumdicker2022`, which are already in the bib. `haller2019recap` supports recapitation.
12. `albertin2015` (p11:462) is cited, with albertin2022, for an uncharacterised protein family in the bobtail squid *Euprymna*. Albertin 2015 is the *Octopus bimaculoides* genome. Albertin 2022 does supply the Euprymna chromosome assembly, so drop 2015 here.
13. `mathger2010`, `ramirez2015` (p3:103-106) are cited for "In squids ... dermal opsins and light-responsive chromatophore systems". Mäthger is cuttlefish and Ramirez & Oakley is octopus; only Kingston 2015 includes a squid. Say "cephalopods" instead.
14. `stefansson2005` (p12:494) is cited for the human 17q21 inversion "whose increase is more recent". Stefansson shows the H1/H2 lineages diverged up to 3 My ago and that H2 is under ongoing positive selection. The paper does not directly state a recent increase. Rephrase.
15. `navarro2000`, `peischl2013` (p6:219-220) are cited for "without gene flux the neutral diversity of each [arrangement] scales with its frequency". The Navarro abstract only mentions equilibrium expressions (full text blocked). Peischl 2013 is an SMC simulation algorithm, so it is a weak source for that analytic statement.
16. `adrion2020relernn` (p8:348-349, p18:799) says ReLERNN "infers the rate from linkage disequilibrium". The ReLERNN abstract says the opposite: "Rather than use summaries of linkage disequilibrium as its input, ReLERNN takes columns from a genotype alignment". Suggested wording: "from patterns of linkage in genotype data".
17. `waples1998`, `hauser2008`, `gagnaire2015` (p2:39-41) are cited for "marine ... geographic datasets have largely been limited to mitochondrial genomes or a few nuclear loci". These works discuss weak structure and marker limits in high-gene-flow marine species, but none states that specific point. Hauser & Carvalho actually stress that genetics revealed *extensive* structure in many marine species.

**Missing / uncited keys**
18. Cited in main.tex but missing from references.bib: **none**.
19. In references.bib but cited nowhere (main.tex, supplement.tex or the \input tables): **`keightley2018estsfs`, `schrider2017`**.

---

## Main: Per-key detail

Format: **key**. Metadata | manuscript location and claim | supporting passage (location) | rating.

**adrion2020relernn** [doi:10.1093/molbev/msaa038](https://doi.org/10.1093/molbev/msaa038). MBE 37:1790-1808 (2020), DOI ok. | p8:348-349 "Recombination was estimated with ReLERNN ... a deep-learning method that infers the rate from linkage disequilibrium"; p18:797-799. | "Here, we describe recombination landscape estimation using recurrent neural networks (ReLERNN), a deep learning method for estimating a genome-wide recombination map ... Rather than use summaries of linkage disequilibrium as its input, ReLERNN takes columns from a genotype alignment" (abstract). | **PARTIAL**: right method, but it does not take LD summaries as input; fix the wording.

**alachiotis2018raisd** [doi:10.1038/s42003-018-0085-8](https://doi.org/10.1038/s42003-018-0085-8). Commun Biol 1:79 (2018), ok. | p3:128 complementary sweep methods; p9:407 "independent LD-footprint method, RAiSD"; p18:809. | "RAiSD ... implements a novel ... parameter-free detection mechanism that relies on multiple signatures of a selective sweep via the enumeration of SNP vectors" (abstract). | **SUPPORTS** (RAiSD combines SFS, LD and diversity signatures; "LD-footprint" is a simplification).

**albertin2015** [doi:10.1038/nature14668](https://doi.org/10.1038/nature14668). Nature 524:220-224, ok. | p2:43 large, repeat-rich, reorganized cephalopod genomes; p3:102 sophisticated nervous systems; p11:462 Euprymna uncharacterised protein family; p13:569 centralized nervous system. | "Nearly 45% of the assembled genome is composed of repetitive elements" (Results); "we found evidence for large-scale genomic rearrangements that are closely associated with transposable element expansions" (abstract). | **PARTIAL**: supports the genome and nervous-system claims, not the Euprymna claim at p11:462 (this is an Octopus genome paper).

**albertin2022** [doi:10.1038/s41467-022-29748-w](https://doi.org/10.1038/s41467-022-29748-w). Nat Commun 13:2427, ok. | p2:43; p3:102; p11:462 (Euprymna). | "The genomes of the soft-bodied (coleoid) cephalopods are highly rearranged relative to other extant molluscs ... feature multi-megabase, tandem arrays of genes associated with brain development" (abstract); "we also produced chromosome-scale assemblies of the Hawaiian bobtail squid, Euprymna scolopes" (Results). | **SUPPORTS**.

**arkhipkin2015** [doi:10.1080/23308249.2015.1026226](https://doi.org/10.1080/23308249.2015.1026226). Rev Fish Sci Aquacult 23:92-252, ok. | p2:51, 55 ommastrephid life history; I. illecebrosus range and fisheries. | Paywalled; the OA repository mirror blocked automated access; Crossref has no abstract. | **CANNOT CHECK** (a world squid-fisheries review is topically appropriate; the ">2,500 km" figure is in any case supported by baker2025's abstract).

**arnason2023** [doi:10.7554/eLife.80781](https://doi.org/10.7554/eLife.80781). eLife 12:e80781, ok. | p3:126 sweepstakes reproduction; p15:638-639 "in Atlantic cod, pervasive linked selection has been proposed to contribute to multiple-merger-like genealogies". | "The Durrett–Schweinsberg coalescent, in which selective sweepstakes result from recurrent and pervasive selective sweeps of new mutations, offers greater explanatory power. Our results show that models of sweepstakes reproduction and multiple-merger coalescents are relevant and necessary" (abstract). | **SUPPORTS**.

**baker2025** [doi:10.1093/icesjms/fsaf201](https://doi.org/10.1093/icesjms/fsaf201). ICES J Mar Sci 82(11):fsaf201, ok. | p2:55 range and fisheries; p2:60-68 "540 specimens collected from over 40 locations ... Thysanoteuthis rhombus ... approximately 51 million years ago ... mean alignment rate of merely ~46%"; p2:74-75 "did not find evidence of large structural variants and explicitly acknowledged reference divergence as a constraint on their detection"; p3:113-116, p4:144, Fig. 2 caption, Table 1; p17:741-746 "633 individuals ... (ENA project PRJEB101095), collected in 2022 and 2023 across eleven NAFO divisions ... That study analyzed 540 of them after removing individuals with more than 25% missing genotypes, a consequence of mapping to the distantly related Thysanoteuthis reference". | Abstract: "338 948 SNPs (9× coverage) for 540 specimens collected during 2022–2023 at >40 sites ... from Newfoundland, Canada to Cape Hatteras ... a lack of significant genetic structure across this geographic range; a distance >2500 km." Results: "Due to mapping against an outgroup reference genome, we found a mean alignment rate of 46% of reads (sd = 9%)." Discussion: "...is still approximately fifty-one million years diverged"; "We also do not find evidence of large-scale structural variants."; "Some caution may be warranted in interpreting these results"; remapping may reveal signals in "genomic regions that are not conserved between Ommastrephidae and diamondback squid". Methods: "we retained 616 individuals that passed quality filters"; filter to "remove individuals or loci with greater than 25% missing genotype information". Data availability: "ENA, Project: PRJEB101095". Table 1 lists NAFO 3K, 3N, 3O, 4R, 4S, 4T, 4W, 4X, 6A, 6B, 6C (publisher full text via web). | **PARTIAL**: the numbers are supported. The 633→540 account and the explicit SV-detection caveat are over-stated (Problem 10).

**beaumont2002** [doi:10.1093/genetics/162.4.2025](https://doi.org/10.1093/genetics/162.4.2025). Genetics 162:2025-2035, ok. | p18:793 neutral baseline fitted "by rejection ABC with regression adjustment". | "This is achieved by fitting a local-linear regression of simulated parameter values on simulated summary statistics" (abstract). | **SUPPORTS**.

**benjamini1995** [doi:10.1111/j.2517-6161.1995.tb02031.x](https://doi.org/10.1111/j.2517-6161.1995.tb02031.x). JRSS-B 57:289-300, ok. | p10:429; p18:815 BH-FDR. | "It calls for controlling the expected proportion of falsely rejected hypotheses — the false discovery rate" (abstract). | **SUPPORTS**.

**berdan2021** [doi:10.1371/journal.pgen.1009411](https://doi.org/10.1371/journal.pgen.1009411). PLoS Genet 17:e1009411, ok. | p12:499-500 "large non-recombining regions can accumulate recessive load". | "Recessive deleterious mutations accumulated on both arrangements with most of them being private to a given arrangement. This led to increasing overdominance, allowing for the maintenance of the inversion polymorphism" (abstract). | **SUPPORTS**.

**berdan2023** [doi:10.1111/jeb.14242](https://doi.org/10.1111/jeb.14242). J Evol Biol 36:1761-1782, ok. | p3:94-96 reconstructing age and frequency history with the barrier; p12:513 selection on arrangements as supergenes. | "the reduced recombination conferred by them strengthens the effects of drift and hitchhiking, which in turn can obscure signatures of selection ... here we review how different mechanisms of selection affect the evolution of inversions" (abstract). | **SUPPORTS** (general review, abstract-level).

**bhatia2013** [doi:10.1101/gr.154831.113](https://doi.org/10.1101/gr.154831.113). Genome Res 23:1514-1521, ok. | p4:160; p18:768 Hudson FST. | "However, we recommend a different estimator motivated by Hudson et al. (1992)." "These estimates are independent of sample sizes even when FST is not identical across populations." (Results). | **SUPPORTS**.

**boelts2025** [doi:10.21105/joss.07754](https://doi.org/10.21105/joss.07754). JOSS 10(108):7754, ok. | p18:788 NPE in sbi. | Title: "sbi reloaded: a toolkit for simulation-based inference workflows" (no abstract in Crossref). | **SUPPORTS** (software citation).

**buchfink2021diamond** [doi:10.1038/s41592-021-01101-x](https://doi.org/10.1038/s41592-021-01101-x). Nat Methods 18:366-368, ok. | p18:814 DIAMOND. | "we introduce an improved version of DIAMOND ... matching the sensitivity of the gold standard BLASTP" (abstract). | **SUPPORTS**.

**burt1991** [doi:10.1046/j.1420-9101.1991.4020259.x](https://doi.org/10.1046/j.1420-9101.1991.4020259.x). J Evol Biol 4:259-277, ok. | p16:683 reduction or loss of recombination in the heterogametic sex is widespread. | "achiasmate meiosis, which has evolved 25–30 times, is always restricted to the heterogametic sex in dioecious species" (abstract). | **SUPPORTS**.

**chang2015plink2** [doi:10.1186/s13742-015-0047-8](https://doi.org/10.1186/s13742-015-0047-8). GigaScience 4:7, ok (article number 7). | p18:770 per-SNP Weir–Cockerham FST in plink2. | "we are developing a second-generation codebase for PLINK" (abstract); the paper does not describe the FST estimator. | **SUPPORTS** (standard software citation).

**charlesworth1993** [doi:10.1093/genetics/134.4.1289](https://doi.org/10.1093/genetics/134.4.1289). Genetics 134:1289-1303, ok. | p9:360; p14:587 BGS reduces diversity in gene-dense regions. | "Selection against deleterious alleles maintained by mutation may cause a reduction in the amount of genetic variability at linked neutral sites" (abstract). | **SUPPORTS**.

**charlesworth2018inv** [doi:10.1534/genetics.117.300426](https://doi.org/10.1534/genetics.117.300426). Genetics 208:377-382, ok. | p3:86 inversion theory focused on local adaptation with migration. | "This model assumes that local selection maintains alleles at two or more loci, despite immigration ... We show that an inversion is favored because it prevents the breakdown of linkage disequilibrium generated by migration" (abstract). | **SUPPORTS**.

**charlesworth2021** [doi:10.1146/annurev-ecolsys-010621-044528](https://doi.org/10.1146/annurev-ecolsys-010621-044528). Annu Rev Ecol Evol Syst 52:177-197, ok. | p9:361; p14:588 BGS. | "We review the major population genetic processes ... selective sweeps, background selection, associative overdominance, and Hill-Robertson interference" (abstract). | **SUPPORTS**.

**chen2018fastp** [doi:10.1093/bioinformatics/bty560](https://doi.org/10.1093/bioinformatics/bty560). Bioinformatics 34:i884-i890, ok. | p17:748 trimming. | "We developed fastp as an ultra-fast FASTQ preprocessor with ... adapter trimming" (abstract). | **SUPPORTS**.

**coffing2025** [doi:10.1016/j.cub.2025.01.005](https://doi.org/10.1016/j.cub.2025.01.005). Curr Biol 35:931-939.e4, ok. | p2:44 cephalopod genomes; p4:141, p15:665 ZZ/ZO system; p16:678-680 female-biased mutation in two-spot octopus; "Z has been conserved as a single syntenic unit across coleoids for more than 300 million years". | "the female octopus has just one sex chromosome, consistent with a ZO karyotype, while the male has two (ZZ) ... the sex chromosome is of an ancient origin, before the radiation of extant cephalopods approximately 480 million years ago, and has been conserved to the present day in all cephalopod genomes available" (Summary). "for synonymous divergence α = 0.127 ... female mutation rate to be ~7.8 times that of the male" (Results). | **SUPPORTS** (paper says ~480 My for cephalopods; ">300 My" is conservative but consistent).

**coffing2026** [doi:10.1093/molbev/msag248](https://doi.org/10.1093/molbev/msag248). MBE 43(10):msag248, ok. | p2:44; p13:571-573 brain and eye development enriched among sweep targets in O. bimaculatus/bimaculoides; p15:664 "Sex-specific recombination has not been measured directly in any cephalopod"; p15:668-669 ratio 0.5 vs 1; p16:678-679 same conclusion (achiasmy, female-biased mutation). | "Protein-coding genes overlapping with selective sweeps were enriched for various functions that included many related to brain and eye development" (abstract). "To our knowledge, sex-specific recombination rates have not been directly measured in any cephalopod." (Discussion). "If females recombine at the same rate as males, then ρZ/ρA=(3/4)×(2/3)=1/2 ... if females lack recombination entirely (achiasmate meiosis) ... the ratio becomes ... =1" (Discussion). "showing further support for a sex-biased mutation rate and an achiasmate meiotic regime in the heterogametic sex" (Introduction). | **SUPPORTS**.

**coindetii2025** [link](https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_977009265.1/). GenBank GCA_977009265.1 verified (see Problem 6). | p3:114 I. coindetii outgroup. | NCBI Datasets: Illex coindetii, xcIllCoin1.1, Chromosome level, submitter Genoscope, 2025-11-08, PRJEB89560. | **SUPPORTS**.

**corbettdetig2012** [doi:10.1371/journal.pgen.1003056](https://doi.org/10.1371/journal.pgen.1003056). PLoS Genet 8:e1003056, ok. | p12:492-493 Drosophila inversions vary along geographic or seasonal gradients. | "Frequency clines of the most common D. melanogaster inversions are independently replicated on many continents ... Recurrent seasonal frequency shifts have been observed in numerous geographically diverse populations" (Introduction). | **SUPPORTS**.

**cranmer2020** [doi:10.1073/pnas.1912789117](https://doi.org/10.1073/pnas.1912789117). PNAS 117:30055-30062, ok. | p6:242-243 NPE. | "We review the rapidly developing field of simulation-based inference" (abstract). | **SUPPORTS**.

**czech2022grenepipe** [doi:10.1093/bioinformatics/btac600](https://doi.org/10.1093/bioinformatics/btac600). Bioinformatics 38:4809-4811, ok. | p17:751. | "grenepipe, an all-in-one Snakemake workflow ... from raw high-throughput sequencing data ... to genotype variant calls" (abstract). | **SUPPORTS**.

**degenotate** [GitHub](https://github.com/harvardinformatics/degenotate). Software, version 1.3 (2024) exists, **URL wrong** (Problem 1). | p18:820 codon degeneracy. | Repository harvardinformatics/degenotate ("annotation of degeneracy of coding sites"). | **SUPPORTS** (fix the URL).

**degiorgio2016** [doi:10.1093/bioinformatics/btw051](https://doi.org/10.1093/bioinformatics/btw051). Bioinformatics 32:1895-1897, ok. | p3:128; p9:408; p18:810 SweepFinder2. | "SweepFinder2, an extension of SweepFinder with increased sensitivity and robustness to the confounding effects of mutation rate variation and background selection" (abstract). | **SUPPORTS**.

**durkan2019** [arXiv:1906.04032](https://arxiv.org/abs/1906.04032). NeurIPS 2019, arXiv:1906.04032 verified (title and authors match). | p18:789 neural spline flows. | arXiv title: "Neural Spline Flows". | **SUPPORTS**.

**eldon2006** [doi:10.1534/genetics.105.052175](https://doi.org/10.1534/genetics.105.052175). Genetics 172:2621-2633, ok. | p3:126, p8:340-341, p15:634-635 sweepstakes and skewed offspring distribution. | "The complexity arises because individuals can potentially have very many offspring ... The data suggest the presence of rare reproduction events in which approximately 8% of the population is replaced by the offspring of a single individual" (abstract). | **SUPPORTS**.

**eldon2015** [doi:10.1534/genetics.114.173807](https://doi.org/10.1534/genetics.114.173807). Genetics 199:841-856, ok. | p8:340-341 sweepstakes reproduction would mimic expansion by inflating rare variants. | "An excess of singletons is a well-known characteristic of both population growth and multiple mergers" (abstract). | **SUPPORTS**.

**faria2019** [doi:10.1016/j.tree.2018.12.005](https://doi.org/10.1016/j.tree.2018.12.005). TREE 34:239-248, ok. | p2:79 supergenes; p11:482 clines and local adaptation. | "Inversion polymorphisms are maintained by different processes, for example, divergent or balancing selection" (abstract). | **SUPPORTS**.

**gagnaire2015** [doi:10.1111/eva.12288](https://doi.org/10.1111/eva.12288). Evol Appl 8:769-786, ok. | p2:39-41 marine datasets limited to mtDNA or few loci. | "In the marine realm, however, decades of molecular studies have met limited success in inferring genetic connectivity, due to the frequent lack of spatial genetic structure in species exhibiting high fecundity and dispersal capabilities" (abstract). | **PARTIAL**: supports weak marine structure, not the mtDNA/few-loci point.

**gilbert2020** [doi:10.1016/j.cub.2019.11.063](https://doi.org/10.1016/j.cub.2019.11.063). Curr Biol 30:101-107.e3, ok. | p3:90; p12:500 associative overdominance from recessive load. | "Purifying selection against recessive variants, however, can also lead to associative overdominance (AOD), due to an apparent heterozygote advantage at linked neutral loci" (abstract). | **SUPPORTS**.

**gone2** (Santiago et al. 2025). Nat Commun 16:6054, ok. | p8:339; p18:764 GONE2. | "We present two software tools ... GONE2, for inferring recent changes in Ne when a genetic map is available, and currentNe2" (abstract). | **SUPPORTS**.

**guerrero2012** [doi:10.1098/rstb.2011.0246](https://doi.org/10.1098/rstb.2011.0246). Phil Trans R Soc B 367:430-438, ok. | p3:95; p6:218-219; p13:554-555 arrangements behave as partially isolated demes. | "we develop coalescent models of neutral sites linked to an inversion polymorphism ... important determinants are the age of the inversion and the rate of gene flux between arrangements" (abstract). | **SUPPORTS**.

**haller2019recap** [doi:10.1111/1755-0998.12968](https://doi.org/10.1111/1755-0998.12968). Mol Ecol Resour 19:552-566, ok. | p18:805 recapitation with msprime. | "With a method that we call 'recapitation', we can leverage the information in the tree sequence to prepend a coalescent simulation of the burn-in period" (Introduction). | **SUPPORTS**.

**haller2019slim** (SLiM 5, 2026). MBE 43:msaf313, ok (key name stale). | p18:804 SLiM forward simulations. | "Here we present SLiM 5, a major extension of SLiM's capabilities" (abstract). | **SUPPORTS**.

**harris2018** [doi:10.1371/journal.pgen.1007859](https://doi.org/10.1371/journal.pgen.1007859). PLoS Genet 14:e1007859, ok. | p14:614-615 hard sweeps arise when adaptation depends on new or rare variants, particularly with a limited mutational target. | "well-fit demographic histories combined with rare hard sweeps serve as the more parsimonious explanation" (abstract); soft sweeps from recurrent mutation need "a large mutational target size, as the competing de novo beneficial mutations must have identical selective effects" (Introduction). | **SUPPORTS**.

**hart2020entap** [doi:10.1111/1755-0998.13106](https://doi.org/10.1111/1755-0998.13106). Mol Ecol Resour 20:591-604, ok. | p18:813. | "EnTAP ... was designed to improve the accuracy, speed, and flexibility of functional gene annotation" (abstract). | **SUPPORTS**.

**hauser2008** [doi:10.1111/j.1467-2979.2008.00299.x](https://doi.org/10.1111/j.1467-2979.2008.00299.x). Fish Fish 9:333-362, ok. | p2:39-41. | "they have contributed to the discovery of extensive genetic population structure in many marine species, overturning the notion of large, essentially homogenous marine populations" (abstract). | **PARTIAL**: about marine genetic approaches, but it stresses structure *found*, not limited loci.

**hedgecock2011** [doi:10.5343/bms.2010.1051](https://doi.org/10.5343/bms.2010.1051). Bull Mar Sci 87:971-1002, ok. | p3:126, p8:341, p15:635 sweepstakes reproductive success definition. | Title: "Sweepstakes reproductive success in highly fecund marine fish and shellfish: a review and commentary" (no abstract retrievable). | **SUPPORTS** (title-level; definitional citation).

**hermisson2005** [doi:10.1534/genetics.104.036947](https://doi.org/10.1534/genetics.104.036947). Genetics 169:2335-2352, ok. | p14:615 hard sweeps from new mutations when the mutational supply is limited. | "Adaptations from the standing genetic variation are favored if either the selective advantage is weak or the selection coefficient and the mutation rate are both high" (abstract). | **SUPPORTS**.

**jay2021** [doi:10.1038/s41588-020-00771-1](https://doi.org/10.1038/s41588-020-00771-1). Nat Genet 53:288-293, ok. | p3:82 Heliconius mimicry; p3:90; p12:496-498 inverted arrangement carries recessive deleterious variation that limits its frequency. | "inversions accumulate recessive deleterious mutations that generate frequency-dependent selection and promote their maintenance at intermediate frequency ... Inversions associate with sharply reduced viability when homozygous, which prevents them from replacing ancestral chromosome arrangements" (abstract). | **SUPPORTS**.

**johri2022** [doi:10.1371/journal.pbio.3001669](https://doi.org/10.1371/journal.pbio.3001669). PLoS Biol 20:e3001669, ok. | p10:403 per-sub-window diversity variation that a uniform-θ coalescent does not reproduce; p14:609-610 scans should accommodate features the null cannot reproduce. | "the starting point for any genomic analysis should be the construction of a biologically relevant baseline model ... This model should include mutation, recombination, and gene conversion ... purifying selection acting on functional regions and its effects on linked [sites]"; "both the mean of, and genomic heterogeneity in, mutation rates" (Recommendations). | **SUPPORTS**.

**jouganous2017moments** [doi:10.1534/genetics.117.200493](https://doi.org/10.1534/genetics.117.200493). Genetics 206:1549-1567, ok. | p18:763 moments. | "we propose a tractable model of ordinary differential equations for the evolution of allele frequencies" (abstract). | **SUPPORTS**.

**kapun2019** [doi:10.1111/mec.14871](https://doi.org/10.1111/mec.14871). Mol Ecol 28:1263-1282, ok. | p12:493 clines and seasonal fluctuation. | "Several studies have shown that inversion polymorphisms can form clines or fluctuate predictably in frequency over seasonal time spans" (abstract). | **SUPPORTS**.

**kelleher2018** [doi:10.1371/journal.pcbi.1006581](https://doi.org/10.1371/journal.pcbi.1006581). PLoS Comput Biol 14:e1006581, ok. | p18:804-805 "recapitated with msprime (Kelleher et al., 2018; Haller et al., 2019)". | "The genetic history of the population is recorded as a succinct tree sequence as introduced in the software package msprime" (abstract). It describes forward-time tree-sequence recording, not msprime or recapitation. | **PARTIAL**: cite kelleher2016/baumdicker2022 for msprime (Problem 11).

**kern2018diploshic** [doi:10.1534/g3.118.200262](https://doi.org/10.1534/g3.118.200262). G3 8:1959-1970, ok. | p3:128; p9:392-393; p14:599; p18:803 diploS/HIC on unphased data. | "we describe a deep learning variant of our method, diploS/HIC, that uses unphased genotypes to accurately classify genomic windows" (abstract). | **SUPPORTS**.

**kingston2015** [doi:10.1242/jeb.117945](https://doi.org/10.1242/jeb.117945). J Exp Biol 218:1596-1602, ok. | p3:105-106 dermal opsins in squids; p13:577 phototransduction components in skin and chromatophores. | "RT-PCR revealed the presence of transcripts encoding rhodopsin and retinochrome within the retinas and skin of the squid Doryteuthis pealeii ... Rhodopsin, retinochrome and Gqα transcripts were also found in ... dissociated chromatophores" (abstract). | **SUPPORTS**.

**kirkpatrick2006** [doi:10.1534/genetics.105.047985](https://doi.org/10.1534/genetics.105.047985). Genetics 173:419-434, ok. | p3:86; p7:289-290; p11:481 local adaptation with migration. | "We study the evolution of inversions that capture locally adapted alleles when two populations are exchanging migrants ... By suppressing recombination between the loci, a new inversion can spread" (abstract). | **SUPPORTS**.

**kirkpatrick2010** [doi:10.1371/journal.pbio.1000501](https://doi.org/10.1371/journal.pbio.1000501). PLoS Biol 8:e1000501 ("How and why chromosome inversions evolve"), ok. | p2:78 inversions maintain supergenes; p12:513. | Primer on inversion evolution, including suppression of recombination that protects co-adapted alleles (OA, no abstract; title-level). | **SUPPORTS**.

**klinger2018** [doi:10.1093/bioinformatics/bty361](https://doi.org/10.1093/bioinformatics/bty361). Bioinformatics 34:3591-3593, ok. | p6:233-234; p18:778 ABC-SMC in pyABC. | "we developed pyABC: a distributed and scalable ABC-Sequential Monte Carlo (ABC-SMC) framework" (abstract). | **SUPPORTS**.

**korneliussen2014angsd** [doi:10.1186/s12859-014-0356-4](https://doi.org/10.1186/s12859-014-0356-4). BMC Bioinf 15:356, ok. | p18:760-761 GL-based SFS/θ/Tajima's D from BAMs. | "perform ... population genetic analyses utilizing the full information in next generation sequencing data by working directly on the raw sequencing data or by using genotype likelihoods" (abstract). | **SUPPORTS**.

**lamichhaney2016** [doi:10.1038/ng.3430](https://doi.org/10.1038/ng.3430). Nat Genet 48:84-88, ok. | p3:82-83 ruff morphs. | "The Satellite and Faeder alleles are both associated with a 4.5-Mb inversion" (abstract). | **SUPPORTS**.

**lenormand2005** [doi:10.1371/journal.pbio.0030063](https://doi.org/10.1371/journal.pbio.0030063). PLoS Biol 3:e63, ok. | p16:683 loss of recombination in the heterogametic sex. | "In species in which meiosis is achiasmate in one sex, it is the heterogametic one. This pattern, known as the Haldane-Huxley rule" (abstract). | **SUPPORTS**.

**li2018minimap2** [doi:10.1093/bioinformatics/bty191](https://doi.org/10.1093/bioinformatics/bty191). Bioinformatics 34:3094-3100, ok. | p18:776. | "Minimap2 is a general-purpose alignment program ... assembly contigs or closely related full chromosomes" (abstract). | **SUPPORTS**.

**li2019lostruct** [doi:10.1534/genetics.118.301747](https://doi.org/10.1534/genetics.118.301747). Genetics 211:289-304, ok. | p5:182; p18:773 local PCA. | "We show how to use local PCA to describe this intermediate-scale heterogeneity in patterns of relatedness" (abstract). | **SUPPORTS**.

**li2022flyingsquid** [doi:10.1016/j.gpb.2022.09.009](https://doi.org/10.1016/j.gpb.2022.09.009). GPB 20:1053-1065, ok. | p5:212 Sthenoteuthis oualaniensis chromosome-level genome used as a second outgroup. | "we obtained high-quality genomes of two purpleback flying squid species (Sthenoteuthis oualaniensis and Sthenoteuthis sp.)"; Fig. 1 "Genome-wide Hi-C map of the 46 pseudo-chromosomes of S. oualaniensis" (Results). | **SUPPORTS**.

**liscovitch2017** [doi:10.1016/j.cell.2017.03.025](https://doi.org/10.1016/j.cell.2017.03.025). Cell 169:191-202.e11, ok. | p3:101 RNA editing. | "RNA editing is particularly common in behaviorally sophisticated coleoid cephalopods, with tens of thousands of evolutionarily conserved sites" (abstract). | **SUPPORTS**.

**liu2020stairway** [doi:10.1186/s13059-020-02196-9](https://doi.org/10.1186/s13059-020-02196-9). Genome Biol 21:280, ok. | p18:763. | "We present Stairway Plot 2 ... capability of handling folded SNP frequency spectra" (abstract). | **SUPPORTS**.

**lowry2010** [doi:10.1371/journal.pbio.1000500](https://doi.org/10.1371/journal.pbio.1000500). PLoS Biol 8:e1000500, ok. | p3:83-84 Mimulus life history; p11:482 clines/local adaptation. | "A set of experiments demonstrates the involvement of a chromosomal inversion in the adaptive transition between annual and perennial ecotypes of the yellow monkeyflower" (summary/title). | **SUPPORTS**.

**manichaikul2010** [doi:10.1093/bioinformatics/btq559](https://doi.org/10.1093/bioinformatics/btq559). Bioinformatics 26:2867-2873, ok. | p4:156; p18:771 KING-robust. | Abstract describes relationship inference robust to population structure (KING). | **SUPPORTS**.

**mathger2010** [doi:10.1098/rsbl.2010.0223](https://doi.org/10.1098/rsbl.2010.0223). Biol Lett 6:600-603, ok. | p3:105 "In squids ... dermal opsins". | "We report that the skin of cuttlefish, Sepia officinalis, contains opsin transcripts suggesting a possible role of distributed light sensing" (abstract). | **PARTIAL**: cuttlefish, not squid.

**matuszewski2015** [doi:10.1534/genetics.115.178574](https://doi.org/10.1534/genetics.115.178574). Genetics 200:1255-1274, ok. | p12:505-506 "a moving selective optimum can hold an allele at an intermediate frequency". | "we derive the distribution of adaptive substitutions from standing genetic variation, that is, the distribution of the phenotypic effects of those alleles from the standing variation that become fixed during adaptation" (abstract). The full text has no statement about maintenance at intermediate frequency. | **DOES NOT SUPPORT**.

**mcdonald1991** [doi:10.1038/351652a0](https://doi.org/10.1038/351652a0). Nature 351:652-654, ok. | p9:365 MK tests. | "we propose a simple statistical test of the neutral protein evolution hypothesis based on a comparison of the number of amino-acid replacement substitutions to synonymous substitutions" (abstract). | **SUPPORTS**.

**mcvicker2009** [doi:10.1371/journal.pgen.1000471](https://doi.org/10.1371/journal.pgen.1000471). PLoS Genet 5:e1000471, ok. | p9:360; p14:587 BGS. | "Neutral sequence diversity ... is substantially reduced near such features ... broadly consistent with 'background selection'" (abstract). | **SUPPORTS**.

**merot2020** [doi:10.1016/j.tree.2020.03.002](https://doi.org/10.1016/j.tree.2020.03.002). TREE 35:561-572, ok. | p13:563 increasing numbers of SVs/inversions found. | "Structural genomic variants (SVs) are ubiquitous and play a major role in adaptation and speciation" (abstract). | **SUPPORTS**.

**messer2013mk** (2013a). PNAS 110:8615-8620, ok. | p9:365 MK tests / α. | "we ... investigate whether such linkage effects impinge on the ability of the McDonald-Kreitman (MK) test to infer the rate of positive selection (α)" (abstract). | **SUPPORTS**.

**messer2013soft** (2013b). TREE 28:659-669, ok. | p14:612-613 large populations with abundant standing variation favour soft sweeps. | "soft sweeps are expected ... (i) in very large populations ... or (iii) when adaptation utilizes alleles present as standing genetic variation" (main text). | **SUPPORTS**.

**min2026** [doi:10.1093/genetics/iyag107](https://doi.org/10.1093/genetics/iyag107). Genetics 233(3):iyag107, ok. | p6:243 NPE for population genetics. | "we apply a recently introduced method, neural posterior estimation ... by training a neural network to estimate the posterior distribution of a population genetics model" (abstract). | **SUPPORTS**.

**mkado** [doi:10.1093/g3journal/jkag142](https://doi.org/10.1093/g3journal/jkag142). G3 16(8):jkag142, ok. | p18:818. | "MKado is a Python toolkit for performing McDonald-Kreitman (MK) tests" (abstract). | **SUPPORTS**.

**msinv** [GitHub](https://github.com/stsmall/msinv). GitHub repository and commit verified to exist; **version/commit concern** (Problem 2). | p6:230; p18:779 "a sequentially Markov coalescent simulator that enforces the recombination barrier between arrangements". | Repository stsmall/msinv (an ARG/SMC inversion simulator). | **SUPPORTS** (software exists; confirm the version).

**navarro2000** [doi:10.1093/genetics/155.2.685](https://doi.org/10.1093/genetics/155.2.685). Genetics 155:685-698, ok. | p3:95; p6:219-220 "without gene flux the neutral diversity of each scales with its frequency"; p13:554. | "we develop analytical expressions to predict nucleotide variability in old inversion polymorphisms that have reached mutation-drift-flux equilibrium" (abstract; full text blocked). | **PARTIAL**: abstract-level only, and the frequency-scaling statement could not be confirmed.

**nielsen2005** [doi:10.1101/gr.4252305](https://doi.org/10.1101/gr.4252305). Genome Res 15:1566-1575, ok. | p9:408 composite-likelihood SFS sweep method. | "a new parametric test, based on composite likelihood, has a high power to detect selective sweeps" (abstract). | **SUPPORTS**.

**odor1998** [link](https://openknowledge.fao.org/handle/20.500.14283/w9000e). FAO Fish Tech Pap 376, ch. 5, pp. 77-104; URL resolves. | p2:51, 55 life history and range; p14:582 visually guided predator with diel vertical migration. | Full text not retrieved (FAO repository). | **CANNOT CHECK** (species account in the standard Illex volume, topically appropriate).

**peischl2013** [doi:10.1038/hdy.2013.38](https://doi.org/10.1038/hdy.2013.38). Heredity 111:200-209, ok. | p6:220 diversity scales with frequency. | "we extend the sequential Markovian coalescent ... to include the effects of polymorphic inversions on patterns of recombination" (abstract). | **PARTIAL**: a simulation-method paper; it does not state the analytic claim.

**pggpu** [doi:10.64898/2026.05.29.728868](https://doi.org/10.64898/2026.05.29.728868). bioRxiv 2026, DOI 10.64898/2026.05.29.728868 resolves via Crossref (posted 2026-05-31; authors match). | p18:767. | "We present pg_gpu, a Python library implementing a comprehensi[ve] ..." (abstract). | **SUPPORTS**.

**pockrandt2020genmap** [doi:10.1093/bioinformatics/btaa222](https://doi.org/10.1093/bioinformatics/btaa222). Bioinformatics 36:3687-3692, ok. | p17:755 mappability. | "We present a fast method GenMap to compute the (k, e)-mappability" (abstract). | **SUPPORTS**.

**poplin2018gatk** [doi:10.1101/201178](https://doi.org/10.1101/201178). bioRxiv 201178, year differs (Problem 4). | p17:750 HaplotypeCaller gVCF joint calling. | "the GATK HaplotypeCaller (HC) and Reference Confidence Model (RCM), that determines genotype likelihoods independently per-sample but performs joint calling across all samples" (abstract). | **SUPPORTS**.

**prangle2017** [doi:10.1214/16-BA1002](https://doi.org/10.1214/16-BA1002). Bayesian Anal 12:289-309, ok. | p18:787 adaptive distance refitted each generation. | "Typically weights are chosen at the start of the algorithm ... However these may not be appropriate in iterative ABC algorithms, where the distribution from which the parameters are proposed is updated" (abstract). | **SUPPORTS**.

**ramirez2015** [doi:10.1242/jeb.110908](https://doi.org/10.1242/jeb.110908). J Exp Biol 218:1513-1520, ok. | p3:105 "In squids ..."; p13:577 skin phototransduction. | "light causes chromatophores to expand in excised pieces of Octopus bimaculoides skin ... octopus skin is intrinsically light sensitive" (abstract). | **PARTIAL**: octopus, not squid, at p3:105; fine at p13:577.

**rudman2022** [doi:10.1126/science.abj7484](https://doi.org/10.1126/science.abj7484). Science 375:eabj7484, ok. | p3:92; p12:504 temporally fluctuating selection. | "The direction and basis of the adaptive response shifted repeatedly over time, consistent with the action of strong and rapidly fluctuating selection" (abstract). | **SUPPORTS**.

**schmidbaur2022** [doi:10.1038/s41467-022-29694-7](https://doi.org/10.1038/s41467-022-29694-7). Nat Commun 13:2172, ok. | p2:43-44; p3:100; p13:570 regulatory reorganization around neural genes. | "coleoid cephalopod genomes have been extensively restructured ... we identify a set of microsyntenies associated with cephalopod innovations (MACIs) broadly enriched in cephalopod nervous system expression" (abstract). | **SUPPORTS**.

**schrider2016** [doi:10.1371/journal.pgen.1005928](https://doi.org/10.1371/journal.pgen.1005928). PLoS Genet 12:e1005928, ok. | p9:392-394 five classes (neutral, hard, soft, linked-hard, linked-soft), 11 sub-windows; p14:599. | "S/HIC classifies each window as a hard sweep, linked to a hard sweep, a soft sweep, linked to a soft sweep, or neutral ... examining values of various summary statistics in 11 different windows" (Fig 1 legend). | **SUPPORTS**.

**schwander2014** [doi:10.1016/j.cub.2014.01.056](https://doi.org/10.1016/j.cub.2014.01.056). Curr Biol 24:R288-R294, ok. | p2:78 supergenes. | "'supergenes' segregate as stable polymorphisms within or between natural populations and influence ecologically relevant traits" (abstract). | **SUPPORTS**.

**schweinsberg2003** [doi:10.1016/S0304-4149(03)00028-0](https://doi.org/10.1016/S0304-4149(03)00028-0). Stoch Proc Appl 106:107-139, ok. | p8:342-343 multiple-merger (Beta) coalescent. | Abstract not retrievable. Title: "Coalescent processes obtained from supercritical Galton–Watson processes". | **CANNOT CHECK** (this is the standard source for the Beta(2−α, α) coalescent, so very likely fine).

**small2026** (no link in bib). `@unpublished`, "Manuscript in preparation", no DOI (Problem 3). | p2:44 cephalopod genomes; p3:111-112 "chromosome-level annotated reference genome ... recently made available"; p4:141, p15:665 ZZ/ZO; p8:333-334 "1.06% heterozygosity estimated from k-mers"; p9:380 long reads of the reference and a second individual; p11:475-476 Z dosage compensation, sex-biased lncRNAs, Zfest; p16:680-681 Z conserved across coleoids; p17:747 Illex_F24 assembly 4.14 Gb, 45 chromosomes plus Z. | Not publicly available. | **CANNOT CHECK**.

**song2022anchorwave** [doi:10.1073/pnas.2113075119](https://doi.org/10.1073/pnas.2113075119). PNAS 119:e2113075119, ok. | p18:777 whole-genome alignment. | "we introduce Anchored Wavefront alignment (AnchorWave), which performs ... collinear anchor identification between genomes" (abstract). | **SUPPORTS**.

**songcocasey2022** [doi:10.1016/j.cub.2022.10.015](https://doi.org/10.1016/j.cub.2022.10.015). Curr Biol 32:5031-5044.e4, ok. | p3:102-103; p14:579 optic lobes. | "we present a comprehensive single-cell molecular atlas of the octopus optic lobe, which is the primary visual processing structure in the cephalopod brain" (abstract). | **SUPPORTS**.

**stefansson2005** [doi:10.1038/ng1508](https://doi.org/10.1038/ng1508). Nat Genet 37:129-137, ok. | p12:494 "the human 17q21 inversion, whose increase is more recent". | "two distinct lineages, H1 and H2, that have diverged for as much as 3 million years ... Here we show that the H2 lineage is undergoing positive selection in the Icelandic population" (abstract). | **PARTIAL**: ongoing selection, but no direct statement of a "more recent increase".

**styfhals2022** [doi:10.1038/s41467-022-35198-1](https://doi.org/10.1038/s41467-022-35198-1). Nat Commun 13:7392, ok. | p3:102; p13:569; p14:580. | "Octopuses are mollusks that have evolved intricate neural systems comparable with vertebrates ... we profile the cell diversity of the paralarval Octopus vulgaris brain" (abstract). | **SUPPORTS**.

**tajima1989** (1989a). Genetics 123:585-595, ok. | p8:334 Tajima's D. | "a statistical method for testing the neutral mutation hypothesis is developed" (abstract). | **SUPPORTS**.

**tajima1989b** [doi:10.1093/genetics/123.3.597](https://doi.org/10.1093/genetics/123.3.597). Genetics 123:597-601, ok. | p8:334-335 D ≈ −2 as a signature of expansion. | "in the case where the population size has changed drastically, the number of segregating sites is influenced by the size of the current population more strongly than is the average number of nucleotide differences" (abstract). | **SUPPORTS** (negative D under expansion follows directly from this).

**tejero2020** [doi:10.21105/joss.02505](https://doi.org/10.21105/joss.02505). JOSS 5:2505, ok. | p18:788 sbi. | Title: "sbi: A toolkit for simulation-based inference". | **SUPPORTS**.

**thompson2014** [doi:10.1038/hdy.2014.20](https://doi.org/10.1038/hdy.2014.20). Heredity 113:1-8, ok. | p2:78 supergenes. | "major chromosomal inversions have been shown to be associated with polymorphisms in butterflies, ants and birds, offering a mechanism for localised reduction in recombination" (abstract). | **SUPPORTS**.

**toni2009** [doi:10.1098/rsif.2008.0172](https://doi.org/10.1098/rsif.2008.0172). J R Soc Interface 6:187-202, ok. | p6:233 ABC-SMC model choice. | "we develop ABC SMC as a tool for model selection" (abstract). | **SUPPORTS**.

**vasimuddin2019bwamem2** [doi:10.1109/IPDPS.2019.00041](https://doi.org/10.1109/IPDPS.2019.00041). IPDPS 2019 pp. 314-324, ok. | p17:749. | BWA-MEM2 proceedings paper (Crossref metadata match). | **SUPPORTS**.

**wang2013** [doi:10.1038/nature11832](https://doi.org/10.1038/nature11832). Nature 493:664-668, ok. | p3:83 fire ant social organization. | "a large region of approximately 13 megabases ... in which recombination is completely suppressed between SB and Sb" (abstract). | **SUPPORTS**.

**waples1998** [doi:10.1093/jhered/89.5.438](https://doi.org/10.1093/jhered/89.5.438). J Hered 89:438-450, ok. | p2:40 marine datasets. | No abstract retrievable. Title: "Separating the wheat from the chaff: patterns of genetic differentiation in high gene flow species". | **PARTIAL** (title-level; about weak differentiation in high-gene-flow species, not specifically about mtDNA/few loci).

**wellenreuther2018** [doi:10.1016/j.tree.2018.04.002](https://doi.org/10.1016/j.tree.2018.04.002). TREE 33:427-440, ok. | p2:79; p13:562. | "inversions are taxonomically widespread, with many being old and large, and that balancing selection is commonly facilitating their maintenance" (abstract). | **SUPPORTS**.

**wittmann2017** [doi:10.1073/pnas.1702994114](https://doi.org/10.1073/pnas.1702994114). PNAS 114(46):E9932-E9941, ok. | p3:91; p12:504 fluctuating selection maintains polymorphism. | "we show that stable polymorphism at many loci is possible if currently favored alleles are sufficiently dominant ... 'segregation lift'" (abstract). | **SUPPORTS**.

**yoshida2015** [doi:10.1093/icb/icv049](https://doi.org/10.1093/icb/icv049). Integr Comp Biol 55:1070-1083, ok. | p14:578-580 optic lobes evolved independently of the vertebrate visual system. | "Coleoid cephalopods show remarkable evolutionary convergence with vertebrates in their neural organization, including (1) eyes and visual system with optic lobes" (abstract). | **SUPPORTS**.

**zhao2016** [doi:10.1534/genetics.116.188912](https://doi.org/10.1534/genetics.116.188912). Genetics 203:1315-1334, ok. | p3:90; p12:499-500 AOD from recessive deleterious load. | "Apparent heterozygote advantage at the neutral locus is, however, always observed when mutations are partially recessive" (abstract). | **SUPPORTS**.

---
Method notes: keys were extracted with `\cite[a-zA-Z]*` from main.tex. Supplement and table keys were checked only to identify uncited bib entries. The scripts and cached abstracts/full texts are in the scratchpad (scripts/, ft/, epmc.json, crossref_raw.json).


---

# Part 2. Supplement-only citations


Scope: 21 citation keys used in supplement.tex + software_table.tex (\input by supplement.tex line 472) that are not cited in main.tex.
Page numbers are PDF pages of supplement.pdf (they match the printed folio). Metadata checked against Crossref (all 19 DOIs resolve), Europe PMC, arXiv and PMLR.

## Supplement: Summary

Ratings (21 keys):
- SUPPORTS: 16 (baril2024earlgrey, baumdicker2022, chao2024lifton, danecek2021samtools, eddy2011hmmer, gower2022demes, kern2016discoal, kovaka2019stringtie, lemos2023, mcvean2005smc, mistry2021pfam, pardo2024sqanti3, quinlan2010bedtools, ragsdale2020unphased, ranwez2018macse, talts2018)
- PARTIAL: 3 (korunes2019, kelleher2016, miyata1987)
- DOES NOT SUPPORT: 0
- CANNOT CHECK (verbatim): 2 (nei1986, which is supported at abstract level only; ohta1971, where no full text was reachable)

No work is missing and no DOI is wrong. Problems:

1. **chao2024lifton, incomplete or outdated metadata.** The bib has no volume or pages. The issue of record is *Genome Research* 35(2):311–325, **2025** (online 27 Dec 2024; PMID 39730188, PMC11874971). Add the volume and pages, and consider changing the year to 2025.
2. **talts2018, no identifier beyond the arXiv note.** Add the arXiv DOI 10.48550/arXiv.1804.06788 (or eprint fields). The citation itself is correct.
3. **lemos2023, incomplete proceedings data.** Add PMLR vol. 202, pp. 19256–19273 (ICML 2023). The title, authors and year are correct.
4. **msprime version inconsistency.** The text on p. 4 ('Inversion coalescent model') says "msprime 1.4.1". Table S1 (p. 33) says "msprime 1.4.2". Both versions exist (1.4.1 released 2026-03-05, 1.4.2 released 2026-05-20). Reconcile them, or say that different steps used different versions.
5. **msinv version inconsistency, in the kelleher2016 row.** Table S1 lists "msinv 0.3.5". The text on p. 3 says "msinv v0.4.1", and the memory notes say all pre-0.4.1 inversion simulations are stale. The table row is out of date.
6. **StringTie name and version mismatch.** The text on p. 3 says "transcripts from StringTie2 (Kovaka et al., 2019)". Table S1 says "StringTie 3.0.0" (released 2025-01-07) with the same citation. Kovaka et al. 2019 is still the citation the developers request for long-read mode, so it is acceptable. Make the text and table agree, e.g. "StringTie (v3.0.0; Kovaka et al., 2019)".
7. **korunes2019, the 300 bp value is not stated by the source.** Korunes & Noor report mean tract lengths of 188–390 bp. 300 bp sits within that range but does not appear in the paper. Suggested wording: "tracts of mean length 300 bp, within the 188–390 bp estimated by Korunes and Noor (2019)".
8. **miyata1987, the formula is adapted rather than taken from the paper.** Miyata et al. derived male-driven-evolution ratios for X, Y and autosomes (mammalian XY). The Z/A expression α=(3R−2)/(4−3R) is the ZW analogue of their model. I checked the algebra: Z/A = 2(2α+1)/[3(α+1)] gives α=(3R−2)/(4−3R), which is correct. Attributing it to Miyata alone is conventional but loose. Consider "following the model of Miyata et al. (1987)" or adding a ZW source (e.g. Ellegren & Fridolfsson 1997; Kirkpatrick & Hall 2004). No abstract or full text was reachable, so the wording of the original could not be checked. Same issue in the Table S12/zsex caption (p. 42).
9. **kelleher2016 cited for msinv (Table S1).** This is the original msprime / sparse-tree paper. It supports "tree-sequence output" only loosely. The canonical tskit citation (kelleher2018) is already in the table on the tskit row. Consider citing kelleher2018 for msinv's tree-sequence output, or drop kelleher2016.
10. **Table S1 caption has draft text left in.** The caption ends "Citations are pointers to confirm at final formatting." Remove it before submission. Also, a source comment in software_table.tex says it is "Included by main.tex", but it is actually \input by supplement.tex. This is cosmetic.
11. **Software rows with no version** (minor, for G3 reproducibility): SQANTI3, LiftOn and Pfam-A (the release, e.g. Pfam 37.x). The 2021 Pfam paper describes Pfam 33.1.
12. **baumdicker2022 author spelling.** The bib has "Kretzschumar". This matches the publisher/Crossref record, which contains the typo; the person is Warren W. Kretzschmar. Optional fix.
13. **ohta1971.** The bib title is complete and correct. No verbatim check was possible (the PMC scan was blocked to scripted access).

Versions checked as plausible against GitHub releases: EarlGrey 7.0.1 (2026-01-06), StringTie 3.0.0, discoal 0.1.7a (tag exists), msprime 1.4.1/1.4.2, demes 0.2.3 (2023-05), bedtools 2.30.0 (2021-01), HMMER 3.4 (2023-08), bcftools 1.21 (2024-09), MACSE 2.07 (the v2 paper covers 2.x).

---

## Supplement: Per-key blocks

### baril2024earlgrey [doi:10.1093/molbev/msae068](https://doi.org/10.1093/molbev/msae068)
- **Metadata:** OK. Baril, Galbraith, Hayward. *Mol Biol Evol* 41(4):msae068, 2024. DOI 10.1093/molbev/msae068. Matches Crossref.
- **Supplement location / claim:** p. 3, para 'The inversion', last sentences: "Unannotated, non-repeat sequence in the inversion was searched with DIAMOND blastx against seven cephalopod proteomes, with transposable elements from Earl Grey (Baril et al., 2024) and transcripts from StringTie2 ... as supporting evidence." Also Table S1, p. 33: "EarlGrey 7.0.1 – TE annotation".
- **Quote (abstract, PMC11003543):** "we present Earl Grey, a fully automated TE annotation pipeline designed for user-friendly curation and annotation of TEs in eukaryotic genome assemblies."
- **Rating:** SUPPORTS. It is the canonical citation for the tool, and v7.0.1 is a real release.

### baumdicker2022 [doi:10.1093/genetics/iyab229](https://doi.org/10.1093/genetics/iyab229)
- **Metadata:** OK. *Genetics* 220(3):iyab229, 2022, with 32 authors that match Crossref. DOI correct. "Kretzschumar" matches the publisher record but the real spelling is Kretzschmar (optional fix).
- **Supplement location / claim:** p. 4, para 'Inversion coalescent model', sentence 9: "...a separate 40 kb of collinear sequence with msprime 1.4.1 (Baumdicker et al., 2022) as the control..." Table S1, p. 33: "msprime 1.4.2 – coalescent sim, mutation overlay".
- **Quote (abstract, PMC9176297):** "Here, we introduce msprime version 1.0, which efficiently implements ancestry and mutation simulations based on the succinct tree sequence data structure and the tskit library."
- **Rating:** SUPPORTS. It is the canonical citation. **Version conflict: 1.4.1 in the text vs 1.4.2 in the table.**

### chao2024lifton [doi:10.1101/gr.279620.124](https://doi.org/10.1101/gr.279620.124)
- **Metadata:** The DOI 10.1101/gr.279620.124 is correct, and the authors and title match. **Missing volume and pages: *Genome Res* 35(2):311–325 (2025).** The year in the bib is 2024 (the online date); the issue year is 2025.
- **Supplement location / claim:** (a) p. 2–3, para 'The inversion': "...HiC_scaffold_1, is the scaffold that receives 265 of the 319 Illex chromosome-2 genes transferred with LiftOn (Chao et al., 2024)". (b) p. 7, para 'Z recombination and sex-biased mutation', para 2, sentence 2: "The Illex illecebrosus annotation was transferred to the S. oualaniensis assembly with LiftOn (Chao et al., 2024)". (c) Table S1.
- **Quote (abstract, PMC11874971):** "Here we introduce LiftOn, a homology-based software tool that integrates DNA and protein alignments to enhance the accuracy of genome-scale annotation and to allow mapping between relatively distant species."
- **Rating:** SUPPORTS. It is the canonical citation, and the cross-genus transfer use fits the stated scope. No version is given in Table S1.

### danecek2021samtools [doi:10.1093/gigascience/giab008](https://doi.org/10.1093/gigascience/giab008)
- **Metadata:** OK. *GigaScience* 10(2):giab008, 2021, 11 authors. DOI correct.
- **Supplement location / claim:** p. 1, para 'Read processing and variant calling', sentence 4: "Variant records were annotated with snpEff and soft-filtered with bcftools (Danecek et al., 2021)." Table S1: bcftools/samtools 1.21.
- **Quote (abstract, PMC7931819):** "SAMtools and BCFtools are widely used programs for processing and analysing high-throughput sequencing data. They include tools for file format conversion and manipulation, sorting, querying, statistics, variant calling, and effect analysis amongst other methods."
- **Rating:** SUPPORTS. It is the canonical citation, and 1.21 is a real release (2024-09).

### eddy2011hmmer [doi:10.1371/journal.pcbi.1002195](https://doi.org/10.1371/journal.pcbi.1002195)
- **Metadata:** OK. *PLoS Comput Biol* 7(10):e1002195, 2011. DOI correct.
- **Supplement location / claim:** Table S1, p. 33: "HMMER 3.4 – Pfam domain screen".
- **Quote (abstract, PMC3197634):** "This accelerated pipeline is implemented in the freely available HMMER3 software package."
- **Rating:** SUPPORTS. It is the standard HMMER3 citation, and 3.4 is a real release (2023-08).

### gower2022demes [doi:10.1093/genetics/iyac131](https://doi.org/10.1093/genetics/iyac131)
- **Metadata:** OK. *Genetics* 222(3):iyac131, 2022, 10 authors. DOI correct.
- **Supplement location / claim:** Table S1, p. 33: "demes 0.2.3 – demography interchange".
- **Quote (abstract, PMC9630982):** "We propose the Demes data model and file format, built on widely used technologies, to alleviate these issues."
- **Rating:** SUPPORTS. It is the canonical citation, and demes-python 0.2.3 is a real release (2023-05).

### kelleher2016 [doi:10.1371/journal.pcbi.1004842](https://doi.org/10.1371/journal.pcbi.1004842)
- **Metadata:** OK. *PLoS Comput Biol* 12(5):e1004842, 2016. DOI correct.
- **Supplement location / claim:** Table S1, p. 33, msinv row: "msinv 0.3.5 – chr2 inversion age/decline – (Small, 2025; Kelleher et al., 2016)". Presumably this credits the tree-sequence format that msinv writes (p. 3: "a coalescent simulator with recombination that writes tree sequences").
- **Quote (abstract, PMC4856371):** "We solve these problems by introducing sparse trees and coalescence records as the key units of genealogical analysis."
- **Rating:** PARTIAL. The paper introduced the record structure that underlies tree sequences, but the citation commonly used for tree-sequence I/O is tskit (kelleher2018), which the table already cites. **The msinv version in the same row (0.3.5) conflicts with the v0.4.1 in the text.**

### kern2016discoal [doi:10.1093/bioinformatics/btw556](https://doi.org/10.1093/bioinformatics/btw556)
- **Metadata:** OK. *Bioinformatics* 32(24):3839–3841, 2016. DOI correct. (Crossref capitalises the title as "Discoal"; the lower-case form in the bib is fine.)
- **Supplement location / claim:** Table S1, p. 33: "discoal 0.1.7a – sweep-sim validation".
- **Quote (abstract, PMC5167068):** "Here we describe discoal, a coalescent simulator able to generate population samples that include selective sweeps in a feature-rich, flexible manner."
- **Rating:** SUPPORTS. It is the canonical citation, and the v0.1.7a tag exists in kr-colab/discoal.

### korunes2019 [doi:10.1111/mec.14921](https://doi.org/10.1111/mec.14921)
- **Metadata:** OK. *Mol Ecol* 28(6):1302–1315, 2019. DOI correct.
- **Supplement location / claim:** p. 3, para 'Inversion coalescent model', sentence 6: "Lineages of different classes cannot coalesce while the barrier exists, except through gene flux, which copies tracts of mean length 300 bp between arrangements (Korunes and Noor, 2019)..."
- **Quote (full text, PMC6475484, Results 'Gene Conversion Rate and Tract Length' and Discussion):** "The maximum-likelihood approach for the hemizygous × chromosome data estimated a tract length of about 390 bp ... In these data, the maximum likelihood estimation of the tract length was 188 bp" and "Our estimates of average tract length ranged from 188–390 bp, which again falls within the range observed in other studies." Abstract: "gene conversion occurs regularly within and around inverted regions of species hybrids, even near inversion breakpoints."
- **Rating:** PARTIAL. The source supports gene conversion between arrangements and a tract length of the right order (188–390 bp). The specific value of 300 bp is the authors' choice and is not stated in the paper, so reword the sentence to say so.

### kovaka2019stringtie [doi:10.1186/s13059-019-1910-1](https://doi.org/10.1186/s13059-019-1910-1)
- **Metadata:** OK. *Genome Biol* 20:278, 2019. DOI correct.
- **Supplement location / claim:** p. 3, para 'The inversion' (same sentence as Earl Grey): "...transcripts from StringTie2 (Kovaka et al., 2019) as supporting evidence." Table S1: "StringTie 3.0.0 – transcript assembly".
- **Quote (abstract, PMC6912988):** "Here we present StringTie2, a reference-guided transcriptome assembler that works with both short and long reads."
- **Rating:** SUPPORTS. **The text says "StringTie2" but the table says 3.0.0.** Make them consistent. The citation is still the right one for long-read StringTie.

### lemos2023 [arXiv:2302.03026](https://arxiv.org/abs/2302.03026)
- **Metadata:** The work exists: arXiv:2302.03026, published in ICML 2023 as PMLR 202:19256–19273. The authors and title are correct. **Add the PMLR volume and pages (no DOI exists).**
- **Supplement location / claim:** p. 6, para 'Neural posterior estimation' (NPE paragraph), sentence 6: "Calibration was tested on 2,000 held-out simulations by simulation-based calibration (Talts et al., 2018) and by the TARP coverage test (Lemos et al., 2023)."
- **Quote (abstract, arXiv):** "In this paper, we introduce `Tests of Accuracy with Random Points' (TARP) coverage testing as a method to estimate coverage probabilities of generative posterior estimators." From PMLR full text, Sec. 2.3, Theorem 3: "Suppose the expected coverage probability of p̂ is equal to 1 − α for a positionable credible region generator ... for all α ... Then p̂(·|x) = p(·|x)."
- **Rating:** SUPPORTS.

### mcvean2005smc [doi:10.1098/rstb.2005.1673](https://doi.org/10.1098/rstb.2005.1673)
- **Metadata:** OK. *Phil Trans R Soc B* 360(1459):1387–1393, 2005. DOI correct.
- **Supplement location / claim:** p. 3, para 'Inversion coalescent model', sentence 2: "It uses the sequentially Markov approximation (McVean and Cardin, 2005), ... in which two lineages coalesce only where they share ancestral material."
- **Quote (abstract, PMC1569517):** "We introduce a simplification of the coalescent process in which coalescence between lineages with no overlapping ancestral material is banned. The resulting process has a simple Markovian structure when generating genealogies sequentially along a sequence..."
- **Rating:** SUPPORTS. The extension to two karyotype classes is credited separately to Peischl et al. (2013).

### mistry2021pfam [doi:10.1093/nar/gkaa913](https://doi.org/10.1093/nar/gkaa913)
- **Metadata:** OK. *Nucleic Acids Res* 49(D1):D412–D419, 2021, 12 authors. DOI correct.
- **Supplement location / claim:** Table S1, p. 33: "Pfam-A – domain screen".
- **Quote (abstract, PMC7779014):** "The Pfam database is a widely used resource for classifying protein sequences into families and domains."
- **Rating:** SUPPORTS. The table gives no Pfam release; the paper describes Pfam 33.1.

### miyata1987 [doi:10.1101/SQB.1987.052.01.094](https://doi.org/10.1101/SQB.1987.052.01.094)
- **Metadata:** OK. *Cold Spring Harb Symp Quant Biol* 52:863–867, 1987. DOI correct (PMID 3454295).
- **Supplement location / claim:** p. 7, para 'Z recombination and sex-biased mutation', para 2: "The male-to-female mutation ratio α=μm/μf was obtained from the ratio R of mean Z to mean autosomal dS as α=(3R−2)/(4−3R) (Miyata et al., 1987), as in Coffing et al. (2025) and Coffing et al. (2026)." Also the Table S12 (zsex) caption, p. 42.
- **Quote:** None available. There is no abstract in PubMed or Europe PMC and the full text is not open access.
- **Rating:** PARTIAL. Miyata et al. is the original male-driven-evolution model with X/Y/autosome ratios. The Z/A form is its ZW analogue, and I checked that the algebra is correct. Suggest "following Miyata et al. (1987)" and/or adding a ZW-specific source. The paper's wording could not be checked.

### nei1986 [doi:10.1093/oxfordjournals.molbev.a040410](https://doi.org/10.1093/oxfordjournals.molbev.a040410)
- **Metadata:** OK. *Mol Biol Evol* 3(5):418–426, 1986. DOI correct. Crossref lacks the volume and pages but the publisher record agrees with the bib.
- **Supplement location / claim:** p. 7, para 'Z recombination and sex-biased mutation', para 2, sentence 3: "Synonymous and nonsynonymous divergence (dS, dN) were computed per gene by the method of Nei and Gojobori (1986) with a Jukes–Cantor correction."
- **Quote (abstract only, Europe PMC):** "Two simple methods for estimating the numbers of synonymous and nonsynonymous nucleotide substitutions are presented."
- **Rating:** CANNOT CHECK (verbatim, beyond the abstract). The paper is the canonical source for the method, and to my knowledge it applies the Jukes–Cantor correction to pS/pN, but I could not reach the full text to quote that.

### ohta1971 [doi:10.1093/genetics/68.4.571](https://doi.org/10.1093/genetics/68.4.571)
- **Metadata:** OK. *Genetics* 68(4):571–580, 1971. DOI correct. The full title matches Crossref.
- **Supplement location / claim:** p. 7, para 'Z recombination and sex-biased mutation', sentence 3: "The decay of σD²=D²/π2 with distance was fitted to the expectation of Ohta and Kimura (1971) plus a constant sample-size term, giving ρ=4Ner..."
- **Quote:** None. The PMC1212677 scanned PDF was blocked to scripted access and there is no abstract.
- **Rating:** CANNOT CHECK. To my knowledge this is the canonical source of σd² ≈ (10+ρ)/(22+13ρ+ρ²), the standard citation for this fit, but I did not verify it here.

### pardo2024sqanti3 [doi:10.1038/s41592-024-02229-2](https://doi.org/10.1038/s41592-024-02229-2)
- **Metadata:** OK. *Nat Methods* 21(5):793–797, 2024, 12 authors. DOI correct.
- **Supplement location / claim:** Table S1, p. 33: "SQANTI3 – Iso-Seq QC" (no version).
- **Quote (abstract, PMC11093726):** "SQANTI3 is a tool designed for the quality control, curation and annotation of long-read transcript models obtained with third-generation sequencing technologies."
- **Rating:** SUPPORTS. No version is given.

### quinlan2010bedtools [doi:10.1093/bioinformatics/btq033](https://doi.org/10.1093/bioinformatics/btq033)
- **Metadata:** OK. *Bioinformatics* 26(6):841–842, 2010. DOI correct.
- **Supplement location / claim:** Table S1, p. 33: "bedtools 2.30.0 – interval operations".
- **Quote (abstract, PMC2832824):** "This article introduces a new software suite for the comparison, manipulation and annotation of genomic features in Browser Extensible Data (BED) and General Feature Format (GFF) format."
- **Rating:** SUPPORTS. 2.30.0 is a real release (2021-01).

### ragsdale2020unphased [doi:10.1093/molbev/msz265](https://doi.org/10.1093/molbev/msz265)
- **Metadata:** OK. *Mol Biol Evol* 37(3):923–932, 2020. DOI correct.
- **Supplement location / claim:** p. 7, para 'Z recombination and sex-biased mutation', sentence 2: "Two-locus statistics (D² and π2, unbiased for unphased genotypes (Ragsdale and Gravel, 2020)) were computed with pg_gpu..."
- **Quote (abstract and Results, PMC7038669):** "Here, we show how to find unbiased estimators for a wide range of two-locus statistics, including D2, for both single and multiple randomly mating populations." Results, 'Estimating LD from Data': "we present an approach to compute unbiased estimators for a large family of two-locus statistics, using either phased or unphased data."
- **Rating:** SUPPORTS.

### ranwez2018macse [doi:10.1093/molbev/msy159](https://doi.org/10.1093/molbev/msy159)
- **Metadata:** OK. *Mol Biol Evol* 35(10):2582–2584, 2018. DOI correct.
- **Supplement location / claim:** p. 7, para 'Z recombination and sex-biased mutation', para 2, sentence 2: "...one transcript per gene was aligned codon by codon with MACSE 2.07 (Ranwez et al., 2018), which allows frameshifts in the transferred sequence..." Also Table S1.
- **Quote (abstract, PMC6188553):** "Multiple Alignment of Coding Sequences (MACSE) is a multiple sequence alignment program that explicitly accounts for the underlying codon structure of protein-coding nucleotide sequences. Its unique characteristic allows building reliable codon alignments even in the presence of frameshifts."
- **Rating:** SUPPORTS. Version 2.07 falls under the v2 paper.

### talts2018 [arXiv:1804.06788](https://arxiv.org/abs/1804.06788)
- **Metadata:** The work exists: arXiv:1804.06788, with five authors and a title that match. The bib is @misc with an arXiv note only. **Add the DOI 10.48550/arXiv.1804.06788 or eprint fields.** The paper was never published in a journal, so arXiv is the correct record.
- **Supplement location / claim:** p. 6, NPE paragraph, sentence 6 (same as lemos2023): calibration by simulation-based calibration. The following sentence reports the rank distributions are uniform (KS p 0.12–0.93).
- **Quote (arXiv v2, Theorem 1, SBC section):** "Let θ̃ ∼ π(θ), ỹ ∼ π(y | θ̃), and {θ1, ..., θL} ∼ π(θ | ỹ) for any joint distribution π(y, θ). The rank statistic of any one-dimensional random variable over θ is uniformly distributed over the integers [0, L]." Abstract: "In this paper we introduce simulation-based calibration (SBC), a general procedure for validating inferences from Bayesian algorithms capable of generating posterior samples."
- **Rating:** SUPPORTS.
