#!/usr/bin/env python3
"""
Collapse a DIPLOID-coded hemizygous VCF (female Z: ZW/ZO -> one Z copy) to TRUE HAPLOID GT,
streaming stdin -> stdout. FORMAT must be GT-only (verified for Z_female.kept.vcf).

Mapping per sample genotype (alleles split on '/' or '|'):
  homozygous a/a (a in 0,1,..) -> a          (the single real Z allele)
  heterozygous a/b (a != b)    -> .          (artifact on a hemizygous chromosome; drop)
  any missing ('.' present)    -> .

scikit-allel's vcf_to_hdf5 reads single-allele GT as ploidy-2 padded with -1 in the 2nd slot,
so ReLERNN's ploidy detector (2nd allele all-missing => haploid) sets nSamps = n_samples (280),
NOT 2*n_samples. Verified empirically before use.
"""
import sys

def hap(gt):
    # gt is the sample's GT string (FORMAT is GT-only)
    a = gt.replace("|", "/").split("/")
    if len(a) == 1:
        # already haploid
        return a[0] if a[0] != "." else "."
    if "." in a:
        return "."
    if a[0] == a[1]:
        return a[0]
    return "."   # heterozygous -> missing (artifact on hemizygous Z)

for line in sys.stdin:
    if line.startswith("#"):
        sys.stdout.write(line)
        continue
    f = line.rstrip("\n").split("\t")
    fmt = f[8]
    if fmt != "GT":
        sys.exit(f"ERROR: FORMAT is '{fmt}', expected GT-only. Aborting to avoid corrupting genotypes.")
    f[9:] = [hap(g) for g in f[9:]]
    sys.stdout.write("\t".join(f) + "\n")
