#!/usr/bin/env python
# pg_gpu diversity on the VARIANT VCF, span-normalized by the accessibility mask
# (callable-site denominator) with per-site missing-data handling. Demonstrates
# that pi/theta_w/Tajima's D are correctly estimated from variants-only + mask —
# no invariant sites needed. chr24 test.
import os
os.environ.setdefault("CUDA_VISIBLE_DEVICES", "0")
import numpy as np
from pg_gpu import HaplotypeMatrix, diversity

BASE = "/sietch_colab/data_share/illex/popgen_data"
VCF = f"{BASE}/analysis/steps/00_callset/filtered/24/variants_filt.vcf.gz"
ACCESS = f"{BASE}/degenotate_illex/accessible_sites.bed"
POPS = f"{BASE}/pggpu_illex/popstats/pops.tsv"

print(f"[load] {VCF}", flush=True)
h = HaplotypeMatrix.from_vcf(VCF, accessible_bed=ACCESS)
print(f"  n_variants   = {h.num_variants}", flush=True)
print(f"  n_haplotypes = {h.num_haplotypes}  ({h.num_haplotypes//2} diploids)", flush=True)
print(f"  n_total_sites (accessible denom, chr24) = {h.n_total_sites}", flush=True)

ms = h.summarize_missing_data()
print(f"  missing-data summary: {ms}", flush=True)

# span-normalized (per callable base) with per-site n_valid (missing_data='include')
pi = float(diversity.pi(h, missing_data="include"))
tw = float(diversity.theta_w(h, missing_data="include"))
td = float(diversity.tajimas_d(h, missing_data="include"))
print("\n=== chr24 diversity (span-normalized per callable base, missing=include) ===", flush=True)
print(f"  pi        = {pi:.6e}", flush=True)
print(f"  theta_w   = {tw:.6e}", flush=True)
print(f"  Tajimas D = {td:.4f}", flush=True)

# per-NAFO-division (from pops.tsv)
try:
    h.load_pop_file(POPS)
    print("\n=== per-NAFO-division pi / Tajima's D ===", flush=True)
    for pop in sorted(h.sample_sets):
        n = len(h.sample_sets[pop])
        p = float(diversity.pi(h, population=pop, missing_data="include"))
        d = float(diversity.tajimas_d(h, population=pop, missing_data="include"))
        print(f"  {pop:>4}: n_hap={n:<5} pi={p:.6e}  TajD={d:+.4f}", flush=True)
except Exception as e:
    print(f"[warn] per-pop step failed: {e}", flush=True)

print("\nDONE", flush=True)
