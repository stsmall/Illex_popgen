"""Per-arrangement pi_AA, pi_BB, dxy, fst in 100-kb windows across chr2:40-100 Mb (inversion + 20 Mb
collinear flanks), replicating msinv empirical_jackknife.windowed() on the full chr2 filtered callset."""
import os; os.environ.setdefault("CUDA_VISIBLE_DEVICES", "0")
import pandas as pd
from pg_gpu import HaplotypeMatrix, windowed_analysis
VCF = "/sietch_colab/data_share/illex/popgen_data/analysis/steps/00_callset/filtered/2/variants_filt.vcf.gz"
POPS = "/sietch_colab/ssmall/projects/msinv_dir/inversion_sims/files/.tmp/illex_chr2/pops.tsv"
OUT = "chr2_40_100_karyo_windows.csv"
h = HaplotypeMatrix.from_vcf(VCF, region="2:40000000-100000000")
h.load_pop_file(POPS)
print(f"loaded {h.num_variants:,} variants x {h.num_haplotypes:,} haplotypes (AA {len(h.sample_sets['AA'])}, BB {len(h.sample_sets['BB'])}) device={h.device}", flush=True)
kw = dict(window_size=100_000, step_size=100_000, missing_data="include")
df = windowed_analysis(h, statistics=["fst", "dxy"], populations=["AA", "BB"], **kw)
for pop in ("AA", "BB"):   # pi must be requested per population (pg_gpu bug)
    p = windowed_analysis(h, statistics=["pi"], populations=[pop], **kw)
    df = df.merge(p[["window_id", "pi"]].rename(columns={"pi": f"pi_{pop}"}), on="window_id", how="left")
df = df.rename(columns={"start": "window_start", "end": "window_stop"})
df.to_csv(OUT, index=False); print("wrote", OUT, df.shape, flush=True)
