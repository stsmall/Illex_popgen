"""Generate K subsample lists of NDIP samples each, without replacement, seeded."""
import sys, numpy as np, os

WD = "/sietch_colab/data_share/illex/popgen_data/analysis/steps/13_diploshic"
K = int(sys.argv[1])
ND = int(sys.argv[2])
samples = [s.strip() for s in open(f"{WD}/all_samples.txt")]
os.makedirs(f"{WD}/subsamples", exist_ok=True)
rng = np.random.default_rng(12345)
for k in range(K):
    sub = rng.choice(samples, size=ND, replace=False)  # WITHOUT replacement
    open(f"{WD}/subsamples/sub_{k:02d}.txt", "w").write("\n".join(sub) + "\n")
print(f"wrote {K} subsample lists of {ND} from {len(samples)} total samples")
