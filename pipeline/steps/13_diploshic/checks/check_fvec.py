import os, glob
WD = "/sietch_colab/data_share/illex/popgen_data/analysis/steps/13_diploshic"
need = ["hard", "soft", "linkedHard", "linkedSoft", "neutral"]
for split in ("train", "test"):
    for n in need:
        p = f"{WD}/fvec/{split}/{n}.fvec"
        assert os.path.exists(p) and os.path.getsize(p) > 0, f"missing/empty: {p}"
        assert len(open(p).readline().split()) > 50, f"too few columns: {p}"  # 11 subwin x several stats
subdirs = glob.glob(f"{WD}/fvec/real/sub*")
assert len(subdirs) >= 50, f"only {len(subdirs)} subdirs"
for sd in subdirs[:3]:
    nfvec = len(glob.glob(f"{sd}/*.fvec"))
    assert nfvec >= 40, f"only {nfvec} fvecs in {sd}"
assert len(glob.glob(f"{WD}/subsamples/sub_*.txt")) == 50
print("fvec check PASSED")
