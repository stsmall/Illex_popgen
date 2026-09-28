"""
80/20 train/test split of makeTrainingSets output.
Reads each .fvec from fvec/train, shuffles rows (seeded),
keeps 80% in train and writes 20% to test.
"""
import os, glob, random

WD = "/sietch_colab/data_share/illex/popgen_data/analysis/steps/13_diploshic"
os.makedirs(f"{WD}/fvec/test", exist_ok=True)

for f in sorted(glob.glob(f"{WD}/fvec/train/*.fvec")):
    lines = open(f).read().splitlines()
    hdr = lines[0]
    rows = lines[1:]
    # Remove trailing empty lines
    rows = [r for r in rows if r.strip()]
    random.Random(0).shuffle(rows)
    k = max(1, len(rows) // 5)   # 20% for test
    train_rows = rows[k:]
    test_rows = rows[:k]
    # Overwrite train file with 80%
    open(f, "w").write("\n".join([hdr] + train_rows) + "\n")
    # Write test file
    test_path = f.replace("/train/", "/test/")
    open(test_path, "w").write("\n".join([hdr] + test_rows) + "\n")
    print(f"  {os.path.basename(f)}: {len(train_rows)} train, {len(test_rows)} test")

print("Train/test split complete.")
