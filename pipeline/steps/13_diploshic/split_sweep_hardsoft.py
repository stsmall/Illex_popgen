"""
Classify each sweep sim as hard vs soft from its TRUE label — the raw discoal
command line — then write concatenated per-(type,subwin) ascertained ms files.

Usage: python split_sweep_hardsoft.py
Output:
  WD/fvec/msc/hard_${s}.msc  -- concatenated hard sweep sims for subwin s
  WD/fvec/msc/soft_${s}.msc  -- concatenated soft sweep sims for subwin s

NOTE: an earlier version re-derived hard/soft by replaying
hash(("sweep",sub,rep)) to reconstruct the RNG. That is broken: Python string
hashing is per-process randomized (PYTHONHASHSEED unset), so the replay seed
never matches the seed used at generation time -> labels agreed with the truth
only ~50% of the time (chance), scrambling the hard/soft training classes.
"""
import sys, os, glob

WD = "/sietch_colab/data_share/illex/popgen_data/analysis/steps/13_diploshic"
os.makedirs(f"{WD}/fvec/msc", exist_ok=True)

# Parse rep number from filename like: sweep_s${sub}_r${rep}.ms
def parse_sub_rep(fname):
    base = os.path.basename(fname)
    # e.g. sweep_s3_r17.ms
    parts = base.replace(".ms", "").split("_")
    sub = int(parts[1][1:])   # s3 -> 3
    rep = int(parts[2][1:])   # r17 -> 17
    return sub, rep

def is_soft(asc_path):
    """True label from the raw discoal command line: soft sweeps carry ' -f '.

    Ascertainment strips the discoal header (asc files start with '//'), so we
    read the corresponding RAW file, which preserves the full command line.
    """
    raw_path = asc_path.replace("/sims/asc/", "/sims/raw/")
    if not os.path.exists(raw_path):
        raise FileNotFoundError(
            f"cannot label {asc_path}: raw file {raw_path} missing (need discoal cmd line)"
        )
    with open(raw_path) as fh:
        cmd = fh.readline()
    return " -f " in cmd

# Gather all asc sweep files
all_sweep = glob.glob(f"{WD}/sims/asc/sweep/sweep_s*_r*.ms")
print(f"Found {len(all_sweep)} asc sweep files")

# Organize by (type, subwin) using the TRUE label from each raw command line
hard_files = {s: [] for s in range(11)}
soft_files = {s: [] for s in range(11)}

for f in sorted(all_sweep):
    sub, rep = parse_sub_rep(f)
    if is_soft(f):
        soft_files[sub].append(f)
    else:
        hard_files[sub].append(f)

# Report counts
for s in range(11):
    print(f"  subwin {s}: {len(hard_files[s])} hard, {len(soft_files[s])} soft")

# Write concatenated msc files per type+subwin
# diploSHIC's openMsOutFileForSequentialReading expects an MS header line:
#   "program numSamples numReps"
# Our asc/*.ms files (after ascertainment) start with "//" directly — no header.
# The diploid re-encoding in ascertain.py emits pairs of haplotypes for ndip diploids.
# From simulate.py: NHAP=700, so ndip=350, and each rep has 700 haplotype rows.
NHAP = 700  # haplotypes per rep (matches simulate.py)

def cat_ms_files(file_list, out_path):
    """Concatenate ms files with proper MS header for diploSHIC."""
    reps = []
    for f in sorted(file_list):
        content = open(f).read().strip()
        if content:
            reps.append(content)
    nreps = len(reps)
    with open(out_path, "w") as out:
        # MS header: program numSamples numReps
        out.write(f"discoal {NHAP} {nreps} 44000\n")
        for rep_content in reps:
            out.write(rep_content + "\n")

for s in range(11):
    hp = f"{WD}/fvec/msc/hard_{s}.msc"
    sp = f"{WD}/fvec/msc/soft_{s}.msc"
    cat_ms_files(hard_files[s], hp)
    cat_ms_files(soft_files[s], sp)
    print(f"  wrote hard_{s}.msc ({len(hard_files[s])} sims), soft_{s}.msc ({len(soft_files[s])} sims)")

# Also write neutral concatenated msc
neutral_files = glob.glob(f"{WD}/sims/asc/neutral/neutral_s-1_r*.ms")
print(f"\nFound {len(neutral_files)} neutral files")
cat_ms_files(neutral_files, f"{WD}/fvec/msc/neutral.msc")
print(f"  wrote neutral.msc ({len(neutral_files)} sims)")
print("Done.")
