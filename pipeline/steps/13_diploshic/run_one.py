"""Single-rep discoal runner: called by run_sims.sh via $PY run_one.py kind sub rep out"""
import sys, numpy as np
sys.path.insert(0, "/sietch_colab/data_share/illex/popgen_data/analysis/steps/13_diploshic")
import simulate as S, subprocess, os

kind, sub, rep, out = sys.argv[1], int(sys.argv[2]), int(sys.argv[3]), sys.argv[4]
rng = np.random.default_rng(hash((kind, sub, rep)) % 2**32)
if kind == "neutral":
    args = S.neutral(rng)
else:
    args, _ = S.sweep(rng, sub)
with open(out, "w") as f:
    subprocess.run([os.environ["DISCOAL"]] + args, stdout=f, check=True)
