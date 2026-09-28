import numpy as np
from seqmodel import EPS, STATE_IDX as SI
def _norm(d):
    tot=sum(d.values()); return {k:v/tot for k,v in d.items()}
def hmm_log_transition(dg, dx, params):
    for _k in ("narrow_entry", "abort_left", "abort_center", "skip_weight"):
        _v = params[_k]
        if not (0.0 <= _v <= 1.0):
            raise ValueError(f"transition param {_k}={_v} must be in [0, 1]")
    if params["abort_left"] + params["skip_weight"] > 1.0:
        raise ValueError(
            f"abort_left + skip_weight = {params['abort_left'] + params['skip_weight']} must be <= 1 "
            "(else LL->C weight is negative)"
        )
    d_eff = max(dg,0.0) + params["beta"]*max(dx,0.0)
    lam = np.array([params["lambda_n"],params["lambda_ll"],params["lambda_c"],params["lambda_lr"]])
    stay = np.clip(np.exp(-d_eff/np.maximum(lam,EPS)), EPS, 1-EPS)
    sw = params["skip_weight"]
    exits = {
        SI["N"]:  _norm({SI["LL"]:1-params["narrow_entry"], SI["C"]:params["narrow_entry"]}),
        SI["LL"]: _norm({SI["C"]:1-params["abort_left"]-sw, SI["N"]:params["abort_left"], SI["LR"]:sw}),
        SI["C"]:  _norm({SI["LR"]:1-params["abort_center"], SI["N"]:params["abort_center"]}),
        SI["LR"]: {SI["N"]:1.0},
    }
    A=np.full((4,4),EPS)
    for s in range(4):
        A[s,s]=stay[s]
        for dest,w in exits[s].items(): A[s,dest]=max((1-stay[s])*w, EPS)
        A[s]/=A[s].sum()
    return np.log(A)
