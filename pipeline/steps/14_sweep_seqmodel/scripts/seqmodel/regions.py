import numpy as np, pandas as pd
from seqmodel import STATE_IDX
from seqmodel.transitions import hmm_log_transition
from seqmodel.decode import _dur_score, _tp

_N = STATE_IDX["N"]


def _region_score(seg, path, emis_log, params, model, dg, dx, gpos):
    """Discriminant S(R) over the region's contiguous window indices `seg`.

    Emission term always; when `params` is given, also the transition term (model='hmm')
    or the transition + log-normal-duration term (model='hsmm') that DESIGN §4 requires.
    Region windows are all non-N (LL/C/LR): the all-neutral path contributes A[N,N] per
    internal HMM step, or 0 for the HSMM (one zero-duration N run, no internal transition).
    """
    assigned = [STATE_IDX[path[k]] for k in seg]
    emis = float(sum(emis_log[k, assigned[m]] - emis_log[k, _N] for m, k in enumerate(seg)))
    if params is None:
        return emis
    if model == "hmm":
        trans = 0.0
        for m in range(1, len(seg)):
            k = seg[m]
            A = hmm_log_transition(dg[k - 1], dx[k - 1], params)
            trans += float(A[assigned[m - 1], assigned[m]] - A[_N, _N])
        return emis + trans
    if model == "hsmm":
        extra = 0.0
        prev = None
        m = 0
        while m < len(seg):
            st = path[seg[m]]
            a = seg[m]
            b = a
            while m + 1 < len(seg) and path[seg[m + 1]] == st:
                m += 1
                b = seg[m]
            span = float(gpos[b] - gpos[a]) if b > a else 1e-6
            extra += _dur_score(st, span, params)
            if prev is not None:
                extra += _tp(prev, st, params)
            prev = st
            m += 1
        return emis + extra
    raise ValueError(f"unknown model {model!r} (expected 'hmm' or 'hsmm')")


def extract_regions(pred_df, path, emis_log, collapsed, *,
                    params=None, model=None, dg=None, dx=None, gpos=None):
    path = list(path); rows = []; n = len(path); i = 0
    while i < n:
        if path[i] in ("N", "GAP"):
            i += 1; continue
        j = i
        while j + 1 < n and path[j + 1] not in ("N", "GAP"):
            j += 1
        seg = list(range(i, j + 1))
        if "C" in [path[k] for k in seg]:
            score = _region_score(seg, path, emis_log, params, model, dg, dx, gpos)
            cen = [k for k in seg if path[k] == "C"]
            lk = [k for k in seg if path[k] in ("LL", "LR")]
            he = float(pred_df["p_hard"].iloc[cen].sum() + pred_df["p_linked_hard"].iloc[lk].sum())
            se = float(pred_df["p_soft"].iloc[cen].sum() + pred_df["p_linked_soft"].iloc[lk].sum())
            rows.append(dict(
                chrom=pred_df["chrom"].iloc[i],
                start=int(pred_df["start"].iloc[i]), end=int(pred_df["end"].iloc[j]),
                center_start=int(pred_df["start"].iloc[cen[0]]), center_end=int(pred_df["end"].iloc[cen[-1]]),
                n_windows=len(seg), n_center=len(cen),
                score=score, min_pN=float(collapsed[i:j + 1, 0].min()), max_pC=float(collapsed[cen, 2].max()),
                hard_evidence=he, soft_evidence=se, hard_fraction=(he / (he + se) if he + se > 0 else np.nan)))
        i = j + 1
    cols = ["chrom", "start", "end", "center_start", "center_end", "n_windows", "n_center",
            "score", "min_pN", "max_pC", "hard_evidence", "soft_evidence", "hard_fraction"]
    return pd.DataFrame(rows, columns=cols)
