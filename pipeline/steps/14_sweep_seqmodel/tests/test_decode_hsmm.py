# tests/test_decode_hsmm.py
import numpy as np, itertools
from seqmodel.decode import viterbi_hsmm
HP=dict(ll_med=2e-3,ll_sig=0.6,c_med=1e-3,c_sig=0.6,lr_med=2e-3,lr_sig=0.6,
        max_dur_windows=6,trans_penalty=-0.5,narrow_penalty=-2.0,abort_penalty=-3.0)
def _emis(seq):
    e=np.full((len(seq),4),np.log(1e-3))
    for i,s in enumerate(seq): e[i,s]=np.log(0.999)
    return e
def _bruteforce(emis_log,gpos,params):
    from seqmodel.decode import _hsmm_score_path
    n=emis_log.shape[0]; best=(-np.inf,None)
    for path in itertools.product(range(4),repeat=n):
        sc=_hsmm_score_path(list(path),emis_log,gpos,params)
        if sc>best[0]: best=(sc,list(path))
    return best
def test_matches_bruteforce_small():
    seq=[0,1,2,3,0]; g=np.array([0,2e-3,3e-3,5e-3,7e-3]); e=_emis(seq)
    _,p=_bruteforce(e,g,HP); path,_=viterbi_hsmm(e,g,HP)
    from seqmodel import STATE_IDX
    assert [STATE_IDX[s] for s in path]==p
def test_single_window_center_ok():
    e=_emis([0,2,0]); g=np.array([0,1e-3,2e-3]); path,_=viterbi_hsmm(e,g,HP)
    assert path[1] in ("C","N")
def test_hsmm_all_neg_inf_raises():
    import numpy as np, pytest
    e = np.full((3, 4), -np.inf); g = np.array([0., 1e-3, 2e-3])
    with pytest.raises(ValueError):
        viterbi_hsmm(e, g, HP)
def test_hsmm_matches_bruteforce_random_trials():
    import numpy as np
    from seqmodel import STATE_IDX
    rng = np.random.default_rng(0)
    for _ in range(30):
        n = 5
        e = np.log(rng.random((n, 4)) + 1e-3)                       # finite log-emissions
        g = np.cumsum(np.concatenate([[0.0], rng.random(n - 1) * 3e-3]))  # strictly increasing gpos
        sc_bf, p_bf = _bruteforce(e, g, HP)
        path, sc = viterbi_hsmm(e, g, HP)
        assert [STATE_IDX[s] for s in path] == p_bf
        assert abs(sc - sc_bf) < 1e-6
