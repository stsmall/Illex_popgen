import numpy as np
from seqmodel.transitions import hmm_log_transition
P=dict(lambda_n=0.01,lambda_ll=1e-3,lambda_c=5e-4,lambda_lr=1e-3,beta=1e-9,
       skip_weight=0.02,narrow_entry=0.05,abort_left=0.1,abort_center=0.05)
def test_rows_sum_to_one():
    A=np.exp(hmm_log_transition(1e-4,1e4,P)); assert np.allclose(A.sum(1),1,atol=1e-9)
def test_zero_dg_does_not_freeze():
    A=np.exp(hmm_log_transition(0.0,1e5,P))  # dg=0 but big physical gap
    assert A[0,0] < 0.9999   # via beta*dx, staying prob < 1 so state can change
def test_large_distance_allows_skip_N_to_C():
    A=np.exp(hmm_log_transition(1.0,1e7,P))  # essentially certain to leave N
    assert A[0,2] > 0        # N->C skip has positive probability
def test_stay_prob_decays_with_distance():
    near=np.exp(hmm_log_transition(1e-5,1e3,P))[2,2]; far=np.exp(hmm_log_transition(1e-2,1e6,P))[2,2]
    assert near > far
def test_invalid_exit_params_raise():
    import pytest
    bad = dict(P); bad["abort_left"] = 0.6; bad["skip_weight"] = 0.6
    with pytest.raises(ValueError):
        hmm_log_transition(1e-4, 1e4, bad)
    bad2 = dict(P); bad2["narrow_entry"] = 1.5
    with pytest.raises(ValueError):
        hmm_log_transition(1e-4, 1e4, bad2)
