import numpy as np, pandas as pd
from seqmodel.emissions import collapse, emission_logprob

def _pred(p):  # p = list of (pn,plh,pls,ps,ph)
    return pd.DataFrame(p, columns=["p_neutral","p_linked_hard","p_linked_soft","p_hard","p_soft"])

def test_collapse_sums_and_columns():
    c=collapse(_pred([(0.9,0.02,0.03,0.03,0.02)]))
    assert c.shape==(1,3) and abs(c.sum()-1)<1e-9
    assert np.allclose(c[0], [0.9,0.05,0.05], atol=1e-9)

def test_ll_lr_share_linked_emission():
    e=emission_logprob(collapse(_pred([(0.2,0.2,0.2,0.2,0.2)])), transform="raw")
    assert e.shape==(1,4)
    assert e[0,1]==e[0,3]  # LL == LR

def test_prior_correction_upweights_neutral():
    c=collapse(_pred([(0.2,0.2,0.2,0.2,0.2)]))  # equal collapsed 0.2/0.4/0.4
    raw=emission_logprob(c,transform="raw"); pc=emission_logprob(c,transform="prior_corrected")
    # prior_corrected divides by training prior => N (÷0.2) rises most vs C (÷0.4)
    assert (pc[0,0]-raw[0,0]) > (pc[0,2]-raw[0,2])
