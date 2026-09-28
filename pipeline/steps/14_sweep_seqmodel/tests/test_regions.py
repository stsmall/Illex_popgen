import numpy as np, pandas as pd
from seqmodel.regions import extract_regions
from seqmodel.emissions import collapse, emission_logprob
def _pred(rows): return pd.DataFrame(rows, columns=["chrom","start","end","p_neutral","p_linked_hard","p_linked_soft","p_hard","p_soft"])
def test_extracts_one_region_with_center():
    rows=[("1",0,4000,.95,.02,.01,.01,.01),("1",4000,8000,.1,.4,.1,.3,.1),
          ("1",8000,12000,.05,.05,.05,.7,.15),("1",12000,16000,.1,.1,.4,.1,.3),
          ("1",16000,20000,.95,.02,.01,.01,.01)]
    df=_pred(rows); c=collapse(df); e=emission_logprob(c,transform="raw")
    path=["N","LL","C","LR","N"]
    R=extract_regions(df,path,e,c)
    assert len(R)==1 and R.n_center.iloc[0]==1 and R.score.iloc[0]>0
    assert 0.0<=R.hard_fraction.iloc[0]<=1.0
def test_no_region_without_center():
    df=_pred([("1",0,4000,.2,.4,.1,.2,.1),("1",4000,8000,.2,.4,.1,.2,.1)])
    c=collapse(df); e=emission_logprob(c,transform="raw")
    assert len(extract_regions(df,["LL","LL"],e,c))==0

def test_hmm_region_score_adds_transition_terms():
    import numpy as np
    from seqmodel.transitions import hmm_log_transition
    from seqmodel import STATE_IDX as SI
    rows=[("1",0,4000,.95,.02,.01,.01,.01),("1",4000,8000,.1,.4,.1,.3,.1),
          ("1",8000,12000,.05,.05,.05,.7,.15),("1",12000,16000,.1,.1,.4,.1,.3),
          ("1",16000,20000,.95,.02,.01,.01,.01)]
    df=_pred(rows); c=collapse(df); e=emission_logprob(c,transform="raw")
    path=["N","LL","C","LR","N"]
    P=dict(lambda_n=0.05,lambda_ll=2e-3,lambda_c=1e-3,lambda_lr=2e-3,beta=1e-9,
           skip_weight=0.01,narrow_entry=0.02,abort_left=0.05,abort_center=0.02)
    dg=np.full(4,5e-4); dx=np.full(4,4000.)
    R0=extract_regions(df,path,e,c)
    R1=extract_regions(df,path,e,c,params=P,model="hmm",dg=dg,dx=dx)
    seg=[1,2,3]; a=[SI[path[k]] for k in seg]; tdiff=0.0
    for m in range(1,len(seg)):
        k=seg[m]; A=hmm_log_transition(dg[k-1],dx[k-1],P)
        tdiff+=A[a[m-1],a[m]]-A[SI["N"],SI["N"]]
    assert tdiff != 0.0
    assert abs((R1.score.iloc[0]-R0.score.iloc[0]) - tdiff) < 1e-9

def test_hsmm_region_score_adds_duration_and_transition():
    import numpy as np
    from seqmodel.decode import _dur_score, _tp
    rows=[("1",0,4000,.95,.02,.01,.01,.01),("1",4000,8000,.1,.4,.1,.3,.1),
          ("1",8000,12000,.05,.05,.05,.7,.15),("1",12000,16000,.1,.1,.4,.1,.3),
          ("1",16000,20000,.95,.02,.01,.01,.01)]
    df=_pred(rows); c=collapse(df); e=emission_logprob(c,transform="raw")
    path=["N","LL","C","LR","N"]
    HP=dict(ll_med=2e-3,ll_sig=0.6,c_med=1e-3,c_sig=0.6,lr_med=2e-3,lr_sig=0.6,
            max_dur_windows=6,trans_penalty=-0.5,narrow_penalty=-2.0,abort_penalty=-3.0)
    gpos=np.array([0.,2e-3,3e-3,5e-3,7e-3])
    R0=extract_regions(df,path,e,c)
    R1=extract_regions(df,path,e,c,params=HP,model="hsmm",gpos=gpos)
    dur=_dur_score("LL",1e-6,HP)+_dur_score("C",1e-6,HP)+_dur_score("LR",1e-6,HP)
    trans=_tp("LL","C",HP)+_tp("C","LR",HP)
    assert abs((R1.score.iloc[0]-R0.score.iloc[0]) - (dur+trans)) < 1e-9
