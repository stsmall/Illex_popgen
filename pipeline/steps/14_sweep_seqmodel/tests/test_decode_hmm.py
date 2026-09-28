# tests/test_decode_hmm.py
import numpy as np, pandas as pd
from seqmodel.decode import viterbi_hmm, decode_hmm_blocked
from seqmodel.transitions import hmm_log_transition
P=dict(lambda_n=0.05,lambda_ll=2e-3,lambda_c=1e-3,lambda_lr=2e-3,beta=1e-9,
       skip_weight=0.01,narrow_entry=0.02,abort_left=0.05,abort_center=0.02)
def _emis(seq):  # strong evidence per state index
    e=np.full((len(seq),4),np.log(1e-4))
    for i,s in enumerate(seq): e[i,s]=np.log(0.999)
    return e
def test_recovers_canonical_sweep_path():
    truth=[0,0,1,1,2,2,3,3,0,0]  # N N LL LL C C LR LR N N
    n=len(truth); dg=np.full(n-1,5e-4); dx=np.full(n-1,4000.)
    path,score=viterbi_hmm(_emis(truth),dg,dx,P)
    assert path==["N","N","LL","LL","C","C","LR","LR","N","N"]
def test_all_neutral_stays_neutral():
    truth=[0]*8; n=8; path,_=viterbi_hmm(_emis(truth),np.full(n-1,1e-3),np.full(n-1,4000.),P)
    assert set(path)=={"N"}
def test_blocked_marks_gap_and_decodes_each_block():
    w=pd.DataFrame({"chrom":"1","start":[0,4000,5_000_000],"end":[4000,8000,5_004_000]})
    m=pd.DataFrame([("1",0,10_000_000,1e-8)],columns=["chrom","start","end","rate"])
    emis=_emis([0,2,0])
    path,score=decode_hmm_blocked(w,m,"1",emis,P,max_phys_gap=1_000_000)
    assert len(path)==3 and path[2] in ("N","C","LL","LR")  # 3rd is its own block, decoded
