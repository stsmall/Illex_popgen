import numpy as np
from seqmodel import EPS
def collapse(pred_df):
    pN=pred_df["p_neutral"].to_numpy(float)
    pL=(pred_df["p_linked_hard"]+pred_df["p_linked_soft"]).to_numpy(float)
    pC=(pred_df["p_hard"]+pred_df["p_soft"]).to_numpy(float)
    P=np.clip(np.column_stack([pN,pL,pC]),EPS,None); return P/P.sum(1,keepdims=True)
def emission_logprob(collapsed, transform="prior_corrected", priors=(0.2,0.4,0.4), temperature=1.0):
    P=np.clip(collapsed,EPS,None).astype(float)
    if transform=="temperature":
        P=P**(1.0/temperature); P=P/P.sum(1,keepdims=True)
    logp=np.log(np.clip(P,EPS,1.0))
    if transform=="prior_corrected":
        logp=logp-np.log(np.asarray(priors))
    # map [N,L,C] -> [N,LL,C,LR]
    return np.column_stack([logp[:,0], logp[:,1], logp[:,2], logp[:,1]])
