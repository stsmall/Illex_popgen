import numpy as np
from seqmodel import STATES
from seqmodel.transitions import hmm_log_transition
from seqmodel.mapcoord import split_blocks, deltas
def viterbi_hmm(emis_log, dg, dx, params):
    n=emis_log.shape[0]
    dp=np.full((n,4),-np.inf); back=np.zeros((n,4),int)
    dp[0]=np.log(np.array([0.985,0.010,0.004,0.001]))+emis_log[0]
    for t in range(1,n):
        A=hmm_log_transition(dg[t-1],dx[t-1],params)
        for s in range(4):
            cand=dp[t-1]+A[:,s]; back[t,s]=int(np.argmax(cand)); dp[t,s]=cand[back[t,s]]+emis_log[t,s]
    st=int(np.argmax(dp[-1])); score=float(dp[-1,st]); path=[st]
    for t in range(n-1,0,-1): st=back[t,st]; path.append(st)
    return [STATES[i] for i in reversed(path)], score
def decode_hmm_blocked(windows_df, map_df, chrom, emis_log, params, max_phys_gap):
    blocks=split_blocks(windows_df,map_df,chrom,max_phys_gap)
    path=["GAP"]*len(windows_df); total=0.0
    for (a,b) in blocks:
        wb=windows_df.iloc[a:b+1]
        dg,dx=deltas(wb,map_df,chrom) if b>a else (np.array([]),np.array([]))
        p,sc=viterbi_hmm(emis_log[a:b+1],dg,dx,params); total+=sc
        path[a:b+1]=p
    return path, total
import math
from seqmodel import STATE_IDX
def _lognorm(span, med, sig):
    span=max(span,1e-12); z=(math.log(span)-math.log(med))/sig
    return -0.5*z*z - math.log(span*sig*math.sqrt(2*math.pi))
def _dur_score(state, span, p):
    return {"N":0.0,"LL":_lognorm(span,p["ll_med"],p["ll_sig"]),
            "C":_lognorm(span,p["c_med"],p["c_sig"]),"LR":_lognorm(span,p["lr_med"],p["lr_sig"])}[state]
_ALLOWED={None:{"N":0.0},
  "N":{"N":0.0,"LL":None,"C":None},"LL":{"C":None,"N":None},"C":{"LR":None,"N":None},"LR":{"N":None}}
def _tp(prev,cur,p):
    if prev is None: return 0.0 if cur=="N" else p["abort_penalty"]
    table={"N":{"N":0.0,"LL":p["trans_penalty"],"C":p["narrow_penalty"]},
           "LL":{"C":p["trans_penalty"],"N":p["abort_penalty"]},
           "C":{"LR":p["trans_penalty"],"N":p["abort_penalty"]},"LR":{"N":p["trans_penalty"]}}
    return table.get(prev,{}).get(cur,-np.inf)
def _hsmm_score_path(idx_path, emis_log, gpos, p):  # reference scorer for tests
    from seqmodel import STATES
    segs=[]; i=0
    while i<len(idx_path):
        j=i
        while j+1<len(idx_path) and idx_path[j+1]==idx_path[i]: j+=1
        segs.append((STATES[idx_path[i]],i,j)); i=j+1
    sc=0.0; prev=None
    for (state,a,b) in segs:
        sc+=_tp(prev,state,p)
        span=(gpos[b]-gpos[a]) if b>a else 1e-6
        sc+=_dur_score(state,span,p)+float(emis_log[a:b+1,STATE_IDX[state]].sum()); prev=state
    return sc
def viterbi_hsmm(emis_log, gpos, params):
    from seqmodel import STATES
    n,ns=emis_log.shape; pref=np.vstack([np.zeros(ns),np.cumsum(emis_log,0)])
    dp=np.full((n,ns),-np.inf); bs=np.full((n,ns),-1,int); bp=np.full((n,ns),-1,int)
    for end in range(n):
        for s in range(ns):
            best=(-np.inf,-1,-1)
            for d in range(1,min(params["max_dur_windows"],end+1)+1):
                start=end-d+1
                seg=float(pref[end+1,s]-pref[start,s])
                span=(gpos[end]-gpos[start]) if end>start else 1e-6
                dur=_dur_score(STATES[s],span,params)
                if start==0:
                    sc=_tp(None,STATES[s],params)+dur+seg; prev=-1
                else:
                    cand=[dp[start-1,q]+_tp(STATES[q],STATES[s],params) for q in range(ns)]
                    prev=int(np.argmax(cand)); sc=cand[prev]+dur+seg
                if sc>best[0]: best=(sc,start,prev)
            dp[end,s],bs[end,s],bp[end,s]=best
    end=n-1; s=int(np.argmax(dp[end])); score=float(dp[end,s]); out=np.full(n,-1,int)
    if not np.isfinite(score):
        raise ValueError("viterbi_hsmm: no valid path (best score is -inf; check for -inf emissions / EPS clipping upstream)")
    while end>=0:
        start=int(bs[end,s])
        if start<0:
            raise ValueError("viterbi_hsmm: traceback reached an unreached state (start=-1)")
        out[start:end+1]=s; prev=bp[end,s]; end=start-1; s=prev
    return [STATES[i] for i in out], score
