"""Does chr2 karyotype track capture depth? Individual BB-allele count (AA=0, AB=1, BB=2) vs recorded depth.
3K was jig-caught at 3-5 m (Baker Table 1); 4X was bulk-sampled with no depth recorded -> excluded.
Tests: (1) across individuals, Spearman + ordinal-free permutation; (2) across hauls (unique division x date x
depth), weighted regression of arrangement frequency on depth; (3) within divisions that sampled >1 depth,
permutation of depth labels within division (depth effect net of geography)."""
import pandas as pd, numpy as np
from scipy import stats
rng=np.random.default_rng(1)
m=pd.read_csv("/sietch_colab/data_share/illex/popgen_data/seq_data/baker_2025/docs/Squid_Meta_Sept2024_Simple.txt",sep="\t",dtype=str)
k=pd.read_csv("/sietch_colab/data_share/illex/popgen_data/analysis/steps/03_karyotype/karyotypes.tsv",sep="\t")
d=m.merge(k[["sample","karyotype"]],left_on="ID",right_on="sample")
d=d[~d.Zone.isin(["4X"])&d.Depth.notna()].copy()
d["depth"]=d.Depth.astype(float); d["b"]=d.karyotype.map({"AA":0,"AB":1,"BB":2}); d=d.dropna(subset=["b"])
out=[]
rho,p=stats.spearmanr(d.depth,d.b)
perm=np.array([stats.spearmanr(rng.permutation(d.depth),d.b)[0] for _ in range(5000)])
out.append(f"(1) individuals n={len(d)}, depth {d.depth.min():.0f}-{d.depth.max():.0f} m, {d.Zone.nunique()} divisions: Spearman rho={rho:+.3f}, p={p:.3f}, permutation p={(abs(perm)>=abs(rho)).mean():.3f}")
h=d.groupby(["Zone","Collection_Date","depth"]).agg(n=("b","size"),q=("b",lambda x:x.sum()/(2*len(x)))).reset_index()
h=h[h.n>=5]
w=h.n; X=np.c_[np.ones(len(h)),h.depth]; W=np.diag(w)
beta=np.linalg.solve(X.T@W@X,X.T@W@h.q); r=h.q-X@beta
sig2=(w*r**2).sum()/(len(h)-2); se=np.sqrt(sig2*np.linalg.inv(X.T@W@X)[1,1]); t=beta[1]/se
out.append(f"(2) hauls with n>=5: {len(h)}; weighted slope of BB frequency on depth = {beta[1]*100:+.4f} per 100 m (SE {se*100:.4f}), t={t:+.2f}, p={2*stats.t.sf(abs(t),len(h)-2):.3f}; freq range {h.q.min():.2f}-{h.q.max():.2f}")
out.append("    hauls: "+"; ".join(f"{r.Zone} {r.depth:.0f}m n={r.n} p={r.q:.2f}" for r in h.itertuples()))
multi=d.groupby("Zone").depth.nunique(); zs=multi[multi>1].index; dd=d[d.Zone.isin(zs)].copy()
def within(dd):
    return sum(stats.spearmanr(g.depth,g.b)[0]*len(g) for _,g in dd.groupby("Zone") if g.depth.nunique()>1)/len(dd)
obs=within(dd)
null=[]
for _ in range(5000):
    x=dd.copy(); x["depth"]=x.groupby("Zone").depth.transform(lambda s: rng.permutation(s.values)); null.append(within(x))
null=np.array(null)
out.append(f"(3) within-division, {len(zs)} divisions with >1 depth ({', '.join(zs)}), n={len(dd)}: n-weighted mean Spearman rho={obs:+.3f}, within-division permutation p={(abs(null)>=abs(obs)).mean():.3f}")
open("karyo_depth_test.txt","w").write("\n".join(out)+"\n"); print("\n".join(out))
