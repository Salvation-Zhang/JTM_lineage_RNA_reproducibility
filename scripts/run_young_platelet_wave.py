from pathlib import Path
import gzip, glob, itertools, re
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy import stats

ROOT=Path(r"D:\TISTA_JTH_manuscript")
EXT=ROOT/"young_platelet_analysis"/"GSE126448"
PRIMARY=ROOT/"primary_raw"
OUT=ROOT/"young_platelet_analysis"/"results"
OUT.mkdir(parents=True,exist_ok=True)

PLATELET=["PPBP","PF4","GP1BA","GP9","ITGA2B","ITGB3","TUBB1","TREML1","RGS18","SDPR","SPARC","CLU"]

def load_two_col_gz(path):
    x=pd.read_csv(path,sep="\t",header=None,index_col=0,compression="gzip",dtype={1:str})
    x=pd.to_numeric(x.iloc[:,0],errors="coerce").fillna(0); x.name=Path(path).name.split('_')[0]
    return x

def zscore_signature(logcpm,genes):
    g=[x for x in genes if x in logcpm.index]
    z=logcpm.loc[g].sub(logcpm.loc[g].mean(axis=1),axis=0).div(logcpm.loc[g].std(axis=1,ddof=1).replace(0,np.nan),axis=0)
    return z.mean(axis=0),g

def exact_signflip(x):
    x=np.asarray(x,float); obs=abs(x.mean()); null=[]
    for s in itertools.product([-1,1],repeat=len(x)): null.append(abs(np.mean(x*np.asarray(s))))
    return np.mean(np.asarray(null)>=obs-1e-12)

def exact_between(a,b):
    x=np.r_[a,b]; k=len(a); obs=abs(np.mean(a)-np.mean(b)); vals=[]
    for ix in itertools.combinations(range(len(x)),k):
        m=np.zeros(len(x),bool); m[list(ix)]=True; vals.append(abs(x[m].mean()-x[~m].mean()))
    e=np.sum(np.asarray(vals)>=obs-1e-12)
    return (e+1)/(len(vals)+1)

# External paired immature/mature platelet RNA-seq.
files=sorted(glob.glob(str(EXT/"GSM*_HTSeqCount_Lib-DBO-*.txt.gz")))
ext=pd.concat([load_two_col_gz(f) for f in files],axis=1).fillna(0)
stages=np.array(["immature" if i%2==0 else "mature" for i in range(len(files))])
pairs=np.repeat(np.arange(len(files)//2),2)
lib=ext.sum(axis=0); elog=np.log2(ext.add(.5).div(lib.add(1),axis=1)*1e6+1)
delta=pd.DataFrame({f"pair{i+1}":elog.iloc[:,2*i]-elog.iloc[:,2*i+1] for i in range(4)})
external=pd.DataFrame({"gene":delta.index,"mean_immature_minus_mature":delta.mean(axis=1),
                       "min_pair_delta":delta.min(axis=1),"max_pair_delta":delta.max(axis=1),
                       "sign_consistency":np.sign(delta).replace(0,np.nan).mean(axis=1).abs()})
external["mean_logcpm"]=elog.mean(axis=1).values
external=external.sort_values("mean_immature_minus_mature",ascending=False)

# Lock two 50-gene tails, requiring direction in all four donor pairs and adequate expression.
young=external[(external.min_pair_delta>0)&(external.mean_logcpm>1)].head(50).gene.tolist()
old=external[(external.max_pair_delta<0)&(external.mean_logcpm>1)].tail(50).gene.tolist()
external.to_csv(OUT/"GSE126448_paired_platelet_age_gene_statistics.tsv",sep="\t",index=False)
pd.DataFrame({"young_gene":pd.Series(young),"mature_gene":pd.Series(old)}).to_csv(OUT/"locked_platelet_age_signatures.tsv",sep="\t",index=False)

# Reconstruct GSE112278 matrix from GEO two-column files.
pfiles=sorted(glob.glob(str(PRIMARY/"GSM*.gz")))
series=[]; meta=[]
for f in pfiles:
    name=Path(f).name
    m=re.search(r'(GSM\d+)_((?:R|NR)\d+)_(pre|1wk|1mon)_',name)
    if not m: continue
    s=load_two_col_gz(f); s.name=m.group(1); series.append(s)
    patient=m.group(2); meta.append({"gsm":m.group(1),"patient":patient,"response":"nonresponder" if patient.startswith("NR") else "responder","time":m.group(3)})
counts=pd.concat(series,axis=1).fillna(0); meta=pd.DataFrame(meta)
lib=counts.sum(axis=0); logcpm=np.log2(counts.add(.5).div(lib.add(1),axis=1)*1e6+1)
platelet,pg=zscore_signature(logcpm,PLATELET)
ys,yg=zscore_signature(logcpm,young)
os,og=zscore_signature(logcpm,old)
contrast=ys-os

# Remove total platelet abundance from the age contrast without using outcome labels.
X=np.column_stack([np.ones(len(platelet)),platelet.values])
b=np.linalg.lstsq(X,contrast.values,rcond=None)[0]
resid=pd.Series(contrast.values-X@b,index=contrast.index,name="platelet_age_residual")
for name,s in [("platelet_abundance",platelet),("young_score",ys),("mature_score",os),("young_minus_mature",contrast),("platelet_age_residual",resid)]: meta[name]=meta.gsm.map(s)
meta["time_order"]=meta.time.map({"pre":0,"1wk":1,"1mon":2})
meta.sort_values(["patient","time_order"]).to_csv(OUT/"GSE112278_platelet_age_scores.tsv",sep="\t",index=False)

def changes(score,t1="pre",t2="1wk"):
    w=meta.pivot(index=["patient","response"],columns="time",values=score).dropna(subset=[t1,t2]); w["delta"]=w[t2]-w[t1]
    return w.reset_index()

rows=[]
for score in ["platelet_abundance","young_score","mature_score","young_minus_mature","platelet_age_residual"]:
    ch=changes(score); a=ch[ch.response=="responder"].delta.values; n=ch[ch.response=="nonresponder"].delta.values
    rows.append({"score":score,"responder_n":len(a),"responder_mean_pre_to_week1":a.mean(),"responder_signflip_p":exact_signflip(a),
                 "nonresponder_n":len(n),"nonresponder_mean_pre_to_week1":n.mean(),"between_exact_p":exact_between(a,n)})
    wm=changes(score,"1wk","1mon"); r=wm[wm.response=="responder"].delta.values
    rows[-1].update({"complete_responder_n":len(r),"responder_mean_week1_to_month1":r.mean(),"week1_to_month1_signflip_p":exact_signflip(r) if len(r) else np.nan})
summary=pd.DataFrame(rows); summary.to_csv(OUT/"platelet_age_longitudinal_tests.tsv",sep="\t",index=False)

# Gene-set overlap and external behavior of locked platelet markers.
pd.DataFrame({"young_signature_n":[len(yg)],"mature_signature_n":[len(og)],"young_overlap_locked_platelet":[','.join(sorted(set(yg)&set(pg)))],
              "mature_overlap_locked_platelet":[','.join(sorted(set(og)&set(pg)))],"correlation_age_contrast_vs_platelet_abundance":[np.corrcoef(contrast,platelet)[0,1]]}).to_csv(OUT/"platelet_age_signature_audit.tsv",sep="\t",index=False)

plt.rcParams.update({"font.family":"DejaVu Sans","font.size":9,"axes.spines.top":False,"axes.spines.right":False,"savefig.dpi":600})
fig,axs=plt.subplots(1,3,figsize=(11,4))
# External paired validation
for i in range(4):
    a=elog.iloc[:,2*i]; b0=elog.iloc[:,2*i+1]
    ext_y=a.loc[yg].mean()-b0.loc[yg].mean(); ext_o=a.loc[og].mean()-b0.loc[og].mean()
    axs[0].plot([0,1],[ext_y,ext_o],"o-",alpha=.7)
axs[0].axhline(0,color="#A0AEC0",lw=.8); axs[0].set_xticks([0,1],["Young signature","Mature signature"]); axs[0].set_ylabel("Immature - mature score"); axs[0].set_title("External paired platelets")
# Primary residual trajectories
colors={"responder":"#2B6CB0","nonresponder":"#C53030"}
for grp in ["responder","nonresponder"]:
    d=meta[meta.response==grp]
    for _,g in d.groupby("patient"):
        g=g.sort_values("time_order"); axs[1].plot(g.time_order,g.platelet_age_residual,"-o",color=colors[grp],alpha=.35,ms=3)
    mn=d.groupby("time_order").platelet_age_residual.mean(); axs[1].plot(mn.index,mn.values,"-o",color=colors[grp],lw=2.5,label=grp.capitalize())
axs[1].axhline(0,color="#A0AEC0",lw=.8); axs[1].set_xticks([0,1,2],["Baseline","Week 1","Month 1"]); axs[1].set_ylabel("Age score residualized for platelet abundance"); axs[1].set_title("Putative young-platelet wave"); axs[1].legend(frameon=False,fontsize=7)
# Compare changes
show=summary[summary.score.isin(["platelet_abundance","young_minus_mature","platelet_age_residual"])]
x=np.arange(len(show)); axs[2].bar(x-.18,show.responder_mean_pre_to_week1,width=.36,label="Responder",color="#2B6CB0"); axs[2].bar(x+.18,show.nonresponder_mean_pre_to_week1,width=.36,label="Nonresponder",color="#C53030")
axs[2].axhline(0,color="#A0AEC0",lw=.8); axs[2].set_xticks(x,["Platelet\nabundance","Young-mature\ncontrast","Age residual"],rotation=10); axs[2].set_ylabel("Week 1 - baseline change"); axs[2].set_title("Separation from abundance"); axs[2].legend(frameon=False,fontsize=7)
for i,a in enumerate(axs): a.text(-.13,1.04,chr(65+i),transform=a.transAxes,fontweight="bold",fontsize=12)
fig.tight_layout()
for extn in ["png","pdf","svg"]: fig.savefig(OUT/f"Young_platelet_wave_exploratory.{extn}",bbox_inches="tight")
plt.close(fig)
print(summary.to_string(index=False)); print('young genes',yg[:10]); print('mature genes',og[:10])
