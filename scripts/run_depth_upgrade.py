from pathlib import Path
import glob, itertools, re
import numpy as np
import pandas as pd
try:
    import matplotlib.pyplot as plt
    HAS_MPL=True
except ModuleNotFoundError:
    HAS_MPL=False

ROOT=Path(r"D:\TISTA_JTH_manuscript"); PRIMARY=ROOT/"primary_raw"; EXT=ROOT/"young_platelet_analysis"/"GSE126448"; OUT=ROOT/"upgrade_analysis"/"results"
OUT.mkdir(parents=True,exist_ok=True); RNG=np.random.default_rng(20260724)
PLATELET=["PPBP","PF4","GP1BA","GP9","ITGA2B","ITGB3","TUBB1","TREML1","RGS18","SDPR","SPARC","CLU"]

def load_two_col_gz(path):
    x=pd.read_csv(path,sep="\t",header=None,index_col=0,compression="gzip"); x=pd.to_numeric(x.iloc[:,0],errors="coerce").fillna(0); x.name=Path(path).name.split('_')[0]; return x
def zscore_signature(logcpm,genes):
    g=[x for x in genes if x in logcpm.index]; z=logcpm.loc[g].sub(logcpm.loc[g].mean(axis=1),axis=0).div(logcpm.loc[g].std(axis=1,ddof=1).replace(0,np.nan),axis=0); return z.mean(axis=0),g
def exact_signflip(x):
    x=np.asarray(x,float); obs=abs(x.mean()); vals=[abs(np.mean(x*np.asarray(s))) for s in itertools.product([-1,1],repeat=len(x))]; return np.mean(np.asarray(vals)>=obs-1e-12)
def exact_between(a,b):
    x=np.r_[a,b]; k=len(a); obs=abs(np.mean(a)-np.mean(b)); vals=[]
    for ix in itertools.combinations(range(len(x)),k):
        m=np.zeros(len(x),bool); m[list(ix)]=True; vals.append(abs(x[m].mean()-x[~m].mean()))
    return (np.sum(np.asarray(vals)>=obs-1e-12)+1)/(len(vals)+1)
def residualize(y,x,fit_mask=None):
    y=np.asarray(y,float); x=np.asarray(x,float); fit_mask=np.ones(len(y),bool) if fit_mask is None else fit_mask; X=np.column_stack([np.ones(len(x)),x]); b=np.linalg.lstsq(X[fit_mask],y[fit_mask],rcond=None)[0]; return y-X@b,b
def paired_table(meta,value,t1,t2):
    d=meta[["gsm","patient","response","time"]].copy(); d["value"]=d.gsm.map(value); w=d.pivot(index=["patient","response"],columns="time",values="value").dropna(subset=[t1,t2]); w["delta"]=w[t2]-w[t1]; return w.reset_index()
def group_design(response): return np.column_stack([np.ones(len(response)),(np.asarray(response)=="responder").astype(float)])
def multivariate_partial_r2(Y,response,platelet_delta):
    X0=group_design(response); X1=np.column_stack([X0,np.asarray(platelet_delta,float)]); E0=Y-X0@np.linalg.lstsq(X0,Y,rcond=None)[0]; E1=Y-X1@np.linalg.lstsq(X1,Y,rcond=None)[0]; s0=np.sum(E0**2); s1=np.sum(E1**2); return 1-s1/s0,s0,s1
def coupling_stats(Y,response,pdelt):
    r=(response=="responder").astype(float); yc=Y-Y.mean(axis=0); pc=pdelt-pdelt.mean()
    den=np.sqrt(np.sum(yc**2,axis=0)*np.sum(pc**2)); cor=np.divide(pc@yc,den,out=np.full(Y.shape[1],np.nan),where=den>0)
    X0=np.column_stack([np.ones(len(r)),r]); X1=np.column_stack([np.ones(len(r)),r,pdelt])
    b0=np.linalg.lstsq(X0,Y,rcond=None)[0][1]; b1=np.linalg.lstsq(X1,Y,rcond=None)[0][1]
    att=np.divide(b0-b1,b0,out=np.full_like(b0,np.nan),where=np.abs(b0)>1e-12)
    return np.nanmedian(cor),np.nanmedian(att)

# Primary cohort
series=[]; rows=[]
for f in sorted(glob.glob(str(PRIMARY/"GSM*.gz"))):
    m=re.search(r'(GSM\d+)_((?:R|NR)\d+)_(pre|1wk|1mon)_',Path(f).name)
    if not m: continue
    s=load_two_col_gz(f); s.name=m.group(1); series.append(s); rows.append({"gsm":m.group(1),"patient":m.group(2),"response":"nonresponder" if m.group(2).startswith("NR") else "responder","time":m.group(3)})
meta=pd.DataFrame(rows); counts=pd.concat(series,axis=1).fillna(0); logcpm=np.log2(counts.add(.5).div(counts.sum(axis=0).add(1),axis=1)*1e6+1); platelet,pg=zscore_signature(logcpm,PLATELET)

# Overlap-free multivariate variance attribution
pre1=meta[meta.time.isin(["pre","1wk"])].pivot(index=["patient","response"],columns="time",values="gsm").dropna()
delta=pd.DataFrame({idx[0]:logcpm[row["1wk"]]-logcpm[row["pre"]] for idx,row in pre1.iterrows()}); resp_pat=[idx[0] for idx in pre1.index if idx[1]=="responder"]
top208=delta[resp_pat].mean(axis=1).nlargest(208).index.tolist(); analysis_genes=[g for g in top208 if g not in set(pg)]; patient_order=[idx[0] for idx in pre1.index]; response=np.array([idx[1] for idx in pre1.index])
Y=delta.loc[analysis_genes,patient_order].T.values; pdelt=np.array([platelet[row["1wk"]]-platelet[row["pre"]] for _,row in pre1.iterrows()]); pr2,sse0,sse1=multivariate_partial_r2(Y,response,pdelt); corr0,att0=coupling_stats(Y,response,pdelt)

# Leave-one-patient-out incremental predictive R2 (fixed outcome-gene set).
pred0=np.zeros_like(Y); pred1=np.zeros_like(Y)
for i in range(len(response)):
    tr=np.arange(len(response))!=i; x0tr=group_design(response[tr]); x0te=group_design(response[[i]])
    x1tr=np.column_stack([x0tr,pdelt[tr]]); x1te=np.column_stack([x0te,pdelt[[i]]])
    pred0[i]=x0te@np.linalg.lstsq(x0tr,Y[tr],rcond=None)[0]
    pred1[i]=x1te@np.linalg.lstsq(x1tr,Y[tr],rcond=None)[0]
cv_sse0=np.sum((Y-pred0)**2); cv_sse1=np.sum((Y-pred1)**2); cv_incremental_r2=1-cv_sse1/cv_sse0

# Fully nested LOPO: reselect responder-induced genes within each training fold.
nested_sse0=0.0; nested_sse1=0.0; nested_gene_n=[]
for i,held in enumerate(patient_order):
    train_pat=[p for p in patient_order if p!=held]; train_resp=[p for p in train_pat if response[patient_order.index(p)]=="responder"]
    fold_top=delta[train_resp].mean(axis=1).nlargest(208).index.tolist(); fold_genes=[g for g in fold_top if g not in set(pg)]; nested_gene_n.append(len(fold_genes))
    tr=np.arange(len(response))!=i; ytr=delta.loc[fold_genes,[patient_order[k] for k in np.where(tr)[0]]].T.values; yte=delta.loc[fold_genes,held].values
    x0tr=group_design(response[tr]); x0te=group_design(response[[i]]); x1tr=np.column_stack([x0tr,pdelt[tr]]); x1te=np.column_stack([x0te,pdelt[[i]]])
    p0=(x0te@np.linalg.lstsq(x0tr,ytr,rcond=None)[0]).ravel(); p1=(x1te@np.linalg.lstsq(x1tr,ytr,rcond=None)[0]).ravel()
    nested_sse0+=np.sum((yte-p0)**2); nested_sse1+=np.sum((yte-p1)**2)
nested_cv_r2=1-nested_sse1/nested_sse0

# Stratified patient bootstrap
ri=np.where(response=="responder")[0]; ni=np.where(response=="nonresponder")[0]; boot=[]
for _ in range(10000):
    ix=np.r_[RNG.choice(ri,len(ri),replace=True),RNG.choice(ni,len(ni),replace=True)]; yb=Y[ix]; rb=response[ix]; pb=pdelt[ix]; pbr2,_,_=multivariate_partial_r2(yb,rb,pb); c,a=coupling_stats(yb,rb,pb); boot.append([pdelt[ix[:len(ri)]].mean(),pdelt[ix[:len(ri)]].mean()-pdelt[ix[len(ri):]].mean(),pbr2,c,a])
boot=np.asarray(boot); names=["responder_platelet_change","between_group_platelet_change_difference","overlap_free_multivariate_partial_r2","overlap_free_median_gene_correlation","overlap_free_median_attenuation"]; obs=[pdelt[ri].mean(),pdelt[ri].mean()-pdelt[ni].mean(),pr2,corr0,att0]; boot_summary=[]
for j,n in enumerate(names):
    lo,hi=np.nanpercentile(boot[:,j],[2.5,97.5]); boot_summary.append({"estimand":n,"estimate":obs[j],"bootstrap_95ci_low":lo,"bootstrap_95ci_high":hi,"bootstrap_replicates":10000})
pd.DataFrame(boot_summary).to_csv(OUT/"patient_bootstrap_and_variance_attribution.tsv",sep="\t",index=False)
pd.DataFrame([{"top_genes_initial_n":208,"direct_platelet_score_overlap_n":208-len(analysis_genes),"overlap_free_gene_n":len(analysis_genes),"reduced_sse":sse0,"platelet_adjusted_sse":sse1,"multivariate_partial_r2":pr2,"lopo_reduced_prediction_sse":cv_sse0,"lopo_platelet_prediction_sse":cv_sse1,"lopo_incremental_predictive_r2":cv_incremental_r2,"nested_lopo_reduced_prediction_sse":nested_sse0,"nested_lopo_platelet_prediction_sse":nested_sse1,"nested_lopo_incremental_predictive_r2":nested_cv_r2,"nested_lopo_gene_n_min":min(nested_gene_n),"nested_lopo_gene_n_max":max(nested_gene_n),"median_gene_correlation":corr0,"median_gene_attenuation":att0}]).to_csv(OUT/"composition_variance_attribution.tsv",sep="\t",index=False)

# Complete-case month 1, no imputation
month=[]
for t1,t2 in [("pre","1mon"),("1wk","1mon")]:
    w=paired_table(meta,platelet.to_dict(),t1,t2)
    for grp in ["responder","nonresponder"]:
        x=w[w.response==grp].delta.values; month.append({"score":"platelet_abundance","contrast":f"{t2}_minus_{t1}","group":grp,"n_complete":len(x),"mean_change":x.mean() if len(x) else np.nan,"exact_test_p":exact_signflip(x) if len(x) else np.nan})
    a=w[w.response=="responder"].delta.values; n=w[w.response=="nonresponder"].delta.values; month.append({"score":"platelet_abundance","contrast":f"{t2}_minus_{t1}","group":"between_group","n_complete":len(a)+len(n),"mean_change":a.mean()-n.mean(),"exact_test_p":exact_between(a,n)})
pd.DataFrame(month).to_csv(OUT/"complete_case_month1.tsv",sep="\t",index=False)

# Young-platelet signature size and residualization sensitivity
efiles=sorted(glob.glob(str(EXT/"GSM*_HTSeqCount_Lib-DBO-*.txt.gz"))); ext=pd.concat([load_two_col_gz(f) for f in efiles],axis=1).fillna(0); elog=np.log2(ext.add(.5).div(ext.sum(axis=0).add(1),axis=1)*1e6+1); ed=pd.DataFrame({f"pair{i+1}":elog.iloc[:,2*i]-elog.iloc[:,2*i+1] for i in range(4)}); estat=pd.DataFrame({"gene":ed.index,"mean_delta":ed.mean(axis=1),"min_delta":ed.min(axis=1),"max_delta":ed.max(axis=1),"mean_logcpm":elog.mean(axis=1)}).sort_values("mean_delta",ascending=False)
sens=[]
for size in [25,50,100]:
    young=estat[(estat.min_delta>0)&(estat.mean_logcpm>1)].head(size).gene.tolist(); mature=estat[(estat.max_delta<0)&(estat.mean_logcpm>1)].tail(size).gene.tolist(); ys,yg=zscore_signature(logcpm,young); ms,mg=zscore_signature(logcpm,mature); con=ys-ms
    for method,mask in {"all_samples":np.ones(len(meta),bool),"baseline_fitted":(meta.time.values=="pre")}.items():
        rv,_=residualize(con.values,platelet.values,mask); resid=pd.Series(rv,index=con.index)
        for score_name,score in [("young_minus_mature",con),("age_residual",resid)]:
            w=paired_table(meta,score.to_dict(),"pre","1wk"); a=w[w.response=="responder"].delta.values; n=w[w.response=="nonresponder"].delta.values; sens.append({"requested_tail_size":size,"young_genes_detected":len(yg),"mature_genes_detected":len(mg),"residualization":method,"score":score_name,"responder_n":len(a),"responder_mean_change":a.mean(),"responder_exact_signflip_p":exact_signflip(a),"nonresponder_n":len(n),"nonresponder_mean_change":n.mean(),"between_exact_p":exact_between(a,n),"sample_correlation_with_platelet_score":np.corrcoef(con,platelet)[0,1]})
pd.DataFrame(sens).to_csv(OUT/"young_platelet_signature_size_sensitivity.tsv",sep="\t",index=False)

# Four-panel upgrade figure (tabular results remain primary if plotting is unavailable).
bs=pd.DataFrame(boot_summary).set_index("estimand"); ss=pd.DataFrame(sens); cc=pd.DataFrame(month)
if not HAS_MPL:
    print(pd.DataFrame(boot_summary).to_string(index=False)); print(pd.DataFrame(month).to_string(index=False)); print(pd.DataFrame(sens).to_string(index=False))
    raise SystemExit(0)
plt.rcParams.update({"font.family":"DejaVu Sans","font.size":9,"axes.spines.top":False,"axes.spines.right":False,"savefig.dpi":600}); fig,axs=plt.subplots(2,2,figsize=(10,7.2))
ax=axs[0,0]; keys=["overlap_free_multivariate_partial_r2","overlap_free_median_gene_correlation","overlap_free_median_attenuation"]; x=np.arange(3); vals=bs.loc[keys,"estimate"].values; lo=bs.loc[keys,"bootstrap_95ci_low"].values; hi=bs.loc[keys,"bootstrap_95ci_high"].values; ax.errorbar(x,vals,yerr=[vals-lo,hi-vals],fmt="o",capsize=4,color="#2B6CB0"); ax.axhline(0,color="#A0AEC0",lw=.8); ax.set_xticks(x,["Partial $R^2$","Median\ncorrelation","Median\nattenuation"]); ax.set_title("Patient-bootstrap uncertainty")
ax=axs[0,1]; q=ss[(ss.score=="age_residual")&(ss.residualization=="baseline_fitted")]; ax.plot(q.requested_tail_size,q.responder_mean_change,"o-",label="Responder",color="#2B6CB0"); ax.plot(q.requested_tail_size,q.nonresponder_mean_change,"o-",label="Nonresponder",color="#C53030"); ax.axhline(0,color="#A0AEC0",lw=.8); ax.set_xlabel("Requested genes per signature tail"); ax.set_ylabel("Week 1 - baseline residual change"); ax.set_title("Platelet-age sensitivity"); ax.legend(frameon=False)
ax=axs[1,0]; q=cc[cc.group.isin(["responder","nonresponder"])]
for i,c in enumerate(["1mon_minus_pre","1mon_minus_1wk"]):
    for j,g in enumerate(["responder","nonresponder"]):
        z=q[(q.contrast==c)&(q.group==g)].iloc[0]; ax.bar(i+(j-.5)*.32,z.mean_change,width=.32,color=["#2B6CB0","#C53030"][j])
ax.axhline(0,color="#A0AEC0",lw=.8); ax.set_xticks([0,1],["Month 1 -\nbaseline","Month 1 -\nweek 1"]); ax.set_ylabel("Platelet-score change"); ax.set_title("Complete-case month 1")
ax=axs[1,1]; ax.axis("off"); ax.text(.02,.96,"Composition-dominance framework",fontweight="bold",fontsize=11,va="top"); steps=["External lineage score","Within-patient coupling","Effect attenuation / partial $R^2$","Lineage + matched-set controls","Patient-exclusion stability","Multiplicity-corrected residual test"]
for i,s in enumerate(steps): ax.text(.08,.82-i*.13,s,va="center",bbox=dict(boxstyle="round,pad=.3",fc="#EDF2F7",ec="#718096"))
for i,a in enumerate(axs.flat): a.text(-.10,1.04,chr(65+i),transform=a.transAxes,fontweight="bold",fontsize=12)
fig.tight_layout()
for e in ["png","pdf","svg"]: fig.savefig(OUT/f"Depth_upgrade_analysis.{e}",bbox_inches="tight")
plt.close(fig)
print(pd.DataFrame(boot_summary).to_string(index=False)); print(pd.DataFrame(month).to_string(index=False)); print(pd.DataFrame(sens).to_string(index=False))
