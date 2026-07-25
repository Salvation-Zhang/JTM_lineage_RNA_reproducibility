from pathlib import Path
import json
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch

DATA = Path(r"F:\TISTA_data\data\processed")
OUT = Path(r"D:\TISTA_JTH_manuscript\figures")
OUT.mkdir(parents=True, exist_ok=True)

plt.rcParams.update({
    "font.family": "DejaVu Sans", "font.size": 9, "axes.titlesize": 11,
    "axes.labelsize": 9, "figure.dpi": 150, "savefig.dpi": 600,
    "axes.spines.top": False, "axes.spines.right": False,
})
BLUE, RED, GREY, GOLD, GREEN = "#2B6CB0", "#C53030", "#718096", "#D69E2E", "#2F855A"

def save(fig, name):
    fig.savefig(OUT / f"{name}.png", bbox_inches="tight")
    fig.savefig(OUT / f"{name}.pdf", bbox_inches="tight")
    fig.savefig(OUT / f"{name}.svg", bbox_inches="tight")
    plt.close(fig)

def panel(ax, label):
    ax.text(-0.12, 1.06, label, transform=ax.transAxes, fontweight="bold", fontsize=12)

def figure1():
    fig, ax = plt.subplots(figsize=(10, 5.8)); ax.axis("off")
    boxes = [
        (0.04, .68, .24, .20, "Primary cohort\nGSE112278\n17 patients, 46 samples", BLUE),
        (.38, .70, .24, .16, "External platelet references\nGSE302674 + GSE262073", GREEN),
        (.72, .68, .24, .20, "Locked 12-gene\nplatelet RNA score", GOLD),
        (.08, .31, .24, .20, "Patient-level longitudinal\nand paired models", BLUE),
        (.38, .31, .24, .20, "Specificity + 5,000 matched\npermutations + LOPO", RED),
        (.68, .31, .28, .20, "Orthogonal audit\npurified T cells + BM\ndonor pseudobulk", GREY),
    ]
    for x,y,w,h,t,c in boxes:
        ax.add_patch(FancyBboxPatch((x,y),w,h,boxstyle="round,pad=0.015",fc=c,ec="none",alpha=.12))
        ax.text(x+w/2,y+h/2,t,ha="center",va="center",fontweight="bold",color=c)
    arrows = [((.62,.78),(.72,.78)),((.16,.68),(.18,.51)),((.84,.68),(.52,.51)),((.32,.41),(.38,.41)),((.62,.41),(.68,.41))]
    for a,b in arrows: ax.add_patch(FancyArrowPatch(a,b,arrowstyle="-|>",mutation_scale=14,color="#4A5568"))
    ax.text(.5,.12,"Central question: does the whole-blood response reflect platelet RNA recovery or composition-resistant reprogramming?",ha="center",fontsize=11,fontweight="bold")
    ax.text(.02,.96,"A",fontweight="bold",fontsize=13)
    save(fig,"Figure_1_study_design")

def figure2():
    scores=pd.read_csv(DATA/"GSE112278_external_validated_platelet_scores.tsv",sep="\t")
    order={"pre":0,"1wk":1,"1mon":2}; labels=["Baseline","Week 1","Month 1"]
    fig,axs=plt.subplots(1,3,figsize=(11,3.8),gridspec_kw={"width_ratios":[1.25,1.25,.9]})
    for ax,(grp,c) in zip(axs[:2],[("responder",BLUE),("nonresponder",RED)]):
        d=scores[scores.response==grp]
        for _,g in d.groupby("patient"):
            g=g.assign(o=g.time.map(order)).sort_values("o")
            ax.plot(g.o,g.external_validated_platelet_score,"-o",color=c,alpha=.55,lw=1,ms=3)
        means=d.assign(o=d.time.map(order)).groupby("o").external_validated_platelet_score.mean()
        ax.plot(means.index,means.values,"-o",color="black",lw=2.5,ms=5,label="Mean")
        ax.axhline(0,color="#CBD5E0",lw=.8); ax.set_xticks([0,1,2],labels); ax.set_title(grp.capitalize())
        ax.set_ylabel("Platelet RNA score (standardized)")
    wide=scores[scores.time.isin(["pre","1wk"])].pivot(index=["patient","response"],columns="time",values="external_validated_platelet_score").dropna()
    wide["delta"]=wide["1wk"]-wide["pre"]
    for i,(grp,c) in enumerate([("responder",BLUE),("nonresponder",RED)]):
        v=wide.xs(grp,level="response").delta.values
        x=np.full(len(v),i)+np.linspace(-.08,.08,len(v)); axs[2].scatter(x,v,color=c,s=28,zorder=3)
        axs[2].errorbar(i,v.mean(),yerr=v.std(ddof=1)/np.sqrt(len(v)),fmt="D",color="black",capsize=4)
    axs[2].axhline(0,color="#A0AEC0",lw=.8); axs[2].set_xticks([0,1],["Responder","Nonresponder"],rotation=15)
    axs[2].set_ylabel("Week 1 − baseline change"); axs[2].set_title("Paired change")
    axs[2].text(.98,.04,"Between-group exact P=0.0381",transform=axs[2].transAxes,ha="right",va="bottom",fontsize=8)
    for i,a in enumerate(axs): panel(a,chr(65+i))
    fig.tight_layout(); save(fig,"Figure_2_platelet_score_trajectories")

def figure3():
    d=pd.read_csv(DATA/"GSE112278_week1_composition_adjusted_gene_audit.tsv",sep="\t")
    top=d.nlargest(208,"responder_mean_delta")
    perm=json.loads((DATA/"GSE112278_expression_matched_gene_set_permutation.json").read_text())
    fig,axs=plt.subplots(1,2,figsize=(9,4))
    axs[0].scatter(top.delta_score_correlation,top.attenuation_fraction,s=12,alpha=.55,color=BLUE,edgecolor="none")
    axs[0].axvline(perm["observed"]["corr"],color=RED,ls="--",lw=1); axs[0].axhline(perm["observed"]["att"],color=RED,ls="--",lw=1)
    axs[0].set(xlabel="Correlation with platelet-score change",ylabel="Response-effect attenuation",title="Top 208 responder-induced genes")
    vals=[perm["null_corr_median"],perm["observed"]["corr"],perm["null_attenuation_median"],perm["observed"]["att"]]
    axs[1].bar([0,1,3,4],vals,color=[GREY,BLUE,GREY,BLUE],width=.75)
    axs[1].set_xticks([.5,3.5],["Median correlation","Median attenuation"]); axs[1].set_ylabel("Statistic")
    axs[1].set_title("Observed vs expression-matched null")
    axs[1].text(.5,.94,"Empirical P=0.00020 for both",transform=axs[1].transAxes,ha="center")
    for i,a in enumerate(axs): panel(a,chr(65+i))
    fig.tight_layout(); save(fig,"Figure_3_coupling_attenuation")

def figure4():
    d=pd.read_csv(DATA/"GSE112278_cell_signature_negative_controls.tsv",sep="\t")
    names=d.signature.replace({"T_cell":"T cell","B_cell":"B cell"}).str.capitalize()
    y=np.arange(len(d)); fig,ax=plt.subplots(figsize=(7,4.5))
    ax.errorbar(d.responder_mean_delta,y-.13,xerr=d.responder_sd/np.sqrt(d.responder_n),fmt="o",color=BLUE,label="Responder")
    ax.errorbar(d.nonresponder_mean_delta,y+.13,xerr=d.nonresponder_sd/np.sqrt(d.nonresponder_n),fmt="o",color=RED,label="Nonresponder")
    ax.axvline(0,color="#A0AEC0",lw=.8); ax.set_yticks(y,names); ax.invert_yaxis(); ax.set_xlabel("Week 1 − baseline score change (mean ± SE)")
    ax.legend(frameon=False,loc="lower right"); ax.set_title("Lineage specificity")
    ax.text(.98,.94,"Platelet between-group exact P=0.0381",transform=ax.transAxes,ha="right",va="top",color=BLUE,fontsize=8)
    panel(ax,"A"); fig.tight_layout(); save(fig,"Figure_4_lineage_specificity")

def figure5():
    d=pd.read_csv(DATA/"GSE112278_LOPO_composition_stability.tsv",sep="\t")
    cols=[("responder_score_change","Responder score change"),("top208_median_correlation","Median correlation"),("top208_median_attenuation","Median attenuation")]
    fig,axs=plt.subplots(1,3,figsize=(10,4),sharex=True); x=np.arange(len(d))
    for i,(col,title) in enumerate(cols):
        axs[i].plot(x,d[col],"o-",ms=3,lw=1,color=[BLUE,GREEN,GOLD][i]); axs[i].set_title(title); axs[i].set_xticks(x,d.held_out,rotation=90,fontsize=6)
        axs[i].axhline(0,color="#CBD5E0",lw=.8); panel(axs[i],chr(65+i))
    axs[0].set_ylabel("Leave-one-patient-out estimate")
    fig.tight_layout(); save(fig,"Figure_5_LOPO_robustness")

def figure6():
    g431=pd.read_csv(DATA/"GSE43177_purified_Tcell_immunometabolic_validation.tsv",sep="\t")
    g469=pd.read_csv(DATA/"GSE46922_purified_Tcell_chronicity_immunometabolic_validation.tsv",sep="\t")
    bm=pd.read_csv(DATA/"GSE196676_HSPC_pseudobulk_locked_pathways.tsv",sep="\t")
    comp=pd.read_csv(DATA/"GSE196676_marker_compartment_locked_pathways.tsv",sep="\t")
    rid="R-HSA-77289"
    effects=[
        g431.loc[g431.id==rid,"mean_patient"].iloc[0]-g431.loc[g431.id==rid,"mean_control"].iloc[0],
        g469.loc[g469.id==rid,"mean_chronic"].iloc[0]-g469.loc[g469.id==rid,"mean_new"].iloc[0],
        bm.loc[bm.reactome_id==rid,"delta_ITP_minus_HC"].iloc[0],
    ]
    labels=["GSE43177\nITP − control","GSE46922\nchronic − new","GSE196676\nITP − control"]
    pvals=[g431.loc[g431.id==rid,"p"].iloc[0],g469.loc[g469.id==rid,"p"].iloc[0],bm.loc[bm.reactome_id==rid,"exact_two_sided_p"].iloc[0]]
    fig,axs=plt.subplots(1,2,figsize=(9.5,4.3))
    axs[0].barh(np.arange(3),effects,color=[GREY,RED,BLUE]); axs[0].set_yticks(np.arange(3),labels); axs[0].axvline(0,color="black",lw=.8); axs[0].invert_yaxis()
    axs[0].set_xlabel("FAO pathway-score difference"); axs[0].set_title("Cross-dataset FAO effects")
    for i,(v,p) in enumerate(zip(effects,pvals)):
        axs[0].text(.98,i,f"P={p:.3g}",transform=axs[0].get_yaxis_transform(),va="center",ha="right",fontsize=8)
    f=comp[comp.reactome_id==rid].copy(); y=np.arange(len(f)); axs[1].barh(y,f.delta_ITP_minus_HC,color=BLUE,alpha=.75); axs[1].set_yticks(y,f.compartment); axs[1].invert_yaxis(); axs[1].axvline(0,color="black",lw=.8)
    axs[1].set_xlabel("ITP − control FAO score"); axs[1].set_title("Bone-marrow donor pseudobulk")
    axs[1].text(.98,.03,"FAO cross-compartment FDR=0.743",transform=axs[1].transAxes,ha="right",fontsize=8)
    for i,a in enumerate(axs): panel(a,chr(65+i))
    fig.tight_layout(); save(fig,"Figure_6_immunometabolic_audit")

if __name__ == "__main__":
    figure1(); figure2(); figure3(); figure4(); figure5(); figure6()
    print(f"Wrote figures to {OUT}")
