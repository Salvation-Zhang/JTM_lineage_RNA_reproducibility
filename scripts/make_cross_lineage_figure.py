from pathlib import Path
import numpy as np
import pandas as pd
from reportlab.pdfgen import canvas
from reportlab.lib.colors import HexColor, black, white
from reportlab.lib.pagesizes import landscape, A4

REPO=Path(__file__).resolve().parents[1]
OUT=REPO/"figures"; OUT.mkdir(parents=True,exist_ok=True)
R=REPO/"results"/"cross_lineage"
summ=pd.read_csv(R/"cross_lineage_composition_dominance_summary.tsv",sep="\t")
mix=pd.read_csv(R/"B_cell_composition_mixing_benchmark_replicates.tsv",sep="\t")
W,H=landscape(A4); c=canvas.Canvas(str(OUT/"Figure_5_cross_lineage_composition_benchmark.pdf"),pagesize=(W,H))
navy=HexColor("#183B56"); blue=HexColor("#2D7DD2"); red=HexColor("#C94C4C"); gold=HexColor("#E6A23C"); gray=HexColor("#667085"); pale=HexColor("#F4F7FA"); green=HexColor("#2E8B57")

def txt(x,y,s,size=8,color=black,bold=False,anchor="left"):
    c.setFillColor(color); c.setFont("Helvetica-Bold" if bold else "Helvetica",size)
    if anchor=="center": c.drawCentredString(x,y,s)
    elif anchor=="right": c.drawRightString(x,y,s)
    else: c.drawString(x,y,s)
def panel(x,y,w,h,label,title):
    c.setFillColor(pale); c.roundRect(x,y,w,h,6,fill=1,stroke=0); txt(x+8,y+h-15,label,12,navy,True); txt(x+28,y+h-14,title,9,navy,True)
def arrow(x1,y1,x2,y2,color=gray):
    c.setStrokeColor(color); c.setLineWidth(1.2); c.line(x1,y1,x2,y2); c.line(x2,y2,x2-5,y2+3); c.line(x2,y2,x2-5,y2-3)

margin=28; gap=12; pw=(W-2*margin-gap)/2; ph=(H-2*margin-gap)/2
panel(margin,H-margin-ph,pw,ph,"A","Intervention-targeted lineage perturbations")
panel(margin+pw+gap,H-margin-ph,pw,ph,"B","Target-lineage RNA-score changes")
panel(margin,margin,pw,ph,"C","Variance attribution and held-out reconstruction")
panel(margin+pw+gap,margin,pw,ph,"D","Composition-only benchmark (100 replicates)")

# A
x=margin+22; y=H-margin-58
items=[("Eltrombopag","ITP","Platelet restoration",blue),("Erythropoietin","Healthy volunteers","Erythroid restoration",red),("Rituximab","Type 1 diabetes","B-cell depletion",gold)]
for i,(drug,cohort,effect,col) in enumerate(items):
    yy=y-i*55; c.setFillColor(white); c.roundRect(x,yy-28,pw-44,40,5,fill=1,stroke=0); txt(x+10,yy,drug,9,col,True); txt(x+105,yy,cohort,8,gray); arrow(x+205,yy+2,x+245,yy+2,col); txt(x+255,yy,effect,8,col,True)
txt(x,y-178,"Externally specified scores; no outcome-driven marker selection",7,gray)

# B
bx=margin+pw+gap+50; by=H-margin-ph+38; bh=ph-80
c.setStrokeColor(gray); c.line(bx,by,bx,by+bh); c.line(bx,by+bh*.48,bx+pw-90,by+bh*.48)
vals=[("Platelet",0.886,blue,"P=0.00195"),("Erythroid",1.187,red,"P=0.00195"),("B cell",-1.912,gold,"P=0.00001")]
scale=bh/4.8
for i,(lab,v,col,pv) in enumerate(vals):
    xx=bx+50+i*85; zero=by+bh*.48; top=zero+v*scale; c.setFillColor(col); c.rect(xx-14,min(zero,top),28,abs(top-zero),fill=1,stroke=0); txt(xx,by-14,lab,7,gray,anchor="center"); txt(xx,top+(5 if v>0 else -12),f"{v:+.3f}",7,col,True,anchor="center"); txt(xx,by-27,pv,6,gray,anchor="center")
txt(bx-8,by+bh-5,"standardized change",7,gray)

# C
cx=margin+55; cy=margin+38; cw=pw-90; ch=ph-78
labels=["Eltrombopag / platelet","EPO / erythroid","Rituximab / B cell"]
partial=[.771,float(summ.loc[summ.dataset=="GSE186294","partial_r2"].iloc[0]),float(summ.loc[summ.dataset=="GSE112594","partial_r2"].iloc[0])]
pred=[.702,float(summ.loc[summ.dataset=="GSE186294","nested_lopo_predictive_r2"].iloc[0]),float(summ.loc[summ.dataset=="GSE112594","nested_lopo_predictive_r2"].iloc[0])]
for i,lab in enumerate(labels):
    yy=cy+ch-25-i*47; txt(cx-8,yy+4,lab,7,gray); c.setFillColor(blue); c.rect(cx+105,yy,partial[i]*170,10,fill=1,stroke=0); c.setFillColor(green); c.rect(cx+105,yy-13,pred[i]*170,8,fill=1,stroke=0); txt(cx+280,yy+2,f"{partial[i]:.3f}",7,blue,True); txt(cx+280,yy-13,f"{pred[i]:.3f}",7,green,True)
txt(cx+105,cy-4,"0",7,gray); txt(cx+275,cy-4,"1.0",7,gray,anchor="right"); txt(cx+105,cy-18,"Partial R2",7,blue,True); txt(cx+180,cy-18,"Nested LOPO held-out R2",7,green,True)

# D
dx=margin+pw+gap+45; dy=margin+42; dh=ph-82; zero=dy
cols=[("Unadjusted","DE_FDR05_unadjusted",red),("B-score adjusted","DE_FDR05_B_score_adjusted",blue)]
maxv=4500
for j,(lab,col,color) in enumerate(cols):
    xx=dx+80+j*145; vals=mix[col].values
    for k,v in enumerate(vals):
        jitter=((k*37)%23-11)*0.75; yy=zero+min(v,maxv)/maxv*dh
        c.setFillColor(HexColor("#D9E1E8") if j==0 else HexColor("#BFD7F2")); c.circle(xx+jitter,yy,1.4,fill=1,stroke=0)
    med=np.median(vals); med_label=int(np.floor(med+0.5)); c.setStrokeColor(color); c.setLineWidth(2.5); c.line(xx-24,zero+med/maxv*dh,xx+24,zero+med/maxv*dh); txt(xx,dy-16,lab,7,color,True,anchor="center"); txt(xx,dy-29,f"median {med_label:,}",7,color,anchor="center")
c.setStrokeColor(gray); c.line(dx,dy,dx,dy+dh); c.line(dx,dy,dx+330,dy)
for v in [0,1000,2000,3000,4000]: txt(dx-7,dy+v/maxv*dh-2,str(v),6,gray,anchor="right")
txt(dx-25,dy+dh+4,"FDR-significant genes",7,gray)
txt(W-30,12,"Scores are lineage-derived RNA-abundance proxies, not measured cell counts.",6,gray,anchor="right")
c.showPage(); c.save(); print(OUT/"Figure_5_cross_lineage_composition_benchmark.pdf")
