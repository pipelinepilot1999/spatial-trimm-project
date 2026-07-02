"""
09_figures.py — figure-ready PNGs for RESULTS.md.
Palette: validated dataviz reference (blue sequential ordinal ramp; blue<->red diverging).
Principles: thin marks, direct value labels, recessive axes, text in ink tokens (not series color).
"""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
import pandas as pd, numpy as np
from scipy.spatial import Delaunay
from collections import defaultdict
from pathlib import Path

OUT = Path("results/figures"); OUT.mkdir(parents=True, exist_ok=True)
SURF="#fcfcfb"; INK="#0b0b0b"; INK2="#52514e"; GRIDC="#e7e6e2"
BLUE={250:"#86b6ef",450:"#2a78d6",650:"#104281"}
plt.rcParams.update({"figure.facecolor":SURF,"axes.facecolor":SURF,"savefig.facecolor":SURF,
    "font.size":11,"text.color":INK,"axes.labelcolor":INK2,"xtick.color":INK2,"ytick.color":INK2,
    "axes.edgecolor":GRIDC,"font.family":"DejaVu Sans"})
def clean(ax):
    for s in ("top","right"): ax.spines[s].set_visible(False)
    ax.spines["left"].set_color(GRIDC); ax.spines["bottom"].set_color(GRIDC)

# ---------- Fig 1: DEG collapse (headline) ----------
d = pd.read_csv("results/de/deg_summary.csv")
labels = ["Naive\n(confounded +\npseudoreplicated)","Confound-\ncontrolled","Pseudobulk\n(honest)"]
vals = d.degs.tolist()
fig,ax=plt.subplots(figsize=(6.4,4.2))
cols=[BLUE[250],BLUE[450],BLUE[650]]
bars=ax.bar(range(3),vals,color=cols,width=0.62,zorder=3)
for i,v in enumerate(vals): ax.text(i,v+4,str(v),ha="center",va="bottom",fontweight="bold",color=INK,fontsize=13)
# annotate the two fixes
ax.annotate("cell-type confound\nremoved  −78%",xy=(0.5,135),ha="center",color=INK2,fontsize=9)
ax.annotate("pseudoreplication\nremoved  −65%",xy=(1.5,40),ha="center",color=INK2,fontsize=9)
ax.set_xticks(range(3)); ax.set_xticklabels(labels,fontsize=9)
ax.set_ylabel("Differentially expressed genes (padj < 0.05)")
ax.set_ylim(0,250); ax.set_yticks([0,50,100,150,200])
ax.set_title("Each statistical fix collapses the DEG count",fontweight="bold",loc="left",color=INK)
ax.grid(axis="y",color=GRIDC,lw=0.8,zorder=0); clean(ax)
plt.tight_layout(); plt.savefig(OUT/"fig1_deg_collapse.png",dpi=150); plt.close()

# ---------- Fig 2: motif enrichment (permutation z) ----------
import re
_sn={0:"IT-ET",1:"L6b",2:"CGE",3:"MGE",4:"Astro",5:"OPC-Oligo",6:"Vascular",7:"Immune"}
m = pd.read_csv("results/statistics/motif_permutation_enrichment_size3.csv")
def _ids(s):
    got=re.findall(r"int\d*\((\d+)\)",s)        # numbers inside np.int64(...)
    if not got: got=re.findall(r"\b([0-7])\b",s)  # fallback: bare single-digit ids
    return [int(i) for i in got]
m["name"]=m["ids"].apply(lambda s:" + ".join(_sn[i] for i in _ids(s)))
top = m.sort_values("z",ascending=False).head(8).iloc[::-1]
fig,ax=plt.subplots(figsize=(9.2,4.6))
ax.barh(range(len(top)),top.z,color=BLUE[450],height=0.66,zorder=3)
for i,z in enumerate(top.z):
    ax.text(z+2,i,f"z={z:.0f}",va="center",color=INK2,fontsize=9)
ax.set_yticks(range(len(top))); ax.set_yticklabels(top.name,fontsize=9.5)
ax.set_xlabel("Enrichment vs abundance-preserving null  (permutation z-score, K=2000)")
ax.set_title("Spatially over-represented motifs, abundance-controlled",
             fontweight="bold",loc="left",color=INK,fontsize=12,pad=10)
ax.set_xlim(0,max(top.z)*1.18); ax.grid(axis="x",color=GRIDC,lw=0.8,zorder=0); clean(ax)
plt.tight_layout(); plt.savefig(OUT/"fig2_motif_enrichment.png",dpi=150,bbox_inches="tight"); plt.close()

# ---------- Fig 3: pairwise neighbor log2 odds heatmap (diverging) ----------
df=pd.read_csv("data/processed/trimnn_input_sec35_cortex.csv")
idm=pd.read_csv("data/processed/cell_type_to_id_sec35_cortex.csv")
ct2id=dict(zip(idm.cell_type,idm.cell_type_id))
short={0:"IT-ET",1:"L6b",2:"CGE",3:"MGE",4:"Astro",5:"OPC-Oligo",6:"Vascular",7:"Immune"}
lab=df.cell_type.map(ct2id).to_numpy(); K=8
pts=df[["X","Y"]].to_numpy(); tri=Delaunay(pts); edges=set()
for s in tri.simplices:
    for a,b in ((s[0],s[1]),(s[1],s[2]),(s[0],s[2])): edges.add((min(a,b),max(a,b)))
ct=np.zeros((K,K))
for a,b in edges:
    x,y=sorted((lab[a],lab[b])); ct[x,y]+=1; ct[y,x]+=1
row=ct.sum(1); tot=ct.sum(); exp=np.outer(row,row)/tot
with np.errstate(divide="ignore",invalid="ignore"):
    lo=np.log2(ct/exp)
diverge=LinearSegmentedColormap.from_list("bwr_v",["#104281","#2a78d6","#f0efec","#e34948","#8f1f1f"])
fig,ax=plt.subplots(figsize=(6.6,5.6))
vmax=np.nanmax(np.abs(lo[np.isfinite(lo)]))
im=ax.imshow(lo,cmap=diverge,vmin=-vmax,vmax=vmax)
ax.set_xticks(range(K)); ax.set_yticks(range(K))
ax.set_xticklabels(short.values(),rotation=45,ha="right",fontsize=9); ax.set_yticklabels(short.values(),fontsize=9)
for i in range(K):
    for j in range(K):
        if np.isfinite(lo[i,j]):
            ax.text(j,i,f"{lo[i,j]:.1f}",ha="center",va="center",fontsize=7,
                    color=INK if abs(lo[i,j])<vmax*0.6 else "#ffffff")
cb=fig.colorbar(im,ax=ax,fraction=0.046,pad=0.04); cb.set_label("log2(observed / expected) adjacency",color=INK2,fontsize=9)
cb.outline.set_edgecolor(GRIDC)
ax.set_title("Cell-type neighbor organization (red=co-cluster, blue=avoid)",fontweight="bold",loc="left",color=INK,fontsize=11)
plt.tight_layout(); plt.savefig(OUT/"fig3_neighbor_logodds.png",dpi=150); plt.close()

print("wrote:", *[p.name for p in sorted(OUT.glob('*.png'))])
