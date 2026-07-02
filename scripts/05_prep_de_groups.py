"""
05_prep_de_groups.py
--------------------
Prepare the grouping columns the DESeq2 script needs, for the OPC-Oligo focal contrast
(cells of type OPC-Oligo that participate in an OPC-Oligo-OPC-Oligo-OPC-Oligo triangle vs not).

Emits data/processed/de_groups_sec35.csv with, per cell:
  cell_label, cell_type, X, Y
  is_focal        : cell_type == OPC-Oligo
  in_focal_motif  : focal cell is a vertex of >=1 homotypic OPC-Oligo triangle
  naive_group     : "motif" if cell is in any OPC-Oligo triangle else "other"  (confounded contrast)
  tile            : spatial grid tile id (for pseudobulk pseudo-replicates)
"""
import pandas as pd, numpy as np
from scipy.spatial import Delaunay
from collections import defaultdict
from pathlib import Path

PROC = Path("data/processed")
FOCAL = "31 OPC-Oligo"
GRID = 4   # 4x4 spatial tiles -> pseudobulk replicates

df = pd.read_csv(PROC / "trimnn_input_sec35_cortex.csv")
meta = pd.read_csv(PROC / "crop_meta_sec35_cortex.csv")
idm = pd.read_csv(PROC / "cell_type_to_id_sec35_cortex.csv")
ct2id = dict(zip(idm.cell_type, idm.cell_type_id))
lab = df.cell_type.map(ct2id).to_numpy()
Tid = ct2id[FOCAL]

# rebuild graph + triangles (same as everywhere)
pts = df[["X", "Y"]].to_numpy(); tri = Delaunay(pts)
edges = set()
for s in tri.simplices:
    for a, b in ((s[0], s[1]), (s[1], s[2]), (s[0], s[2])): edges.add((min(a, b), max(a, b)))
nb = defaultdict(set)
for a, b in edges: nb[a].add(b); nb[b].add(a)
tris = [(a, b, c) for a, b in edges for c in (nb[a] & nb[b]) if c > b]
tris = np.array(tris)

# cells in >=1 homotypic focal triangle
homo = tris[(lab[tris] == Tid).all(axis=1)]
in_motif_nodes = set(homo.ravel())

# spatial tiles (equal-count quantile bins in X and Y so tiles hold comparable cell numbers)
xb = pd.qcut(df.X, GRID, labels=False, duplicates="drop")
yb = pd.qcut(df.Y, GRID, labels=False, duplicates="drop")
tile = xb.astype(str) + "_" + yb.astype(str)

out = pd.DataFrame({
    "cell_label": meta.cell_label.values,
    "cell_type": df.cell_type.values,
    "X": df.X.values, "Y": df.Y.values,
    "is_focal": lab == Tid,
    "in_focal_motif": [i in in_motif_nodes and lab[i] == Tid for i in range(len(df))],
    "naive_group": ["motif" if i in in_motif_nodes else "other" for i in range(len(df))],
    "tile": tile.values,
})
out.to_csv(PROC / "de_groups_sec35.csv", index=False)

# report pseudobulk viability: focal cells per tile x group
foc = out[out.is_focal]
pv = foc.groupby(["tile", "in_focal_motif"]).size().unstack(fill_value=0)
pv.columns = ["out_motif", "in_motif"]
print("OPC-Oligo cells per tile (pseudobulk replicates):")
print(pv.to_string())
usable = pv[(pv.out_motif >= 8) & (pv.in_motif >= 8)]
print(f"\ntiles usable as paired pseudo-replicates (>=8 cells each group): {len(usable)}")
print(f"focal total: {foc.shape[0]}  in_motif={foc.in_focal_motif.sum()}  out={(~foc.in_focal_motif).sum()}")
