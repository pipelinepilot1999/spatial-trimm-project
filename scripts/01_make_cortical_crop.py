"""
01_make_cortical_crop.py
------------------------
Turn the whole-brain ABC Atlas MERFISH annotation into a TrimNN-ready input for ONE
cortical patch of ONE coronal section.

WHY a spatial crop and not a random subsample:
  TrimNN builds a Delaunay graph where an edge = "these two cells are physically adjacent".
  Randomly deleting cells and then triangulating the survivors invents edges between cells
  that were never neighbours -> every motif becomes an artdefact of the deletion.
  A contiguous crop keeps EVERY cell inside a window, so local adjacency stays true.
  Bonus: a cortical patch naturally contains only ~8-12 cell classes (like TrimNN's demo),
  so we never have to arbitrarily collapse the 34 whole-brain classes.

Outputs (data/processed/):
  trimnn_input_<tag>.csv     columns X,Y,cell_type   <- fed straight to csv2gml.py
  crop_meta_<tag>.csv        cell_label,X,Y,class,subclass  <- provenance + downstream join key
"""
import pandas as pd
import numpy as np
from pathlib import Path

SRC = Path("data/raw/cell_metadata_with_cluster_annotation.csv")
OUT = Path("data/processed"); OUT.mkdir(parents=True, exist_ok=True)

SECTION = "C57BL6J-638850.35"   # chosen empirically: cleanest cortical patch, ~8 classes (demo-like)
N_CELLS = 4000                  # crop size (bounds CPU runtime; demo was 743 -> ~6 min at size-3)
MIN_CELLS_PER_CLASS = 20        # drop ultra-rare classes in the crop (can't form reliable motifs)
TAG = "sec35_cortex"
# The 4 cortical classes, used only to *locate* a cortex-dense seed for the crop centre.
CORTEX = {"01 IT-ET Glut", "02 NP-CT-L6b Glut", "06 CTX-CGE GABA", "07 CTX-MGE GABA"}

# --- load just this section ---
df = pd.read_csv(SRC, usecols=["cell_label", "brain_section_label", "x", "y", "class", "subclass"])
sec = df[df.brain_section_label == SECTION].copy()
print(f"section {SECTION}: {len(sec)} cells, {sec['class'].nunique()} classes")

# --- find a cortex-dense seed via a 2D histogram of cortical cells ---
ctx = sec[sec["class"].isin(CORTEX)]
H, xe, ye = np.histogram2d(ctx.x, ctx.y, bins=40)
i, j = np.unravel_index(np.argmax(H), H.shape)
seed_x = 0.5 * (xe[i] + xe[i + 1])
seed_y = 0.5 * (ye[j] + ye[j + 1])
print(f"cortical seed at (x={seed_x:.1f}, y={seed_y:.1f}); densest bin has {H[i,j]:.0f} cortical cells")

# --- take the N cells nearest the seed = a compact contiguous patch (keeps ALL types within it) ---
sec["d2"] = (sec.x - seed_x) ** 2 + (sec.y - seed_y) ** 2
crop = sec.nsmallest(N_CELLS, "d2").copy()
radius = np.sqrt(crop["d2"].max())

# drop ultra-rare classes (a handful of cells that only add noise + inflate the type count)
keep = crop["class"].value_counts()
keep = keep[keep >= MIN_CELLS_PER_CLASS].index
dropped = crop[~crop["class"].isin(keep)]
crop = crop[crop["class"].isin(keep)].copy()
print(f"\ncrop: {len(crop)} cells within radius {radius:.2f} of seed "
      f"(dropped {len(dropped)} cells in {dropped['class'].nunique()} rare classes)")
print(f"crop classes: {crop['class'].nunique()}")
print(crop["class"].value_counts().to_string())

# --- write TrimNN input (exact column names it expects) + provenance file ---
trimnn = crop.rename(columns={"x": "X", "y": "Y", "class": "cell_type"})[["X", "Y", "cell_type"]]
trimnn.to_csv(OUT / f"trimnn_input_{TAG}.csv", index=False)
crop.rename(columns={"x": "X", "y": "Y"})[["cell_label", "X", "Y", "class", "subclass"]] \
    .to_csv(OUT / f"crop_meta_{TAG}.csv", index=False)
print(f"\nwrote {OUT}/trimnn_input_{TAG}.csv  and  crop_meta_{TAG}.csv")
