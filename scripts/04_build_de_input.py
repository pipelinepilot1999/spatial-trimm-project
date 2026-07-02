"""
04_build_de_input.py
--------------------
Plumbing: subset the section-.35 expression h5ad to our crop cells and emit files DESeq2 (R) can read.
Does NOT decide the DE contrast or pseudobulk scheme — those are analytical choices made separately.

Outputs (data/processed/):
  de_counts_sec35.csv    raw integer counts, GENES x CELLS (DESeq2 orientation)
  de_coldata_sec35.csv   per-cell metadata: cell_label, cell_type (class), X, Y
  de_genes_sec35.csv     gene table: ensembl_id (+ symbol if available)  -> panel = enrichment background
"""
import anndata as ad
import pandas as pd
import numpy as np
from pathlib import Path

PROC = Path("data/processed")
a = ad.read_h5ad("data/raw/C57BL6J-638850.35-raw.h5ad")
crop = pd.read_csv(PROC / "crop_meta_sec35_cortex.csv").set_index("cell_label")

# keep crop cells present in the matrix, in a stable order
cells = [c for c in crop.index if c in set(a.obs_names)]
a = a[cells].copy()
crop = crop.loc[cells]
print(f"subset: {a.shape[0]} cells x {a.shape[1]} genes")

# counts: genes x cells (DESeq2 wants features as rows, samples as columns)
X = a.X
X = np.asarray(X.todense()) if hasattr(X, "todense") else np.asarray(X)
counts = pd.DataFrame(X.T.astype(int), index=a.var_names, columns=cells)
assert (counts.values >= 0).all() and np.allclose(counts.values, counts.values.astype(int))
counts.to_csv(PROC / "de_counts_sec35.csv")

# coldata
crop.reset_index()[["cell_label", "class", "X", "Y"]] \
    .rename(columns={"class": "cell_type"}) \
    .to_csv(PROC / "de_coldata_sec35.csv", index=False)

# gene table (symbol if the h5ad carries one)
genes = pd.DataFrame({"ensembl_id": a.var_names})
for col in ("gene_symbol", "gene_name", "symbol", "name"):
    if col in a.var.columns:
        genes["symbol"] = a.var[col].values
        break
genes.to_csv(PROC / "de_genes_sec35.csv", index=False)

print("wrote de_counts_sec35.csv (%d x %d), de_coldata_sec35.csv, de_genes_sec35.csv"
      % counts.shape)
print("var columns available:", list(a.var.columns))
