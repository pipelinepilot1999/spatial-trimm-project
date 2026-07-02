"""
08_communication_squidpy.py
---------------------------
Cell-cell communication (ligand-receptor) analysis for the section-.35 cortical crop.

This is the RUNNABLE substitute for the deferred CellChat stage (08_cellchat.R): CellChat installs
only from GitHub and the environment guardrail blocked it. squidpy.gr.ligrec is a CellPhoneDB-style
permutation test (are a ligand in cluster A and its receptor in cluster B co-expressed more than a
label-shuffled null?) and installs cleanly from PyPI.

CAVEATS (carry into the write-up):
  - LR communication inferred from CO-EXPRESSION = a PREDICTION, not proof of signaling.
  - The LR database (via omnipath) is HUMAN gene symbols; our data is mouse. We map mouse->human by
    the 1:1-ortholog uppercase convention (Prkcq->PRKCQ). This catches most 1:1 orthologs but misses
    renamed/expanded families -> under-counts.
  - Only LR pairs whose BOTH partners are on the 550-gene MERFISH panel can ever be detected.
"""
import scanpy as sc, anndata as ad, pandas as pd, numpy as np, squidpy as sq
from pathlib import Path

PROC = Path("data/processed"); OUT = Path("results/communication"); OUT.mkdir(parents=True, exist_ok=True)

counts = pd.read_csv(PROC / "de_counts_sec35.csv", index_col=0)      # genes(ensembl) x cells
genes  = pd.read_csv(PROC / "de_genes_sec35.csv")
col    = pd.read_csv(PROC / "de_coldata_sec35.csv").set_index("cell_label")

# build AnnData: cells x genes, var = UPPERCASED mouse symbol (human-ortholog proxy)
sym = genes.set_index("ensembl_id").symbol.reindex(counts.index)
keep = sym.notna() & ~sym.str.upper().duplicated()
X = counts.loc[keep].T                                                # cells x genes
X.columns = sym[keep].str.upper().values
a = ad.AnnData(X.values.astype(float),
               obs=col.loc[X.index, ["cell_type"]].copy(),
               var=pd.DataFrame(index=X.columns))
a.obs["cell_type"] = a.obs["cell_type"].astype("category")

# standard normalization
sc.pp.normalize_total(a, target_sum=1e4); sc.pp.log1p(a)
print(f"AnnData: {a.n_obs} cells x {a.n_vars} genes; {a.obs.cell_type.nunique()} cell types")

res = sq.gr.ligrec(
    a, cluster_key="cell_type", use_raw=False,
    n_perms=1000, threshold=0.05, seed=0, copy=True,
    interactions_params={"resources": "CellPhoneDB"},
)
pvals = res["pvalues"]; means = res["means"]
n_tested = pvals.notna().any(axis=1).sum()
print(f"\nLR pairs testable on this panel (both partners present): {n_tested} of {len(pvals)}")
print(f"min p-value observed: {np.nanmin(pvals.values):.4f}")

records = []
for (lig, rec), row in pvals.iterrows():
    for (src, tgt), p in row.items():
        if pd.notna(p) and p < 0.05:
            records.append(dict(ligand=lig, receptor=rec, source=src, target=tgt,
                                mean=means.loc[(lig, rec), (src, tgt)], pval=p))
out = pd.DataFrame(records)
if len(out):
    out = out.sort_values(["pval", "mean"], ascending=[True, False])
out.to_csv(OUT / "significant_LR_pairs_squidpy.csv", index=False)
# also save all testable pairs' best result for transparency
tested = pvals[pvals.notna().any(axis=1)]
tested.to_csv(OUT / "all_testable_LR_pvalues.csv")

print(f"\nsignificant LR interactions (p<0.05): {len(out)}")
if len(out):
    print(f"unique LR pairs: {out[['ligand','receptor']].drop_duplicates().shape[0]}")
    print("\ntop LR interactions:")
    print(out.head(15).to_string(index=False))
else:
    print("-> no LR interactions reach p<0.05. Consistent with a 550-gene panel + ortholog-mapping")
    print("   loss severely truncating the detectable LR space (documented limitation).")
print(f"\nsaved -> {OUT}/")
