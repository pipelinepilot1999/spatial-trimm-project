"""
03_motif_statistics.py
----------------------
Independent statistical validation of motif non-randomness. This does NOT trust TrimNN's
"predicted occurrence" (which we showed is abundance-driven, not enrichment vs a null).

>>> THIS IS THE PROJECT'S SCIENTIFIC CORE — review the method choices below and make them
>>> your own before defending them. Every design decision is called out in a NOTE comment.

Three analyses:

  (1) EXACT triangle counts + cross-check vs TrimNN's neural prediction.
      A size-3 motif IS a labeled triangle (3-clique). We can count these EXACTLY from the
      graph, so we don't need TrimNN's approximation at all for size-3 — and comparing the two
      tells us how good TrimNN's predictor is on this data.

  (2) PERMUTATION NULL for enrichment (the honest test of "spatial preference").
      NOTE (null model): we permute cell-type LABELS across the fixed set of node positions /
      fixed graph. This holds constant: the graph topology, the number of cells, and the exact
      cell-type composition. It destroys only the SPATIAL ARRANGEMENT of types. So a significant
      result means "these types are arranged non-randomly in space", controlling for abundance.
      This is the standard spatial null (cf. squidpy neighborhood enrichment).
      We report observed count, null mean/sd, z-score, one-sided empirical p (enrichment),
      and Benjamini-Hochberg FDR across all tested motifs.

  (3) FISHER'S EXACT + CRAMER'S V on pairwise neighbor (edge-level) organization.
      NOTE: Fisher assumes independent observations; graph edges are NOT independent, so we treat
      the permutation z (same null as #2, at edge level) as the trusted p-value and report Fisher
      alongside as the classical-but-approximate view. Cramer's V gives ONE global effect size for
      "how strongly do cell types organize along edges" (0 = none, 1 = perfect).
"""
import numpy as np
import pandas as pd
from pathlib import Path
from scipy.spatial import Delaunay
from scipy.stats import fisher_exact
from itertools import combinations_with_replacement
from collections import Counter, defaultdict

RNG_SEED = 0                 # NOTE: fixed seed -> reproducible permutation null
K_PERM = 2000                # NOTE: number of label permutations (empirical p resolution ~1/K)
PROC = Path("data/processed")
OUT = Path("results/statistics"); OUT.mkdir(parents=True, exist_ok=True)

# ---------- load crop + label map ----------
df = pd.read_csv(PROC / "trimnn_input_sec35_cortex.csv")
id_map = pd.read_csv(PROC / "cell_type_to_id_sec35_cortex.csv")
ct2id = dict(zip(id_map.cell_type, id_map.cell_type_id))
id2ct = {v: k.replace("01 ", "").replace("02 ", "").replace("06 ", "").replace("07 ", "")
             .replace("30 ", "").replace("31 ", "").replace("33 ", "").replace("34 ", "")
         for k, v in ct2id.items()}
labels = df.cell_type.map(ct2id).to_numpy()
n = len(labels)
K_types = len(ct2id)

# ---------- rebuild the SAME graph csv2gml builds (Delaunay, no prune) ----------
pts = df[["X", "Y"]].to_numpy()
tri = Delaunay(pts)
# adjacency as a set of undirected edges
edges = set()
for s in tri.simplices:
    for a, b in ((s[0], s[1]), (s[1], s[2]), (s[0], s[2])):
        edges.add((min(a, b), max(a, b)))
neigh = defaultdict(set)
for a, b in edges:
    neigh[a].add(b); neigh[b].add(a)
print(f"graph: {n} nodes, {len(edges)} edges")

# ---------- (1) EXACT labeled-triangle enumeration ----------
# enumerate 3-cliques: for each edge (a,b), common neighbours c>... form a triangle
triangles = []            # list of (i,j,k) node indices
for a, b in edges:
    for c in neigh[a] & neigh[b]:
        if c > b:         # canonical order a<b<c avoids counting each triangle 3x  (a<b already)
            triangles.append((a, b, c))
triangles = np.array(triangles)
print(f"exact 3-cliques (triangles) in graph: {len(triangles)}")

def motif_key(lab_triple):
    return tuple(sorted(lab_triple))                 # unordered multiset of cell-type ids

def count_labeled_triangles(lab):
    tri_labels = lab[triangles]                      # (T,3)
    c = Counter(tuple(sorted(t)) for t in tri_labels)
    return c

obs_counts = count_labeled_triangles(labels)

# ---------- (2) PERMUTATION NULL ----------
rng = np.random.default_rng(RNG_SEED)
motif_keys = sorted(obs_counts.keys())
null = {m: np.empty(K_PERM) for m in motif_keys}
perm = labels.copy()
for p in range(K_PERM):
    rng.shuffle(perm)
    c = count_labeled_triangles(perm)
    for m in motif_keys:
        null[m][p] = c.get(m, 0)

rows = []
for m in motif_keys:
    obs = obs_counts[m]
    nd = null[m]
    mu, sd = nd.mean(), nd.std()
    z = (obs - mu) / sd if sd > 0 else np.nan
    p_up = (np.sum(nd >= obs) + 1) / (K_PERM + 1)     # one-sided enrichment (add-one)
    rows.append(dict(motif="+".join(id2ct[i] for i in m), ids=m,
                     obs=obs, null_mean=round(mu, 2), z=round(z, 2), p_enrich=p_up))
res = pd.DataFrame(rows)
# Benjamini-Hochberg FDR on the enrichment p-values
res = res.sort_values("p_enrich").reset_index(drop=True)
res["fdr"] = (res["p_enrich"] * len(res) / (res.index + 1)).cummin().clip(upper=1)
res = res.sort_values("z", ascending=False)

# ---------- cross-check vs TrimNN prediction ----------
tp = pd.read_csv("results/motifs_sec35_size3/Predicted_occurrence_size3.csv")
tp.columns = [c.strip() for c in tp.columns]
import ast
tp["ids"] = tp["label"].apply(lambda s: tuple(sorted(ast.literal_eval(s))))
tp["pred"] = pd.to_numeric(tp[[c for c in tp.columns if "occ" in c.lower()][0]], errors="coerce")
tp = tp.groupby("ids", as_index=False)["pred"].sum()
merged = res.merge(tp, on="ids", how="left")
corr = merged[["obs", "pred"]].corr().iloc[0, 1]

# ---------- (3) FISHER + CRAMER'S V on edge-level organization ----------
E = np.array(sorted(edges))
et = np.sort(labels[E], axis=1)                       # (E,2) endpoint types, sorted
ct = np.zeros((K_types, K_types), dtype=int)          # symmetric edge-type contingency
for a, b in et:
    ct[a, b] += 1
    if a != b: ct[b, a] += 1
# global Cramer's V vs independence-expected
row = ct.sum(1); tot = ct.sum()
exp = np.outer(row, row) / tot
chi2 = np.nansum((ct - exp) ** 2 / np.where(exp > 0, exp, np.nan))
cramers_v = np.sqrt(chi2 / (tot * (K_types - 1)))
# per-pair Fisher (classical view; permutation z is the trusted one)
fisher_rows = []
for A, B in combinations_with_replacement(range(K_types), 2):
    ab = ct[A, B]
    a_any = ct[A].sum(); b_any = ct[B].sum()
    tbl = [[ab, a_any - ab], [b_any - ab, tot - a_any - b_any + ab]]
    try:
        orr, pf = fisher_exact(tbl)
    except Exception:
        orr, pf = np.nan, np.nan
    fisher_rows.append(dict(pair=f"{id2ct[A]}~{id2ct[B]}", ab_edges=ab, odds=round(orr, 2), p=pf))
fisher = pd.DataFrame(fisher_rows).sort_values("p")

# ---------- report + save ----------
print(f"\n[cross-check] exact vs TrimNN predicted occurrence: Pearson r = {corr:.3f}")
print(f"[global] Cramer's V (edge cell-type organization) = {cramers_v:.3f}  (0=random, 1=perfect)")
print("\n=== top enriched motifs (permutation null, FDR<0.05) ===")
sig = res[res.fdr < 0.05]
print(sig.head(12)[["motif", "obs", "null_mean", "z", "p_enrich", "fdr"]].to_string(index=False))
print(f"\n{len(sig)}/{len(res)} motifs significant at FDR<0.05")
print("\n=== strongest pairwise neighbor associations (Fisher) ===")
print(fisher.head(8).to_string(index=False))

res.to_csv(OUT / "motif_permutation_enrichment_size3.csv", index=False)
fisher.to_csv(OUT / "pairwise_neighbor_fisher.csv", index=False)
pd.DataFrame({"cramers_v": [cramers_v], "exact_vs_trimnn_pearson_r": [corr],
             "n_triangles": [len(triangles)], "K_perm": [K_PERM]}).to_csv(
             OUT / "global_summary.csv", index=False)
print(f"\nsaved -> {OUT}/")
