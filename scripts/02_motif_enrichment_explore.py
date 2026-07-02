"""
02_motif_enrichment_explore.py
------------------------------
Quick exploratory check that motivated the rigorous stats in 03. Shows that TrimNN's raw
`predicted_occurrence` ranking is ABUNDANCE-driven, and that a first-order abundance correction
(observed / multinomial-expected) already inverts the ranking toward real spatial preference.

This is a FIRST-ORDER approximation only (assumes independence). The trustworthy version is the
permutation null in 03_motif_statistics.py — use THAT for any claim. This script exists to make the
"abundance artifact" finding reproducible and to explain why 03 is needed.

Reads: results/motifs_sec35_size3/Predicted_occurrence_size3.csv, data/processed/crop_meta + id map
Writes: results/statistics/abundance_vs_enrichment_size3.csv
"""
import pandas as pd, numpy as np, ast
from math import factorial
from collections import Counter
from pathlib import Path

OUT = Path("results/statistics"); OUT.mkdir(parents=True, exist_ok=True)
id2ct = {0:"IT-ET Glut",1:"NP-CT-L6b Glut",2:"CTX-CGE GABA",3:"CTX-MGE GABA",
         4:"Astro",5:"OPC-Oligo",6:"Vascular",7:"Immune"}
short = {"01 IT-ET Glut":0,"02 NP-CT-L6b Glut":1,"06 CTX-CGE GABA":2,"07 CTX-MGE GABA":3,
         "30 Astro-Epen":4,"31 OPC-Oligo":5,"33 Vascular":6,"34 Immune":7}

meta = pd.read_csv("data/processed/crop_meta_sec35_cortex.csv")
p = (meta["class"].map(short).value_counts().sort_index() / len(meta)).to_dict()

df = pd.read_csv("results/motifs_sec35_size3/Predicted_occurrence_size3.csv")
df.columns = [c.strip() for c in df.columns]
occ = [c for c in df.columns if "occ" in c.lower()][0]
df["labels"] = df["label"].apply(ast.literal_eval)
df["occ"] = pd.to_numeric(df[occ], errors="coerce").fillna(0)
mult = lambda L: factorial(3) // int(np.prod([factorial(v) for v in Counter(L).values()]))
df["exp_w"] = df["labels"].apply(lambda L: mult(L) * np.prod([p[i] for i in L]))
df["expected"] = df["exp_w"] * df["occ"].sum() / df["exp_w"].sum()
df["enrichment"] = df["occ"] / df["expected"].replace(0, np.nan)
df["motif"] = df["labels"].apply(lambda L: "+".join(sorted(id2ct[i] for i in L)))

df["rank_raw"] = df["occ"].rank(ascending=False).astype(int)
df["rank_enrich"] = df["enrichment"].rank(ascending=False)
keep = ["motif", "occ", "expected", "enrichment", "rank_raw", "rank_enrich"]
df[keep].sort_values("enrichment", ascending=False).to_csv(OUT / "abundance_vs_enrichment_size3.csv", index=False)

print("=== TOP 5 by RAW occurrence (what TrimNN reports) ===")
print(df.sort_values("occ", ascending=False).head(5)[["motif", "occ"]].to_string(index=False))
print("\n=== TOP 5 by ENRICHMENT (abundance-corrected) ===")
print(df.sort_values("enrichment", ascending=False).head(5)[["motif", "occ", "expected", "enrichment"]].to_string(index=False))
print("\n-> ranking inverts. Caveat: ratios unstable at low counts; see 03_motif_statistics.py.")
