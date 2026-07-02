# Limitations & how each is addressed

Honesty about failure modes is the point of this project. Each limitation below is stated, and
where possible **addressed in code**, with the residual caveat that remains.

## 1. Pseudoreplication (ADDRESSED — key differentiator)
Individual cells in a tissue are spatially autocorrelated; treating each as an independent
replicate in DESeq2 inflates significance.
**Addressed:** pseudobulk aggregation — cells summed into spatial-tile pseudo-samples, DESeq2 on
those (`06_de_deseq2.R`, analysis C). Effect was large: **222 → 17 DEGs** across the fixes.
**Residual caveat:** the pseudo-replicates are spatial tiles from *one* section, NOT independent
biological replicates. Pseudobulk fixes the *statistical* independence problem but cannot
manufacture biological replication — see #7.

## 2. Cell-type confounding (ADDRESSED)
"Motif cells vs non-motif cells" mostly contrasts different cell TYPES, so DEGs just re-discover
cell-type markers, not microenvironment effects.
**Addressed:** focal-type contrast — compare one cell type's cells that are IN vs OUT of the motif
(`06_de_deseq2.R`, analyses B/C). Removed ~78% of naive "DEGs" (222 → 48).
**Residual caveat:** conditions on the focal cell type chosen (OPC-Oligo here); other types may
behave differently.

## 3. TrimNN "over-representation" is abundance-driven, not enrichment (ADDRESSED)
TrimNN Function 2 ranks/selects motifs by raw predicted occurrence; on real data its top
"overrepresented" motif was simply the most abundant cell types (`Astro+Glut+Glut`).
**Addressed:** independent permutation-null enrichment (`03_motif_statistics.py`) — shuffle labels
on the fixed graph, so enrichment controls for abundance. This inverts the ranking to real spatial
biology (L6b layering z=133, glial clustering).
**Residual caveat:** TrimNN's neural prediction correlates only r=0.52 with exact triangle counts —
its approximation is imperfect on this data; we rely on exact counts for size-3.

## 4. Statistical-significance inflation at large N (ADDRESSED by reporting)
With ~7,900 triangles the permutation test is so powerful that 107/108 motifs are "significant"
at FDR<0.05.
**Addressed:** we lead with EFFECT SIZE (z-score, Fisher odds, Cramér's V), not p-values.
**Residual caveat:** effect-size thresholds are judgement calls, not bright lines.

## 5. CellChat = prediction, not proof; panel-bound (DEFERRED)
CellChat infers signaling from ligand/receptor CO-EXPRESSION — a hypothesis, not evidence of
actual signaling. With a ~550-gene panel, only LR pairs whose both partners are on the panel can
be found, truncating the communication space.
**Status:** script ready (`08_cellchat.R`) but NOT run — CellChat installs only from GitHub and the
environment guardrail blocked the unauthorized external source. Runs on explicit approval.

## 6. Enrichment: shared-gene artifacts + panel background (ADDRESSED, underpowered)
Genome-wide GO/Reactome background would call pathways enriched merely because the MERFISH panel
was curated toward them.
**Addressed:** enrichment universe = the 550 panel genes (`07_enrichment.R`).
**Residual caveat:** with an honest 17-gene pseudobulk DEG set, enrichment is underpowered — no
significant terms. Enriched ≠ mechanistic regardless.

## 7. Single section, no biological replicate (INHERENT)
One coronal section from one animal. Spatial tiles are pseudo-replicates, not biological ones.
**Consequence:** all conclusions are exploratory / hypothesis-generating, not confirmatory.

## 8. Cortical crop, not whole section (SCOPE)
Analysis is one contiguous ~4k-cell cortical patch (chosen to bound CPU runtime and match the demo
tissue). Motifs/DE reflect that patch; boundary edges of the Delaunay graph are included un-pruned.

## 9. Upstream tool bug found & fixed
`csv2gml -prune` silently deleted all edges on millimetre-scale coordinates — see
`docs/csv2gml_prune_bug.md`. Fixed by omitting `-prune`.
