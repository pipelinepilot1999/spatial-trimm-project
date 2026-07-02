# Bug note: `csv2gml.py -prune` silently deletes all edges on physical-unit coordinates

**Where:** upstream TrimNN `csv2gml.py` (graph construction from the input CSV).

**Symptom:** Running `csv2gml.py -target crop.csv -out crop.gml -prune False` produced a graph
with all 3,968 nodes but **0 edges**. Motif search on an edgeless graph is meaningless (no
neighborhoods exist), so this would have silently invalidated every downstream result.

**Root cause — two chained problems:**

1. **argparse `type=bool` bug.** `csv2gml.py` declares `-prune` as `type=bool`. argparse then
   evaluates `bool("False")`, and in Python **any non-empty string is truthy**, so
   `-prune False` sets `prune = True`. Pruning ran even though we asked for it off.

2. **The prune threshold collapses on small (physical-unit) coordinates.** The prune path
   classifies Delaunay edges as outliers using a `scipy.stats.lognorm.ppf(0.99, ...)` threshold
   fit to edge *lengths*. TrimNN's demo uses large integer pixel coordinates (edge lengths ~20–100),
   where the threshold behaves. Our ABC Atlas MERFISH coordinates are in **millimetres**
   (edge lengths ~0.01–0.05). On that scale the threshold computation misfires and marks
   **every** edge as an outlier → all triangles filtered → 0 edges.

**Verification that the tool is otherwise fine:** running `csv2gml.py` on the shipped
`demo_data.csv` (pixel coordinates) with pruning produced 2,177 edges (≈ the authors' shipped
2,211). The Delaunay step on our data was also correct in isolation (7,883 simplices, 11,850
undirected edges). Only the prune path on mm-scale coordinates was pathological.

**Fix used:** omit `-prune` entirely (default `False` → uses the full Delaunay triangulation).
Rebuild produced the expected 3,968 nodes / 11,850 edges.

**Implication / TODO:** if edge noise-reduction is wanted later (long artefact edges at the
crop boundary), do NOT rely on `-prune`; implement a unit-aware pruning step (e.g. drop edges
longer than a fixed micron distance, or the 99th percentile of the actual edge-length
distribution) and treat it as a sensitivity check, not a default.
