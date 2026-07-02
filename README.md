# Spatial Cell-Motif Discovery & Downstream Characterization

Rebuild of a [TrimNN](https://github.com/yuyang-0825/TrimNN)-based pipeline for discovering
recurrent **spatial cell-type motifs** (Cellular Community motifs) in spatial transcriptomics,
then interrogating whether those motifs are **molecularly, functionally, communicatively, and
statistically real** — built from scratch on a new dataset as a defensible, honest portfolio.

## The question this asks
A cell's behavior depends not only on its own expression but on *which cell types surround it*.
scRNA-seq dissociates tissue and loses that spatial context; spatial transcriptomics restores it.
TrimNN finds neighborhood patterns of cell types (motifs) that recur more than chance would predict.
This pipeline then asks, of each motif:
- **Molecular** — do motif cells differ in expression from non-motif cells? (DESeq2, **with a pseudobulk
  comparison to avoid pseudoreplication** — the key methodological differentiator)
- **Communicative** — are the neighboring cell types plausibly signaling? (CellChat, CellChatDB.mouse)
- **Functional** — what pathways are enriched? (GO / clusterProfiler + Reactome)
- **Statistical** — is the motif's co-occurrence non-random? (Fisher's exact + Cramér's V)

## Dataset
Mouse-brain MERFISH, single section, subsampled to bound CPU runtime. Chosen to match TrimNN's demo
tissue, ship with cell-type annotations, and work with `CellChatDB.mouse` out of the box.

## Documented limitations (the honesty is the point)
- **Pseudoreplication** — individual cells are spatially correlated; addressed via pseudobulk DE.
- **Cell-type confounding** — motif-vs-nonmotif DEGs may just re-detect cell-type identity; check/stratify.
- **CellChat = co-expression prediction**, not proof of signaling; database-bound.
- **Enrichment shared-gene artifacts** — enriched ≠ mechanistic.
- **Single section, no biological replicate** — conclusions are exploratory.

## Status
- [x] Environment + TrimNN CPU-only install (conda `TrimNNEnv`, python 3.9)
- [x] Proof-of-life on TrimNN `demo_data` (Functions 1 & 2 run end-to-end on CPU)
- [x] Acquire + spatially crop a mouse-brain MERFISH section (ABC Atlas `MERFISH-C57BL6J-638850`, section .35)
- [x] Motif discovery on real data — size-3 and size-4 (greedy)
- [x] Statistical validation — permutation null + Fisher's exact + Cramér's V
- [x] Differential expression — naive vs confound-controlled vs **pseudobulk**
- [x] GO / Reactome enrichment (panel-gene background)
- [x] Cell-cell communication — squidpy `ligrec` (honest null: panel lacks complete LR pairs)
- [~] CellChat proper — script + input ready, run deferred (external-install guardrail; needs approval)

## Results (section .35 cortical crop, 3,968 cells, 8 classes)

> 📊 **Figure walkthrough with plots: [`RESULTS.md`](RESULTS.md)**

**Headline (motif discovery):** TrimNN's top "overrepresented" motif is an *abundance artifact*
(`Astro+Glut+Glut`). Correcting for abundance with a permutation null recovers real cortical
biology — **layer-6b neuron clustering (z=133)**, glial clustering, and neuron↔glia spatial
exclusion (`IT-ET~Astro` odds 0.09). Global cell-type spatial organization: Cramér's V = 0.28.

**Methodological headline (DE):** for the OPC-Oligo motif, DEG counts collapse as each bias is
removed —

| Analysis | DEGs (padj<0.05) | fix applied |
|---|---|---|
| Naive (motif vs rest) | **222** | — (confounded + pseudoreplicated) |
| Confound-controlled | **48** | cell-type confound removed (−78%) |
| **Pseudobulk** | **17** | pseudoreplication removed (−65% more) |

That 222 → 17 collapse is the project's argument: naive spatial DE is severely inflated.

**Cross-validation:** TrimNN's greedy size-4 search independently grew from the `Immune³`
(microglia) motif — the same motif our from-scratch permutation null flags as most over-represented
by ratio. Two independent methods agreeing on the biologically interesting pattern.

Full accounting of assumptions and residual caveats in [`docs/LIMITATIONS.md`](docs/LIMITATIONS.md).

## Layout
```
TrimNN/     upstream tool, pinned as a git submodule (not my code)
scripts/    my pipeline code
data/       datasets (raw data gitignored)
results/    my outputs, including demo proof-of-life
docs/       notes, method rationale, limitation write-ups
```

## Reproducing the environment
```
conda create -n TrimNNEnv python=3.9 -y && conda activate TrimNNEnv
conda install pytorch==1.13.1 torchvision==0.14.1 torchaudio==0.13.1 cpuonly -c pytorch -y
pip install dgl==1.1.2 -f https://data.dgl.ai/wheels/repo.html
cd TrimNN && pip install -r requirements.txt
```
