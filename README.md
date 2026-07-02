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
envs/       conda environment files (3 envs)
scripts/    my pipeline code (numbered 00–09, run in order)
data/       raw/ (gitignored, fetched by 00) · processed/ (committed inputs)
results/    my outputs (motifs, stats, DE, enrichment, communication, figures)
docs/       method rationale + limitation write-ups
```

## Reproduce it

**1. Clone with the submodule**
```
git clone --recurse-submodules <repo-url> && cd spatial-trimm-project
```

**2. Build the three environments**
```
conda env create -f envs/TrimNNEnv.yml    # + DGL/requirements, see file header
conda env create -f envs/spatial.yml
conda env create -f envs/renv.yml
```

**3. Fetch the raw data** (gitignored; ~1.2 GB from the public ABC Atlas S3 bucket)
```
bash scripts/00_download_data.sh
```

**4. Run the pipeline** (each stage saves an output; `export DGLBACKEND=pytorch` for TrimNN)

| Stage | Command | Env |
|---|---|---|
| Cortical crop → TrimNN input | `python scripts/01_make_cortical_crop.py` | spatial |
| Build gml graph | `python TrimNN/csv2gml.py -target … -out …` (no `-prune`!) | TrimNNEnv |
| Motif discovery (size 3 / 4) | `python TrimNN/TrimNN.py -function specific_size … / all_size …` | TrimNNEnv |
| Abundance-vs-enrichment (exploratory) | `python scripts/02_motif_enrichment_explore.py` | spatial |
| Statistical validation (permutation, Fisher, Cramér's V) | `python scripts/03_motif_statistics.py` | spatial |
| DE input | `python scripts/04_build_de_input.py` | spatial |
| DE group prep | `python scripts/05_prep_de_groups.py` | spatial |
| DE + pseudobulk | `Rscript scripts/06_de_deseq2.R` | renv |
| GO / Reactome | `Rscript scripts/07_enrichment.R` | renv |
| Communication (squidpy) | `python scripts/08_communication_squidpy.py` | spatial |
| CellChat (deferred; needs GitHub install) | `Rscript scripts/08_cellchat.R` | renv + CellChat |
| Figures | `python scripts/09_figures.py` | spatial |

## License
MIT ([`LICENSE`](LICENSE)) for the pipeline code. TrimNN is third-party under its own MIT license.
