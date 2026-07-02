# 08_cellchat.R
# Cell-cell communication in the section-.35 cortical crop via CellChat (CellChatDB.mouse,
# Secreted Signaling). Group cells by their 8 classes; infer ligand-receptor signaling.
#
# CAVEATS to keep in the write-up:
#  - CellChat infers communication from CO-EXPRESSION of ligand/receptor genes; it is a
#    PREDICTION, not proof of signaling.
#  - We only have a 550-gene MERFISH panel, so only LR pairs whose BOTH partners are on the
#    panel can ever be detected -> communication space is severely truncated (panel-bound).
#
# Runnable once CellChat is installed (from GitHub jinworks/CellChat).

suppressMessages({library(CellChat); library(Matrix)})
PROC<-"data/processed"; OUT<-"results/cellchat"; dir.create(OUT, showWarnings=FALSE, recursive=TRUE)

counts <- as.matrix(read.csv(file.path(PROC,"de_counts_sec35.csv"), row.names=1, check.names=FALSE))
genes  <- read.csv(file.path(PROC,"de_genes_sec35.csv"), stringsAsFactors=FALSE)
meta   <- read.csv(file.path(PROC,"de_coldata_sec35.csv"), stringsAsFactors=FALSE)
rownames(meta) <- meta$cell_label; meta <- meta[colnames(counts), ]

# CellChatDB matches on gene SYMBOLS -> relabel rows ensembl->symbol, drop unmapped/dups
sym <- genes$symbol[match(rownames(counts), genes$ensembl_id)]
ok  <- !is.na(sym) & !duplicated(sym)
counts <- counts[ok, ]; rownames(counts) <- sym[ok]

# CellChat wants normalized (not raw) data
data.input <- normalizeData(counts)                     # library-size normalization
meta$cell_type <- factor(meta$cell_type)

cc <- createCellChat(object = data.input, meta = meta, group.by = "cell_type")
cc@DB <- subsetDB(CellChatDB.mouse, search = "Secreted Signaling")
cc <- subsetData(cc)
cc <- identifyOverExpressedGenes(cc)
cc <- identifyOverExpressedInteractions(cc)
cc <- computeCommunProb(cc)
cc <- filterCommunication(cc, min.cells = 10)
cc <- computeCommunProbPathway(cc)
cc <- aggregateNet(cc)

saveRDS(cc, file.path(OUT,"cellchat_sec35.rds"))
lr <- subsetCommunication(cc)                            # significant LR pairs
write.csv(lr, file.path(OUT,"significant_LR_pairs.csv"), row.names=FALSE)
cat(sprintf("signaling pathways detected: %d\n", length(cc@netP$pathways)))
print(cc@netP$pathways)
cat(sprintf("significant LR interactions: %d\n", nrow(lr)))
cat("done -> results/cellchat/\n")
