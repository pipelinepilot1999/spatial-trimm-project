# 06_de_deseq2.R
# Differential expression for the OPC-Oligo focal motif, three ways, to make the two
# methodological problems (and their fixes) VISIBLE:
#
#   A. NAIVE          motif cells vs all other cells, per-cell DESeq2
#                     -> confounded by cell type AND pseudoreplicated (the flawed original)
#   B. CONFOUND-CTRL  OPC-Oligo in-motif vs OPC-Oligo out-of-motif, per-cell DESeq2
#                     -> fixes confounding, STILL pseudoreplicated
#   C. PSEUDOBULK     same focal contrast, cells aggregated per spatial tile x group
#                     -> fixes pseudoreplication (tiles = genuine replicates), design ~tile+group
#
# The story is in the DEG counts: A (huge, inflated) -> B -> C (honest, few).
# NOTE: MERFISH panel counts are sparse -> size factors via type="poscounts".

suppressMessages(library(DESeq2))
PROC <- "data/processed"; OUT <- "results/de"; dir.create(OUT, showWarnings=FALSE, recursive=TRUE)

counts <- as.matrix(read.csv(file.path(PROC,"de_counts_sec35.csv"), row.names=1, check.names=FALSE))
g <- read.csv(file.path(PROC,"de_groups_sec35.csv"), stringsAsFactors=FALSE)
g$is_focal        <- g$is_focal == "True"
g$in_focal_motif  <- g$in_focal_motif == "True"
rownames(g) <- g$cell_label
g <- g[colnames(counts), ]                                  # align to count columns
stopifnot(identical(rownames(g), colnames(counts)))

deg_count <- function(res, a=0.05) sum(!is.na(res$padj) & res$padj < a)

run_cell <- function(cnt, cond, ref, label){                # per-cell DESeq2
  cond <- factor(cond); cond <- relevel(cond, ref=ref)
  cd <- data.frame(cond=cond, row.names=colnames(cnt))
  dds <- DESeqDataSetFromMatrix(cnt, cd, ~cond)
  dds <- dds[rowSums(counts(dds)) >= 10, ]
  sizeFactors(dds) <- estimateSizeFactorsForMatrix(counts(dds), type="poscounts")
  dds <- DESeq(dds, quiet=TRUE, fitType="local")
  res <- results(dds, name=resultsNames(dds)[length(resultsNames(dds))])
  write.csv(as.data.frame(res[order(res$padj),]), file.path(OUT, paste0(label,".csv")))
  cat(sprintf("  %-28s samples=%4d genes=%3d  DEGs(padj<0.05)=%d\n",
              label, ncol(cnt), nrow(dds), deg_count(res)))
  deg_count(res)
}

# ---- A: naive (all cells) ----
cat("A. NAIVE (confounded + pseudoreplicated):\n")
A <- run_cell(counts, g$naive_group, ref="other", "A_naive_motif_vs_other")

# ---- B: confound-controlled, per cell ----
cat("B. CONFOUND-CONTROLLED, per-cell (still pseudoreplicated):\n")
foc <- g$is_focal
B <- run_cell(counts[, foc], ifelse(g$in_focal_motif[foc],"in","out"), ref="out",
              "B_focal_cell_in_vs_out")

# ---- C: pseudobulk (aggregate focal cells per tile x group) ----
cat("C. PSEUDOBULK (tiles = replicates, ~tile+group):\n")
fg <- g[foc, ]
key <- paste(fg$tile, ifelse(fg$in_focal_motif,"in","out"), sep="__")
pb <- t(rowsum(t(counts[, foc]), group=key))                # genes x pseudo-samples
samp <- data.frame(row.names=colnames(pb),
                   tile  = sub("__.*","",colnames(pb)),
                   group = sub(".*__","",colnames(pb)))
# keep only tiles present in BOTH groups with >=8 cells each (paired replicates)
ncell <- table(key)
samp$n <- as.integer(ncell[rownames(samp)])
tiles_both <- intersect(samp$tile[samp$group=="in" & samp$n>=8],
                        samp$tile[samp$group=="out"& samp$n>=8])
keep <- samp$tile %in% tiles_both
pb <- pb[, keep]; samp <- samp[keep, ]
samp$tile <- factor(samp$tile); samp$group <- relevel(factor(samp$group), ref="out")
cat(sprintf("  pseudo-samples=%d across %d tiles\n", ncol(pb), length(tiles_both)))
dds <- DESeqDataSetFromMatrix(pb, samp, ~tile + group)
dds <- dds[rowSums(counts(dds)) >= 10, ]
dds <- DESeq(dds, quiet=TRUE)
resC <- results(dds, name="group_in_vs_out")
write.csv(as.data.frame(resC[order(resC$padj),]), file.path(OUT,"C_pseudobulk_in_vs_out.csv"))
C <- deg_count(resC)
cat(sprintf("  %-28s samples=%4d genes=%3d  DEGs(padj<0.05)=%d\n",
            "C_pseudobulk_in_vs_out", ncol(pb), nrow(dds), C))

cat(sprintf("\n=== DEG counts (padj<0.05):  A naive=%d  ->  B confound-ctrl=%d  ->  C pseudobulk=%d ===\n",
            A, B, C))
write.csv(data.frame(analysis=c("A_naive","B_confound_ctrl_cell","C_pseudobulk"),
                     degs=c(A,B,C)), file.path(OUT,"deg_summary.csv"), row.names=FALSE)
