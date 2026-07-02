# 07_enrichment.R
# GO (BP) + Reactome enrichment on the OPC-Oligo motif DEGs.
#
# KEY METHOD CHOICE (differentiator): the enrichment UNIVERSE is the 550-gene MERFISH panel,
# NOT the whole genome. A targeted panel is a biased gene set; testing against a genome-wide
# background would call pathways "enriched" merely because the panel was curated toward them.
# Using the panel as background asks the honest question: given only these 550 genes could be
# detected, are the DEGs enriched for a pathway beyond what the panel itself is?
#
# Run on both DEG sets (B confound-controlled n~48, C pseudobulk n~17) for comparison.

suppressMessages({library(clusterProfiler); library(ReactomePA); library(org.Mm.eg.db)})
PROC<-"data/processed"; DEDIR<-"results/de"; OUT<-"results/enrichment"
dir.create(OUT, showWarnings=FALSE, recursive=TRUE)

genes <- read.csv(file.path(PROC,"de_genes_sec35.csv"), stringsAsFactors=FALSE)
# universe = panel genes as ENTREZ
uni <- bitr(genes$ensembl_id, "ENSEMBL","ENTREZID", org.Mm.eg.db)$ENTREZID |> unique()
cat(sprintf("panel universe: %d genes -> %d mapped to Entrez\n", nrow(genes), length(uni)))

enrich_set <- function(csv, tag){
  d <- read.csv(csv, row.names=1)
  sig <- rownames(d)[!is.na(d$padj) & d$padj < 0.05]
  eg <- bitr(sig, "ENSEMBL","ENTREZID", org.Mm.eg.db)$ENTREZID |> unique()
  cat(sprintf("\n[%s] %d DEGs -> %d Entrez\n", tag, length(sig), length(eg)))
  if(length(eg) < 3){ cat("  too few genes for enrichment\n"); return(invisible()) }
  go <- enrichGO(eg, OrgDb=org.Mm.eg.db, ont="BP", universe=uni,
                 pAdjustMethod="BH", pvalueCutoff=0.1, qvalueCutoff=0.2, readable=TRUE)
  re <- enrichPathway(eg, organism="mouse", universe=uni,
                      pAdjustMethod="BH", pvalueCutoff=0.1, qvalueCutoff=0.2, readable=TRUE)
  if(!is.null(go) && nrow(go)>0){ write.csv(as.data.frame(go), file.path(OUT,paste0(tag,"_GO_BP.csv")))
    cat("  GO BP top terms:\n"); print(head(as.data.frame(go)[,c("Description","p.adjust","Count")],5)) }
  else cat("  GO BP: no terms at q<0.2\n")
  if(!is.null(re) && nrow(re)>0){ write.csv(as.data.frame(re), file.path(OUT,paste0(tag,"_Reactome.csv")))
    cat("  Reactome top terms:\n"); print(head(as.data.frame(re)[,c("Description","p.adjust","Count")],5)) }
  else cat("  Reactome: no terms at q<0.2\n")
}

enrich_set(file.path(DEDIR,"B_focal_cell_in_vs_out.csv"), "B_confound_ctrl")
enrich_set(file.path(DEDIR,"C_pseudobulk_in_vs_out.csv"), "C_pseudobulk")
cat("\ndone -> results/enrichment/\n")
