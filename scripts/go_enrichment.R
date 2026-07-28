## Rule: go_enrichment
## GO (Biological Process) over-representation analysis on the significant
## differentially expressed genes.

log <- file(snakemake@log[[1]], open = "wt")
sink(log)
sink(log, type = "message")

library(clusterProfiler)
library(org.Hs.eg.db)

res <- read.csv(snakemake@input[["results"]], row.names = 1)

padj_cutoff <- snakemake@params[["padj_cutoff"]]
lfc_cutoff  <- snakemake@params[["lfc_cutoff"]]

sig_genes <- rownames(res)[which(res$padj < padj_cutoff & abs(res$log2FoldChange) > lfc_cutoff)]

ego <- enrichGO(
  gene          = sig_genes,
  OrgDb         = org.Hs.eg.db,
  keyType       = snakemake@params[["key_type"]],
  ont           = snakemake@params[["ontology"]],
  pAdjustMethod = snakemake@params[["p_adjust_method"]],
  pvalueCutoff  = padj_cutoff
)

write.csv(as.data.frame(ego), snakemake@output[["go_results"]])
