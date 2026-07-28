## Rule: plot_heatmap
## Heatmap of the top N differentially expressed genes (by adjusted p-value).

log <- file(snakemake@log[[1]], open = "wt")
sink(log)
sink(log, type = "message")

library(DESeq2)
library(pheatmap)

vsd <- readRDS(snakemake@input[["vsd"]])
res <- read.csv(snakemake@input[["results"]], row.names = 1)
coldata <- read.csv(snakemake@input[["coldata"]], row.names = 1, check.names = FALSE)

top_n <- snakemake@params[["top_n"]]
top_genes <- rownames(res)[seq_len(min(top_n, nrow(res)))]

pheatmap(assay(vsd)[top_genes, ],
         annotation_col = coldata,
         show_rownames = TRUE,
         filename = snakemake@output[["plot"]])
