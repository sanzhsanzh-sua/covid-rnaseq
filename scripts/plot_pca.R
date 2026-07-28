## Rule: plot_pca
## PCA plot of variance-stabilized counts, to check for batch effects/outliers.

log <- file(snakemake@log[[1]], open = "wt")
sink(log)
sink(log, type = "message")

library(DESeq2)
library(ggplot2)

vsd <- readRDS(snakemake@input[["vsd"]])

p <- plotPCA(vsd, intgroup = "condition")
ggsave(snakemake@output[["plot"]], plot = p, width = 6, height = 5)
