## Rule: plot_volcano
## Volcano plot of log2 fold change vs. adjusted p-value.

log <- file(snakemake@log[[1]], open = "wt")
sink(log)
sink(log, type = "message")

library(EnhancedVolcano)
library(ggplot2)

res <- read.csv(snakemake@input[["results"]], row.names = 1)

p <- EnhancedVolcano(res,
  lab = rownames(res),
  x = "log2FoldChange",
  y = "padj",
  pCutoff = snakemake@params[["p_cutoff"]],
  FCcutoff = snakemake@params[["fc_cutoff"]],
  title = snakemake@params[["title"]]
)
ggsave(snakemake@output[["plot"]], plot = p, width = 7, height = 6)
