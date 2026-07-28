## Rule: run_deseq2
## Runs the DESeq2 differential expression analysis and writes out the
## results table plus the DESeqDataSet / variance-stabilized objects used
## by the downstream plotting rules.

log <- file(snakemake@log[[1]], open = "wt")
sink(log)
sink(log, type = "message")

library(DESeq2)

counts_mat <- read.csv(snakemake@input[["counts"]], row.names = 1, check.names = FALSE)
coldata <- read.csv(snakemake@input[["coldata"]], row.names = 1, check.names = FALSE)
coldata$condition <- factor(coldata$condition)

contrast <- unlist(snakemake@params[["contrast"]])  # c("condition", "COVID-19", "healthy")

dds <- DESeqDataSetFromMatrix(
  countData = round(as.matrix(counts_mat)),
  colData   = coldata,
  design    = ~ condition
)

dds <- DESeq(dds)
res <- results(dds, contrast = contrast)
res <- res[order(res$padj), ]

write.csv(as.data.frame(res), snakemake@output[["results"]])
saveRDS(dds, snakemake@output[["dds"]])

vsd <- vst(dds, blind = FALSE)
saveRDS(vsd, snakemake@output[["vsd"]])
