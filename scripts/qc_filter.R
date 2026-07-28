## Rule: qc_filter
## Reports sequencing depth per sample and removes low-expression genes
## before differential expression testing.

log <- file(snakemake@log[[1]], open = "wt")
sink(log)
sink(log, type = "message")

counts_mat <- read.csv(snakemake@input[["counts"]], row.names = 1, check.names = FALSE)
min_count_sum <- snakemake@params[["min_count_sum"]]

qc_con <- file(snakemake@output[["qc_summary"]], open = "wt")
writeLines("Sequencing depth per sample (colSums):", qc_con)
writeLines(capture.output(print(summary(colSums(counts_mat)))), qc_con)
close(qc_con)

counts_filtered <- counts_mat[rowSums(counts_mat) > min_count_sum, ]

write.csv(counts_filtered, snakemake@output[["counts_filtered"]], row.names = TRUE)
