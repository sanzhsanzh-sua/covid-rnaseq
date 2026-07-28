## Rule: download_series_matrix
## Downloads the GEO series matrix and writes out the sample phenotype table.

log <- file(snakemake@log[[1]], open = "wt")
sink(log)
sink(log, type = "message")

library(GEOquery)

gse_id <- snakemake@params[["gse_id"]]
data_dir <- snakemake@params[["data_dir"]]

dir.create(data_dir, showWarnings = FALSE, recursive = TRUE)

gse <- getGEO(gse_id, GSEMatrix = TRUE, destdir = data_dir)
pheno <- pData(gse[[1]])

write.csv(pheno, snakemake@output[["pheno"]], row.names = TRUE)
