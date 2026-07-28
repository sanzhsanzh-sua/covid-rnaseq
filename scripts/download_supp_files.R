## Rule: download_supp_files
## Downloads the GEO supplementary files (raw per-sample counts) and unpacks
## the .tar archive if GEO bundled the samples that way.
##
## NOTE: GSE152418 stores raw counts as a supplementary file rather than
## FASTQ, so this workflow uses it directly, without a read alignment step.

log <- file(snakemake@log[[1]], open = "wt")
sink(log)
sink(log, type = "message")

library(GEOquery)

gse_id <- snakemake@params[["gse_id"]]
data_dir <- snakemake@params[["data_dir"]]
supp_dir <- snakemake@output[["supp_dir"]]

dir.create(data_dir, showWarnings = FALSE, recursive = TRUE)

## getGEOSuppFiles() downloads into <baseDir>/<gse_id>/, i.e. supp_dir
getGEOSuppFiles(gse_id, baseDir = data_dir)

## GEO series with multiple samples often bundle them as a single .tar
tar_file <- list.files(supp_dir, pattern = "\\.tar$", full.names = TRUE)
if (length(tar_file) > 0) {
  untar(tar_file[1], exdir = supp_dir)
}
