## Rule: build_counts_matrix
## Builds a single genes x samples counts matrix from the GEO supplementary
## file(s). GEO series store raw counts either as one combined matrix file
## or as one file per sample (<GSM_ID>_<name>.counts.txt.gz) -- for example
## GSE152418 uses a single combined "*RawCounts.txt.gz" matrix, so both
## layouts are handled here.

log <- file(snakemake@log[[1]], open = "wt")
sink(log)
sink(log, type = "message")

supp_dir <- snakemake@input[["supp_dir"]]

count_files <- list.files(
  supp_dir,
  pattern = "[Cc]ounts.*\\.txt\\.gz$",
  full.names = TRUE,
  recursive = TRUE
)

if (length(count_files) == 0) {
  stop(
    "No counts file(s) found in ", supp_dir,
    " -- check the supplementary file naming for this GEO series."
  )
}

if (length(count_files) == 1) {
  ## Single file: already a full genes x samples matrix.
  counts_mat <- read.table(count_files[1], header = TRUE, row.names = 1,
                            check.names = FALSE, sep = "\t")
} else {
  ## Multiple files: one per sample, each with a single count column.
  read_one <- function(f) {
    read.table(f, header = FALSE, row.names = 1,
               col.names = c("gene", basename(f)))
  }
  count_list <- lapply(count_files, read_one)
  counts_mat <- do.call(cbind, count_list)
  colnames(counts_mat) <- gsub("_.*", "", basename(count_files))  # GSM ID as sample name
}

write.csv(counts_mat, snakemake@output[["counts"]], row.names = TRUE)
