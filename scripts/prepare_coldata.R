## Rule: prepare_coldata
## Builds the sample metadata table (condition: e.g. COVID-19 vs Healthy)
## from the GEO phenotype data, aligned to the sample order of the counts
## matrix.

log <- file(snakemake@log[[1]], open = "wt")
sink(log)
sink(log, type = "message")

pheno <- read.csv(snakemake@input[["pheno"]], row.names = 1, check.names = FALSE)
counts_mat <- read.csv(snakemake@input[["counts"]], row.names = 1, check.names = FALSE)

condition_column <- snakemake@params[["condition_column"]]

if (!condition_column %in% colnames(pheno)) {
  stop(
    "Column '", condition_column, "' not found in the phenotype table. ",
    "Inspect ", snakemake@input[["pheno"]],
    " and update 'condition_column' in config.yaml."
  )
}

## The counts matrix may use either the GEO accession (GSM ID, when GEO
## stores one file per sample) or the sample title (when GEO stores a single
## combined matrix, as for GSE152418) as its column names -- match on
## whichever the counts matrix actually uses.
if (all(colnames(counts_mat) %in% rownames(pheno))) {
  sample_ids <- rownames(pheno)
} else if ("title" %in% colnames(pheno) && all(colnames(counts_mat) %in% pheno$title)) {
  sample_ids <- pheno$title
} else {
  stop(
    "Could not match counts matrix sample names to the phenotype table ",
    "(checked GEO accession and 'title')."
  )
}

coldata <- data.frame(
  row.names = sample_ids,
  condition = factor(pheno[[condition_column]])
)

## Align sample order between counts matrix and metadata
coldata <- coldata[colnames(counts_mat), , drop = FALSE]

write.csv(coldata, snakemake@output[["coldata"]], row.names = TRUE)
