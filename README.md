# covid-rnaseq-snakemake

[![Snakemake dry-run](https://github.com/sanzhsanzh-sua/covid-rnaseq/actions/workflows/dry-run.yml/badge.svg)](https://github.com/sanzhsanzh-sua/covid-rnaseq/actions/workflows/dry-run.yml)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

Snakemake workflow for differential expression analysis of GSE152418
(COVID-19 PBMC RNA-seq): download from GEO -> counts matrix -> QC ->
DESeq2 -> visualization -> GO enrichment.

## Pipeline

```
download_series_matrix ─┐
download_supp_files ────┴─> build_counts_matrix ─┐
                                                   ├─> prepare_coldata ─┐
                                                   └─> qc_filter ───────┴─> run_deseq2 ─┬─> plot_pca
                                                                                        ├─> plot_volcano
                                                                                        ├─> plot_heatmap
                                                                                        └─> go_enrichment
```

- **download_series_matrix** / **download_supp_files** — fetch the GEO
  series matrix (sample metadata) and supplementary raw-counts file(s).
- **build_counts_matrix** — assembles a genes x samples counts matrix.
  GSE152418 ships this as a single combined file; the script also handles
  the per-sample-file layout some GEO series use instead.
- **prepare_coldata** — builds the sample condition table (COVID-19 /
  Healthy / Convalescent), aligned to the counts matrix.
- **qc_filter** — reports sequencing depth per sample and drops
  low-expression genes.
- **run_deseq2** — differential expression testing.
- **plot_pca** / **plot_volcano** / **plot_heatmap** — visualization.
- **go_enrichment** — GO (Biological Process) over-representation analysis
  on the significant DE genes.

## Setup

Requires [conda](https://docs.conda.io/) or a conda-compatible tool
(mamba, miniforge).

```bash
conda env create -f environment.yaml
conda activate covid-rnaseq
```

## Run

```bash
snakemake --cores 4
```

Dry-run to preview the DAG without executing anything:

```bash
snakemake -n
```

Outputs are written to `data/` (downloads), `results/` (tables, RDS
objects), and `plots/` (PNG figures); per-rule logs go to `logs/`.

## Configuration

All parameters live in `config.yaml`:

- `gse_id` — GEO series accession.
- `condition_column` — phenotype column holding the sample group. Inspect
  `data/<gse_id>_pheno.csv` after the first run if you point this workflow
  at a different GEO series and need to adjust it.
- `contrast` — DESeq2 contrast `[factor, numerator level, denominator level]`;
  levels must exactly match the values in `condition_column`.
- `min_count_sum` — low-expression gene filter threshold.
- `heatmap_top_n` — number of top DE genes shown in the heatmap.
- `volcano` — `EnhancedVolcano` plot parameters.
- `go_enrichment` — `clusterProfiler::enrichGO` parameters, including
  `key_type` (gene ID type used by the counts matrix, e.g. `ENSEMBL` or
  `SYMBOL`).

## Note on the original script

`covid_rnaseq_analysis.R` is the original standalone script this workflow
was derived from; it is kept for reference. The values pre-filled in
`config.yaml` (contrast level `"Healthy"`, GO `key_type: "ENSEMBL"`, the
combined-matrix handling in `build_counts_matrix.R`) reflect GSE152418's
actual data layout, confirmed by an end-to-end run — not just the
assumptions in the original script's comments.
