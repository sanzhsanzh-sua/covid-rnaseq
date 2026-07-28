# ============================================================
# COVID-19 PBMC RNA-seq differential expression workflow (GSE152418)
# Pipeline: download -> counts matrix -> metadata -> QC -> DESeq2 ->
#           visualization -> GO enrichment
#
# All R packages required by the scripts/ below are provided by the conda
# environment defined in environment.yaml (run with `snakemake --use-conda`).
# ============================================================

configfile: "config.yaml"

DATA_DIR = config["data_dir"]
RESULTS_DIR = config["results_dir"]
PLOTS_DIR = config["plots_dir"]
GSE_ID = config["gse_id"]


rule all:
    input:
        f"{RESULTS_DIR}/results_DESeq2.csv",
        f"{RESULTS_DIR}/GO_enrichment.csv",
        f"{RESULTS_DIR}/qc_summary.txt",
        f"{PLOTS_DIR}/pca_plot.png",
        f"{PLOTS_DIR}/volcano_plot.png",
        f"{PLOTS_DIR}/heatmap_top30.png",


# --- 1a. Download the GEO series matrix (sample phenotype/metadata) ---
rule download_series_matrix:
    output:
        pheno=f"{DATA_DIR}/{GSE_ID}_pheno.csv",
    params:
        gse_id=GSE_ID,
        data_dir=DATA_DIR,
    log:
        "logs/download_series_matrix.log",
    script:
        "scripts/download_series_matrix.R"


# --- 1b. Download (and unpack) the GEO supplementary raw-counts files ---
rule download_supp_files:
    output:
        supp_dir=directory(f"{DATA_DIR}/{GSE_ID}"),
    params:
        gse_id=GSE_ID,
        data_dir=DATA_DIR,
    log:
        "logs/download_supp_files.log",
    script:
        "scripts/download_supp_files.R"


# --- 2. Assemble the per-sample count files into a single counts matrix ---
rule build_counts_matrix:
    input:
        supp_dir=f"{DATA_DIR}/{GSE_ID}",
    output:
        counts=f"{RESULTS_DIR}/counts_raw.csv",
    log:
        "logs/build_counts_matrix.log",
    script:
        "scripts/build_counts_matrix.R"


# --- 3. Build sample metadata (condition: COVID-19 vs healthy) ---
rule prepare_coldata:
    input:
        pheno=f"{DATA_DIR}/{GSE_ID}_pheno.csv",
        counts=f"{RESULTS_DIR}/counts_raw.csv",
    output:
        coldata=f"{RESULTS_DIR}/coldata.csv",
    params:
        condition_column=config["condition_column"],
        gse_id=GSE_ID,
    log:
        "logs/prepare_coldata.log",
    script:
        "scripts/prepare_coldata.R"


# --- 4. QC: sequencing depth summary + low-expression gene filter ---
rule qc_filter:
    input:
        counts=f"{RESULTS_DIR}/counts_raw.csv",
    output:
        counts_filtered=f"{RESULTS_DIR}/counts_filtered.csv",
        qc_summary=f"{RESULTS_DIR}/qc_summary.txt",
    params:
        min_count_sum=config["min_count_sum"],
    log:
        "logs/qc_filter.log",
    script:
        "scripts/qc_filter.R"


# --- 5. DESeq2 differential expression ---
rule run_deseq2:
    input:
        counts=f"{RESULTS_DIR}/counts_filtered.csv",
        coldata=f"{RESULTS_DIR}/coldata.csv",
    output:
        results=f"{RESULTS_DIR}/results_DESeq2.csv",
        dds=f"{RESULTS_DIR}/dds.rds",
        vsd=f"{RESULTS_DIR}/vsd.rds",
    params:
        contrast=config["contrast"],
    threads: config.get("deseq2_threads", 1)
    log:
        "logs/run_deseq2.log",
    script:
        "scripts/run_deseq2.R"


# --- 6a. Visualization: PCA plot ---
rule plot_pca:
    input:
        vsd=f"{RESULTS_DIR}/vsd.rds",
    output:
        plot=f"{PLOTS_DIR}/pca_plot.png",
    log:
        "logs/plot_pca.log",
    script:
        "scripts/plot_pca.R"


# --- 6b. Visualization: volcano plot ---
rule plot_volcano:
    input:
        results=f"{RESULTS_DIR}/results_DESeq2.csv",
    output:
        plot=f"{PLOTS_DIR}/volcano_plot.png",
    params:
        p_cutoff=config["volcano"]["p_cutoff"],
        fc_cutoff=config["volcano"]["fc_cutoff"],
        title=config["volcano"]["title"],
    log:
        "logs/plot_volcano.log",
    script:
        "scripts/plot_volcano.R"


# --- 6c. Visualization: heatmap of top differentially expressed genes ---
rule plot_heatmap:
    input:
        vsd=f"{RESULTS_DIR}/vsd.rds",
        results=f"{RESULTS_DIR}/results_DESeq2.csv",
        coldata=f"{RESULTS_DIR}/coldata.csv",
    output:
        plot=f"{PLOTS_DIR}/heatmap_top30.png",
    params:
        top_n=config["heatmap_top_n"],
    log:
        "logs/plot_heatmap.log",
    script:
        "scripts/plot_heatmap.R"


# --- 7. Functional analysis: GO enrichment ---
rule go_enrichment:
    input:
        results=f"{RESULTS_DIR}/results_DESeq2.csv",
    output:
        go_results=f"{RESULTS_DIR}/GO_enrichment.csv",
    params:
        padj_cutoff=config["go_enrichment"]["padj_cutoff"],
        lfc_cutoff=config["go_enrichment"]["lfc_cutoff"],
        ontology=config["go_enrichment"]["ontology"],
        key_type=config["go_enrichment"]["key_type"],
        p_adjust_method=config["go_enrichment"]["p_adjust_method"],
    log:
        "logs/go_enrichment.log",
    script:
        "scripts/go_enrichment.R"
