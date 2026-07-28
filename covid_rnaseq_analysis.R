## ============================================================
## Дифференциальная экспрессия генов: GSE152418 (COVID-19 PBMC RNA-seq)
## Пайплайн: скачивание -> QC -> DESeq2 -> визуализация -> GO enrichment
## ============================================================

## --- 0. Установка пакетов (запустить один раз) ---
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(c(
  "GEOquery",
  "DESeq2",
  "clusterProfiler",
  "org.Hs.eg.db",
  "EnhancedVolcano",
  "pheatmap"
), update = FALSE, ask = FALSE)

library(GEOquery)
library(DESeq2)
library(clusterProfiler)
library(org.Hs.eg.db)
library(EnhancedVolcano)
library(pheatmap)

## --- 1. Скачивание данных ---
## GSE152418 хранит raw counts как supplementary file (не FASTQ),
## поэтому используем его напрямую, без выравнивания.

gse_id <- "GSE152418"
dir.create("data", showWarnings = FALSE)

gse <- getGEO(gse_id, GSEMatrix = TRUE, destdir = "data")

## Скачиваем supplementary файл с матрицей counts
getGEOSuppFiles(gse_id, baseDir = "data")

## Путь к архиву с raw counts (проверьте имя файла после скачивания —
## GEO иногда меняет названия supplementary файлов)
supp_dir <- file.path("data", gse_id)
list.files(supp_dir)

## Пример распаковки, если файл .tar (типично для GEO с несколькими образцами)
tar_file <- list.files(supp_dir, pattern = "\\.tar$", full.names = TRUE)
if (length(tar_file) > 0) {
  untar(tar_file[1], exdir = supp_dir)
}

## --- 2. Загрузка counts matrix ---
## ПРИМЕЧАНИЕ: структура файлов у GSE152418 — по одному файлу на образец
## (формат: <GSM_ID>_<name>.counts.txt.gz). Ниже — сборка единой матрицы.

count_files <- list.files(supp_dir, pattern = "counts.*\\.txt\\.gz$",
                           full.names = TRUE, recursive = TRUE)

read_one <- function(f) {
  df <- read.table(f, header = FALSE, row.names = 1,
                    col.names = c("gene", basename(f)))
  df
}

count_list <- lapply(count_files, read_one)
counts_mat <- do.call(cbind, count_list)
colnames(counts_mat) <- gsub("_.*", "", basename(count_files))  # GSM ID как имя колонки

## --- 3. Метаданные образцов (условие: COVID vs здоровые) ---
pheno <- pData(gse[[1]])

## Адаптируйте под реальные названия колонок после просмотра pheno:
## head(pheno) — найдите колонку с группой (например "disease state:ch1")
coldata <- data.frame(
  row.names = rownames(pheno),
  condition = factor(pheno$`disease state:ch1`)  # проверить точное имя колонки!
)

## Выровнять порядок образцов между counts и coldata
coldata <- coldata[colnames(counts_mat), , drop = FALSE]

## --- 4. QC: базовая проверка ---
summary(colSums(counts_mat))          # глубина секвенирования по образцам
counts_mat <- counts_mat[rowSums(counts_mat) > 10, ]  # фильтр низкоэкспрессируемых генов

## --- 5. DESeq2: дифференциальная экспрессия ---
dds <- DESeqDataSetFromMatrix(
  countData = round(counts_mat),
  colData   = coldata,
  design    = ~ condition
)

dds <- DESeq(dds)
res <- results(dds, contrast = c("condition", "COVID-19", "healthy"))  # проверить точные уровни фактора
res <- res[order(res$padj), ]

write.csv(as.data.frame(res), "results_DESeq2.csv")

## --- 6. Визуализация ---

## PCA для проверки батч-эффектов / выбросов
vsd <- vst(dds, blind = FALSE)
plotPCA(vsd, intgroup = "condition")
ggplot2::ggsave("pca_plot.png", width = 6, height = 5)

## Volcano plot
EnhancedVolcano(res,
  lab = rownames(res),
  x = "log2FoldChange",
  y = "padj",
  pCutoff = 0.05,
  FCcutoff = 1,
  title = "COVID-19 vs Healthy (PBMC)"
)
ggplot2::ggsave("volcano_plot.png", width = 7, height = 6)

## Heatmap топ-30 дифференциально экспрессированных генов
top_genes <- rownames(res)[1:30]
pheatmap(assay(vsd)[top_genes, ],
         annotation_col = coldata,
         show_rownames = TRUE,
         filename = "heatmap_top30.png")

## --- 7. Функциональный анализ (GO enrichment) ---
sig_genes <- rownames(res)[which(res$padj < 0.05 & abs(res$log2FoldChange) > 1)]

ego <- enrichGO(
  gene          = sig_genes,
  OrgDb         = org.Hs.eg.db,
  keyType       = "SYMBOL",
  ont           = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05
)

write.csv(as.data.frame(ego), "GO_enrichment.csv")

## ============================================================
## Что нужно проверить/адаптировать под ваши данные:
## 1. Точное имя колонки в pheno с группой (disease state / condition)
## 2. Точные названия уровней фактора (COVID-19 / healthy — могут отличаться)
## 3. Формат файлов counts после скачивания (может измениться со стороны GEO)
## ============================================================
