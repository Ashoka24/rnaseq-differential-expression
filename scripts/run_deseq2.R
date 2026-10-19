#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(DESeq2)
  library(ggplot2)
})

dir.create("results/figures", showWarnings = FALSE, recursive = TRUE)
count_url <- "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE116139&file=GSE116139_BulkRNAseqCounts.txt.gz&format=file"
input <- file.path(tempdir(), "GSE116139_BulkRNAseqCounts.txt.gz")
download.file(count_url, input, mode = "wb", quiet = TRUE)

meta <- data.frame(
  sample = c("GSM3211196","GSM3211194","GSM3211169","GSM3211171","GSM3211174","GSM3211175",
             "GSM3211177","GSM3211178","GSM3211186","GSM3211185","GSM3211180","GSM3211179",
             "GSM3211183","GSM3211182","GSM3211189","GSM3211190","GSM3211193","GSM3211191"),
  subject = factor(c(19,19,34,34,40,40,42,42,43,43,47,47,49,49,53,53,56,56)),
  tissue = factor(c("mucosa","blood","mucosa","blood","mucosa","blood","mucosa","blood",
                    "mucosa","blood","mucosa","blood","mucosa","blood","mucosa","blood",
                    "mucosa","blood"), levels = c("blood","mucosa"))
)
rownames(meta) <- meta$sample

raw <- read.delim(gzfile(input), check.names = FALSE, comment.char = "", stringsAsFactors = FALSE)
raw_names <- colnames(raw)
gsm <- c("GSM3211169","GSM3211170","GSM3211171","GSM3211172","GSM3211173","GSM3211174",
"GSM3211175","GSM3211176","GSM3211177","GSM3211178","GSM3211179","GSM3211180",
"GSM3211181","GSM3211182","GSM3211183","GSM3211184","GSM3211185","GSM3211186",
"GSM3211187","GSM3211188","GSM3211189","GSM3211190","GSM3211191","GSM3211192",
"GSM3211193","GSM3211194","GSM3211195","GSM3211196")
libraries <- c("lib10819","lib10820","lib10821","lib10823","lib11819","lib11820","lib11821",
"lib11822","lib11823","lib11824","lib12778","lib12779","lib12780","lib12783","lib12784",
"lib12785","lib12786","lib12787","lib12788","lib12791","lib12792","lib12793","lib12794",
"lib12795","lib12796","lib7438","lib7439","lib7440")
if (length(gsm) != 28 || length(libraries) != 28 || anyDuplicated(gsm) || anyDuplicated(libraries))
  stop("Expected a one-to-one mapping of 28 GEO samples to count-matrix libraries")
if (!all(libraries %in% raw_names)) stop("Mapped libraries absent from the count matrix")

mapping <- data.frame(gsm = gsm, library_name = libraries, stringsAsFactors = FALSE)
write.table(mapping, "results/library_mapping.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
meta$library_id <- libraries[match(meta$sample, gsm)]
gene_col <- intersect(c("gene_id","Geneid","gene","Gene"), raw_names)[1]
if (is.na(gene_col)) gene_col <- raw_names[1]
library_idx <- match(meta$library_id, raw_names)
if (anyNA(library_idx)) stop("Mapped libraries absent from count matrix")
counts <- raw[, c(match(gene_col, raw_names), library_idx), drop = FALSE]
colnames(counts)[1] <- "gene_id"
colnames(counts)[-1] <- meta$sample
counts$gene_id <- as.character(counts$gene_id)
counts <- counts[counts$gene_id != "" & !duplicated(counts$gene_id), , drop = FALSE]
rownames(counts) <- counts$gene_id
counts$gene_id <- NULL
counts[] <- lapply(counts, function(x) suppressWarnings(as.numeric(as.character(x))))
if (anyNA(as.matrix(counts)) || any(as.matrix(counts) < 0)) stop("Invalid count matrix values")
counts <- round(as.matrix(counts))
meta <- meta[colnames(counts), , drop = FALSE]
stopifnot(all(rownames(meta) == colnames(counts)))
keep <- rowSums(counts >= 10) >= 9
counts <- counts[keep, , drop = FALSE]

dds <- DESeqDataSetFromMatrix(countData = counts, colData = meta, design = ~ subject + tissue)
dds <- DESeq(dds)
res <- as.data.frame(results(dds, contrast = c("tissue","mucosa","blood")))
res$gene_id <- rownames(res)
res <- res[order(res$padj, -abs(res$log2FoldChange)),]
write.table(res, "results/differential_expression.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
norm <- counts(dds, normalized = TRUE)
write.table(data.frame(gene_id=rownames(norm), norm, check.names=FALSE), "results/normalized_counts.tsv", sep="\t", quote=FALSE, row.names=FALSE)
sig <- subset(res, !is.na(padj) & padj < 0.05 & abs(log2FoldChange) >= 1)
vsd <- vst(dds, blind=FALSE)
ggsave("results/figures/01_pca.svg", plotPCA(vsd, intgroup=c("tissue","subject")), width=7, height=5)
res_plot <- res[!is.na(res$padj),]
res_plot$significant <- ifelse(res_plot$padj < 0.05 & abs(res_plot$log2FoldChange) >= 1, "Significant", "Not significant")
p <- ggplot(res_plot, aes(log2FoldChange, -log10(padj), color=significant)) + geom_point(alpha=.5, size=1) + theme_minimal() + labs(title="Mucosa vs blood: differential expression", x="log2 fold change", y="-log10 adjusted p-value") + theme(legend.title=element_blank())
ggsave("results/figures/02_volcano_plot.svg", p, width=7, height=5)
up_count <- sum(sig$log2FoldChange > 0); down_count <- sum(sig$log2FoldChange < 0)
metrics <- data.frame(metric=c("samples","subjects","genes_after_filtering","significant_genes","significant_up","significant_down"), value=c(nrow(meta),length(unique(meta$subject)),nrow(counts),nrow(sig),up_count,down_count))
write.table(metrics, "results/metrics.tsv", sep="\t", quote=FALSE, row.names=FALSE)
