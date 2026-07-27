# Build the DESeq2 object, filter low-count genes, fit the model, and emit
# a VST-based PCA as a first-pass diagnostic for batch or outlier structure.

log <- file(snakemake@log[[1]], open = "wt")
sink(log, type = "output")
sink(log, type = "message")

suppressPackageStartupMessages({
  library(DESeq2)
  library(ggplot2)
})

txi <- readRDS(snakemake@input[["txi"]])
coldata <- read.delim(snakemake@input[["samples"]], stringsAsFactors = TRUE)
rownames(coldata) <- coldata$sample
coldata <- coldata[colnames(txi$counts), , drop = FALSE]

design <- as.formula(snakemake@params[["design"]])
dds <- DESeqDataSetFromTximport(txi, colData = coldata, design = design)

# Pre-filtering: independent of the design, so it does not bias the test.
min_count <- as.integer(snakemake@params[["min_count"]])
min_samples <- as.integer(snakemake@params[["min_samples"]])
keep <- rowSums(counts(dds) >= min_count) >= min_samples
cat("genes before filtering:", nrow(dds), "\n")
dds <- dds[keep, ]
cat("genes after filtering:", nrow(dds), "\n")

dds <- DESeq(dds)
saveRDS(dds, snakemake@output[["dds"]])

vsd <- vst(dds, blind = TRUE)
saveRDS(vsd, snakemake@output[["vst"]])

group <- snakemake@params[["pca_group"]]
pca <- plotPCA(vsd, intgroup = group, returnData = TRUE)
pct <- round(100 * attr(pca, "percentVar"))

p <- ggplot(pca, aes(PC1, PC2, colour = .data[[group]], label = name)) +
  geom_point(size = 3) +
  geom_text(vjust = -1, size = 3, show.legend = FALSE) +
  labs(
    x = paste0("PC1: ", pct[1], "% variance"),
    y = paste0("PC2: ", pct[2], "% variance"),
    colour = group
  ) +
  theme_bw()

ggsave(snakemake@output[["pca"]], p, width = 7, height = 5, dpi = 200)
sessionInfo()
