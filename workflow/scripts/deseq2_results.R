# Extract one contrast, apply LFC shrinkage, and write the results table and
# an MA plot. The table holds every modeled gene: genes without an adjusted
# p-value stay in it, and padj_status says why the value is missing.

log <- file(snakemake@log[[1]], open = "wt")
sink(log, type = "output")
sink(log, type = "message")

suppressPackageStartupMessages({
  library(DESeq2)
  library(ggplot2)
})

source(file.path(snakemake@scriptdir, "lib", "results_table.R"))

dds <- readRDS(snakemake@input[["dds"]])
t2g <- read.delim(snakemake@input[["tx2gene"]], stringsAsFactors = FALSE)
sym <- t2g[!duplicated(t2g$gene_id), c("gene_id", "gene_symbol")]
contrast <- unlist(snakemake@params[["contrast"]])
alpha <- as.numeric(snakemake@params[["alpha"]])
shrink <- snakemake@params[["shrink"]]

stopifnot(length(contrast) == 3)
cat("contrast:", paste(contrast, collapse = " / "), "\n")

wald <- results(dds, contrast = contrast, alpha = alpha)
final <- wald

if (shrink != "none") {
  # apeglm needs a coefficient name, so relevel to the denominator first.
  if (shrink == "apeglm") {
    factor_name <- contrast[1]
    dds[[factor_name]] <- relevel(dds[[factor_name]], ref = contrast[3])
    dds <- nbinomWaldTest(dds, quiet = TRUE)
    coef_name <- paste0(factor_name, "_", contrast[2], "_vs_", contrast[3])
    if (coef_name %in% resultsNames(dds)) {
      final <- lfcShrink(dds, coef = coef_name, type = "apeglm", res = wald)
    } else {
      cat("apeglm coefficient not found; falling back to ashr\n")
      final <- lfcShrink(dds, contrast = contrast, type = "ashr", res = wald)
      shrink <- "ashr (fallback from apeglm)"
    }
  } else {
    final <- lfcShrink(dds, contrast = contrast, type = shrink, res = wald)
  }
}
cat("shrinkage:", shrink, "\n")

df <- assemble_results(wald, final, sym)
stopifnot(nrow(df) == nrow(dds))
write.table(df, snakemake@output[["tsv"]], sep = "\t", quote = FALSE, row.names = FALSE)

cat("genes in the model:", nrow(df), "\n")
print(table(padj_status = df$padj_status))
cat("significant at padj <", alpha, ":", sum(df$padj < alpha, na.rm = TRUE), "\n")

png(snakemake@output[["ma"]], width = 1400, height = 1000, res = 200)
plotMA(final, alpha = alpha, main = paste(contrast[2], "vs", contrast[3]))
dev.off()

sessionInfo()
