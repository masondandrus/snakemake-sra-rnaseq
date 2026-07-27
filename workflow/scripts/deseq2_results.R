# Extract one contrast, apply LFC shrinkage, write a sorted results table
# and an MA plot.

log <- file(snakemake@log[[1]], open = "wt")
sink(log, type = "output")
sink(log, type = "message")

suppressPackageStartupMessages({
  library(DESeq2)
  library(ggplot2)
})

dds <- readRDS(snakemake@input[["dds"]])
contrast <- unlist(snakemake@params[["contrast"]])
alpha <- as.numeric(snakemake@params[["alpha"]])
shrink <- snakemake@params[["shrink"]]

stopifnot(length(contrast) == 3)
cat("contrast:", paste(contrast, collapse = " / "), "\n")

res <- results(dds, contrast = contrast, alpha = alpha)

if (shrink != "none") {
  # apeglm needs a coefficient name, so relevel to the denominator first.
  if (shrink == "apeglm") {
    factor_name <- contrast[1]
    dds[[factor_name]] <- relevel(dds[[factor_name]], ref = contrast[3])
    dds <- nbinomWaldTest(dds, quiet = TRUE)
    coef_name <- paste0(factor_name, "_", contrast[2], "_vs_", contrast[3])
    if (coef_name %in% resultsNames(dds)) {
      res <- lfcShrink(dds, coef = coef_name, type = "apeglm", res = res)
    } else {
      cat("apeglm coefficient not found; falling back to ashr\n")
      res <- lfcShrink(dds, contrast = contrast, type = "ashr", res = res)
    }
  } else {
    res <- lfcShrink(dds, contrast = contrast, type = shrink, res = res)
  }
}

df <- as.data.frame(res)
df$gene_id <- rownames(df)
df <- df[order(df$padj, na.last = NA), c("gene_id", setdiff(names(df), "gene_id"))]
write.table(df, snakemake@output[["tsv"]], sep = "\t", quote = FALSE, row.names = FALSE)

cat("significant at padj <", alpha, ":", sum(df$padj < alpha, na.rm = TRUE), "\n")

png(snakemake@output[["ma"]], width = 1400, height = 1000, res = 200)
plotMA(res, alpha = alpha, main = paste(contrast[2], "vs", contrast[3]))
dev.off()

sessionInfo()
