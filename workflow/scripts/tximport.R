# Summarise salmon transcript quantifications to gene level.
#
# countsFromAbundance = "lengthScaledTPM" produces counts corrected for both
# library size and average transcript length, which keeps DESeq2's negative
# binomial assumptions intact while removing length bias across samples.

log <- file(snakemake@log[[1]], open = "wt")
sink(log, type = "output")
sink(log, type = "message")

suppressPackageStartupMessages({
  library(tximport)
  library(readr)
})

samples <- read.delim(snakemake@input[["samples"]], stringsAsFactors = FALSE)
quants <- snakemake@input[["quants"]]

# Order the quant files to match the sample sheet rather than trusting glob order.
names(quants) <- basename(dirname(quants))
stopifnot(all(samples$run %in% names(quants)))
quants <- quants[samples$run]

t2g <- read.delim(snakemake@input[["tx2gene"]], stringsAsFactors = FALSE)
t2g <- t2g[, c("tx_id", "gene_id")]

txi <- tximport(
  quants,
  type = "salmon",
  tx2gene = t2g,
  ignoreTxVersion = TRUE,
  dropInfReps = TRUE,
  countsFromAbundance = "lengthScaledTPM"
)

colnames(txi$counts) <- samples$sample
colnames(txi$abundance) <- samples$sample

write_counts <- function(mat, path) {
  df <- data.frame(gene_id = rownames(mat), mat, check.names = FALSE)
  readr::write_tsv(df, path)
}

write_counts(round(txi$counts, 3), snakemake@output[["counts"]])
write_counts(round(txi$abundance, 3), snakemake@output[["tpm"]])
saveRDS(txi, snakemake@output[["rds"]])

cat("genes:", nrow(txi$counts), " samples:", ncol(txi$counts), "\n")
sessionInfo()
