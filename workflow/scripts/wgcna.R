# Weighted gene co-expression network: modules, eigengenes, and each gene's
# correlation with its own module eigengene (kME).
#
# The network is built once from all samples and knows nothing about any
# contrast. The association between modules and a contrast is fitted per
# contrast in axis_internal.R, with the DESeq2 design.
#
# Below `min_samples` the script writes empty tables, records skipped = TRUE
# with the reason, and exits 0.

log <- file(snakemake@log[[1]], open = "wt")
sink(log, type = "output")
sink(log, type = "message")

suppressPackageStartupMessages({
  library(DESeq2)
  library(WGCNA)
  library(readr)
  library(jsonlite)
})

options(stringsAsFactors = FALSE)
enableWGCNAThreads(nThreads = as.integer(snakemake@threads))

vsd         <- readRDS(snakemake@input[["vst"]])
min_samples <- as.integer(snakemake@params[["min_samples"]])
n_top       <- as.integer(snakemake@params[["n_top_genes"]])
min_mod     <- as.integer(snakemake@params[["min_module_size"]])
r2_cut      <- as.numeric(snakemake@params[["scale_free_r2"]])
soft_power  <- snakemake@params[["power"]]

n <- ncol(vsd)
cat("samples:", n, " minimum for WGCNA:", min_samples, "\n")

utc_now <- function() format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC")

if (n < min_samples) {
  reason <- sprintf("n = %d samples, below the minimum of %d", n, min_samples)
  cat(reason, "\n")
  write_tsv(data.frame(gene_id = character(), module = character(), kme = numeric()),
            snakemake@output[["modules"]])
  write_tsv(data.frame(sample = character()), snakemake@output[["eigengenes"]])
  write_tsv(data.frame(Power = numeric(), SFT.R.sq = numeric(), mean.k. = numeric()),
            snakemake@output[["soft_threshold"]])
  png(snakemake@output[["plot"]], width = 1400, height = 700, res = 150)
  plot.new()
  text(0.5, 0.55, "WGCNA skipped", cex = 1.6)
  text(0.5, 0.40, reason, cex = 0.9)
  dev.off()
  write_json(list(skipped = TRUE, reason = reason, n_samples = n,
                  min_samples = min_samples, generated_utc = utc_now()),
             snakemake@output[["provenance"]], auto_unbox = TRUE, pretty = TRUE)
  sessionInfo()
  quit(save = "no", status = 0)
}

mat <- assay(vsd)
vars <- apply(mat, 1, var)
keep <- head(order(vars, decreasing = TRUE), min(n_top, nrow(mat)))
datExpr <- t(mat[keep, ])
cat("genes into network:", ncol(datExpr), "\n")

gsg <- goodSamplesGenes(datExpr, verbose = 0)
if (!gsg$allOK) {
  cat("dropping", sum(!gsg$goodGenes), "genes and",
      sum(!gsg$goodSamples), "samples failing QC\n")
  datExpr <- datExpr[gsg$goodSamples, gsg$goodGenes]
}

# Soft threshold: the lowest power whose scale-free fit reaches r2_cut. When no
# power does, the fallback follows the table in the WGCNA FAQ for signed
# networks, by number of samples.
faq_power <- function(n_samples) {
  if (n_samples < 20) 18 else if (n_samples < 30) 16 else if (n_samples < 40) 14 else 12
}

sft <- pickSoftThreshold(datExpr, powerVector = c(1:10, seq(12, 30, 2)),
                         networkType = "signed", RsquaredCut = r2_cut, verbose = 0)
write_tsv(sft$fitIndices, snakemake@output[["soft_threshold"]])

if (identical(soft_power, "auto") || is.null(soft_power)) {
  power <- sft$powerEstimate
  power_source <- sprintf("lowest power with scale-free R^2 >= %s", r2_cut)
  if (is.na(power)) {
    power <- faq_power(nrow(datExpr))
    power_source <- sprintf("no power reached R^2 >= %s; WGCNA FAQ value for %d samples",
                            r2_cut, nrow(datExpr))
  }
} else {
  power <- as.integer(soft_power)
  power_source <- "set in the configuration"
}
cat("soft threshold power:", power, "(", power_source, ")\n")

net <- blockwiseModules(
  datExpr, power = power, networkType = "signed",
  TOMType = "signed", minModuleSize = min_mod,
  reassignThreshold = 0, mergeCutHeight = 0.25,
  numericLabels = FALSE, pamRespectsDendro = FALSE,
  saveTOMs = FALSE, verbose = 0, maxBlockSize = 20000
)

MEs <- orderMEs(net$MEs)
kme_all <- signedKME(datExpr, MEs)
own <- cbind(seq_len(ncol(datExpr)), match(paste0("kME", net$colors), colnames(kme_all)))
modules <- data.frame(gene_id = colnames(datExpr), module = unname(net$colors),
                      kme = kme_all[own])
write_tsv(modules, snakemake@output[["modules"]])
write_tsv(data.frame(sample = rownames(datExpr), MEs, check.names = FALSE),
          snakemake@output[["eigengenes"]])

cat("modules found (grey = unassigned):\n")
print(table(net$colors))

png(snakemake@output[["plot"]], width = 1600, height = 900, res = 150)
plotDendroAndColors(net$dendrograms[[1]], net$colors[net$blockGenes[[1]]],
                    "Module", dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05, main = "Gene clustering")
dev.off()

write_json(
  list(skipped = FALSE, n_samples = nrow(datExpr), n_genes = ncol(datExpr),
       network_type = "signed", soft_power = power, soft_power_source = power_source,
       scale_free_r2_cut = r2_cut,
       n_modules = length(setdiff(unique(net$colors), "grey")),
       n_unassigned = sum(net$colors == "grey"),
       min_module_size = min_mod, generated_utc = utc_now()),
  snakemake@output[["provenance"]], auto_unbox = TRUE, pretty = TRUE
)

sessionInfo()
