# Weighted gene co-expression network analysis: modules, module-trait
# correlation, and hub genes.
#
# SAMPLE SIZE
# -----------
# WGCNA needs enough samples to estimate a correlation matrix that means
# anything. The package authors put the floor at 15 and recommend 20+; below
# that, the topology overlap is dominated by noise and you get modules that
# look plausible and reproduce nothing.
#
# So this script refuses rather than guesses. Below `min_samples` it writes an
# empty hub table with skipped = TRUE in the provenance and exits 0. The
# convergence step then runs on the remaining axes and records that the hub
# axis was unavailable. A missing axis is honest; a fabricated one is not.

log <- file(snakemake@log[[1]], open = "wt")
sink(log, type = "output")
sink(log, type = "message")

suppressPackageStartupMessages({
  library(DESeq2)
  library(WGCNA)
  library(dplyr)
  library(readr)
  library(jsonlite)
})

options(stringsAsFactors = FALSE)
enableWGCNAThreads(nThreads = as.integer(snakemake@threads))

vsd         <- readRDS(snakemake@input[["vst"]])
dds         <- readRDS(snakemake@input[["dds"]])
min_samples <- as.integer(snakemake@params[["min_samples"]])
n_top       <- as.integer(snakemake@params[["n_top_genes"]])
min_mod     <- as.integer(snakemake@params[["min_module_size"]])
kme_cut     <- as.numeric(snakemake@params[["kme_threshold"]])
trait_col   <- snakemake@params[["trait"]]
soft_power  <- snakemake@params[["power"]]

n <- ncol(vsd)
cat("samples:", n, " minimum for WGCNA:", min_samples, "\n")

write_empty <- function(reason) {
  write_tsv(
    tibble::tibble(gene_id = character(), module = character(),
                   kme = numeric(), module_trait_cor = numeric(),
                   module_trait_p = numeric()),
    snakemake@output[["hubs"]]
  )
  write_tsv(
    tibble::tibble(gene_id = character(), module = character()),
    snakemake@output[["modules"]]
  )
  png(snakemake@output[["plot"]], width = 1400, height = 700, res = 150)
  plot.new()
  text(0.5, 0.55, "WGCNA skipped", cex = 1.6)
  text(0.5, 0.40, reason, cex = 0.9)
  dev.off()
  write_json(
    list(skipped = TRUE, reason = reason, n_samples = n,
         min_samples = min_samples,
         generated_utc = format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC")),
    snakemake@output[["provenance"]], auto_unbox = TRUE, pretty = TRUE
  )
}

if (n < min_samples) {
  reason <- sprintf("n = %d samples, below the minimum of %d", n, min_samples)
  cat("\n", reason, "\n", sep = "")
  cat("Writing empty hub table; the convergence step will drop this axis.\n")
  write_empty(reason)
  sessionInfo()
  quit(save = "no", status = 0)
}

# ---------------------------------------------------------------------------
# Input matrix: VST counts, restricted to the most variable genes. WGCNA on the
# full transcriptome is slow and mostly noise; the top few thousand by variance
# is the standard compromise.
# ---------------------------------------------------------------------------
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

# ---------------------------------------------------------------------------
# Soft threshold. Pick the lowest power reaching scale-free R^2 >= 0.8, or take
# the configured value.
# ---------------------------------------------------------------------------
if (identical(soft_power, "auto") || is.null(soft_power)) {
  sft <- pickSoftThreshold(datExpr, powerVector = c(1:10, seq(12, 20, 2)),
                           networkType = "signed", verbose = 0)
  power <- sft$powerEstimate
  if (is.na(power)) {
    power <- if (ncol(datExpr) < 20) 18 else 12
    cat("no power reached R^2 >= 0.8; falling back to", power, "\n")
  }
} else {
  power <- as.integer(soft_power)
}
cat("soft threshold power:", power, "\n")

net <- blockwiseModules(
  datExpr, power = power, networkType = "signed",
  TOMType = "signed", minModuleSize = min_mod,
  reassignThreshold = 0, mergeCutHeight = 0.25,
  numericLabels = FALSE, pamRespectsDendro = FALSE,
  saveTOMs = FALSE, verbose = 0, maxBlockSize = 20000
)

modules <- tibble::tibble(gene_id = colnames(datExpr), module = net$colors)
write_tsv(modules, snakemake@output[["modules"]])
cat("modules found:", length(unique(net$colors)), "\n")
print(table(net$colors))

# ---------------------------------------------------------------------------
# Module-trait correlation. The trait is the experimental factor, coerced to a
# numeric contrast. Hubs are only meaningful relative to a module that tracks
# something; a hub of an unrelated module is just a well-connected gene.
# ---------------------------------------------------------------------------
coldata <- as.data.frame(colData(dds))[rownames(datExpr), , drop = FALSE]
trait_raw <- coldata[[trait_col]]
if (is.null(trait_raw)) {
  stop(sprintf("trait column '%s' not found in colData", trait_col))
}
trait <- as.numeric(as.factor(trait_raw))
cat("trait:", trait_col, " levels:", paste(levels(as.factor(trait_raw)), collapse = ", "), "\n")

MEs <- orderMEs(net$MEs)
mt_cor <- cor(MEs, trait, use = "p")
mt_p <- corPvalueStudent(mt_cor, nrow(datExpr))
mt <- tibble::tibble(
  module = sub("^ME", "", rownames(mt_cor)),
  module_trait_cor = as.numeric(mt_cor),
  module_trait_p = as.numeric(mt_p)
) %>% arrange(module_trait_p)
cat("\nmodule-trait correlation:\n")
print(as.data.frame(mt))

# kME: correlation of each gene with its own module eigengene.
kme_all <- signedKME(datExpr, MEs)
colnames(kme_all) <- sub("^kME", "", colnames(kme_all))

hubs <- modules %>%
  rowwise() %>%
  mutate(kme = if (module %in% colnames(kme_all))
    kme_all[gene_id, module] else NA_real_) %>%
  ungroup() %>%
  left_join(mt, by = "module") %>%
  filter(module != "grey", !is.na(kme), kme >= kme_cut) %>%
  arrange(module_trait_p, desc(kme))

write_tsv(hubs, snakemake@output[["hubs"]])
cat("\nhub genes at kME >=", kme_cut, ":", nrow(hubs), "\n")

png(snakemake@output[["plot"]], width = 1600, height = 900, res = 150)
par(mfrow = c(1, 2), mar = c(6, 5, 3, 2))
plotDendroAndColors(net$dendrograms[[1]],
                    net$colors[net$blockGenes[[1]]],
                    "Module", dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Gene clustering")
barplot(mt$module_trait_cor, names.arg = mt$module, las = 2,
        col = ifelse(mt$module_trait_p < 0.05, "#2c7fb8", "grey70"),
        ylab = "correlation with trait", main = paste("Module vs", trait_col))
dev.off()

write_json(
  list(skipped = FALSE, n_samples = n, n_genes = ncol(datExpr),
       soft_power = power, n_modules = length(unique(net$colors)),
       kme_threshold = kme_cut, n_hubs = nrow(hubs), trait = trait_col,
       generated_utc = format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC")),
  snakemake@output[["provenance"]], auto_unbox = TRUE, pretty = TRUE
)

sessionInfo()
