# Rank-based gene set enrichment on the full DESeq2 result.
#
# Why GSEA rather than overrepresentation of the significant list: with a small
# design the padj < 0.05 set is often too short for ORA to say anything, while
# the ranked statistic still carries coordinated signal across a pathway.
#
# The ranking metric is the Wald statistic from the unshrunken test. Every
# modeled gene with a statistic and a symbol is ranked, whether or not it has
# an adjusted p-value. lib/ranking.R counts the genes left out, by reason, and
# the ranked list is written next to the result.

log <- file(snakemake@log[[1]], open = "wt")
sink(log, type = "output")
sink(log, type = "message")

suppressPackageStartupMessages({
  library(fgsea)
  library(msigdbr)
  library(dplyr)
  library(readr)
  library(ggplot2)
})

set.seed(1)

res <- read.delim(snakemake@input[["results"]], stringsAsFactors = FALSE)
species <- snakemake@params[["species"]]
min_size <- as.integer(snakemake@params[["min_size"]])
max_size <- as.integer(snakemake@params[["max_size"]])

source(file.path(snakemake@scriptdir, "lib", "ranking.R"))

ranked <- rank_genes(res)
cat("ranking:\n")
print(ranked$excluded)
write_tsv(ranked$ranks, snakemake@output[["ranks"]])
stats <- setNames(ranked$ranks$stat, ranked$ranks$gene_symbol)

# ---------------------------------------------------------------------------
# Gene sets. msigdbr changed its argument names between versions, so try the
# current signature first and fall back rather than pinning to one release.
# ---------------------------------------------------------------------------
get_sets <- function(coll, subcoll = NULL) {
  out <- try(
    if (is.null(subcoll)) {
      msigdbr(species = species, collection = coll)
    } else {
      msigdbr(species = species, collection = coll, subcollection = subcoll)
    },
    silent = TRUE
  )
  if (inherits(out, "try-error")) {
    out <- if (is.null(subcoll)) {
      msigdbr(species = species, category = coll)
    } else {
      msigdbr(species = species, category = coll, subcategory = subcoll)
    }
  }
  split(out$gene_symbol, out$gs_name)
}

collections <- list(
  hallmark = get_sets("H"),
  go_bp    = get_sets("C5", "GO:BP")
)

run_one <- function(pathways, label) {
  cat("\n---", label, "---\n")
  cat("pathways:", length(pathways), "\n")
  fg <- fgsea(
    pathways = pathways,
    stats = stats,
    minSize = min_size,
    maxSize = max_size,
    eps = 0,
    nproc = as.integer(snakemake@threads)
  )
  fg <- fg[order(fg$padj, -abs(fg$NES)), ]
  fg$collection <- label
  cat("padj < 0.05:", sum(fg$padj < 0.05, na.rm = TRUE), "\n")
  fg
}

all_res <- lapply(names(collections), function(n) run_one(collections[[n]], n))
combined <- do.call(rbind, all_res)

# leadingEdge is a list column; flatten so it survives a TSV round trip
combined$leadingEdge <- vapply(
  combined$leadingEdge, function(x) paste(x, collapse = ";"), character(1)
)

write_tsv(as.data.frame(combined), snakemake@output[["tsv"]])

# ---------------------------------------------------------------------------
# Plot: top pathways by adjusted p-value, signed by NES
# ---------------------------------------------------------------------------
top <- combined %>%
  as.data.frame() %>%
  filter(!is.na(padj)) %>%
  group_by(collection) %>%
  slice_min(padj, n = 12, with_ties = FALSE) %>%
  ungroup() %>%
  mutate(
    pathway = gsub("^(HALLMARK|GOBP)_", "", pathway),
    pathway = gsub("_", " ", pathway),
    pathway = tolower(pathway),
    pathway = reorder(pathway, NES)
  )

if (nrow(top) > 0) {
  p <- ggplot(top, aes(NES, pathway, fill = padj < 0.05)) +
    geom_col() +
    facet_wrap(~collection, scales = "free_y", ncol = 1) +
    scale_fill_manual(values = c("TRUE" = "#2c7fb8", "FALSE" = "grey70"),
                      name = "padj < 0.05") +
    labs(x = "Normalised enrichment score", y = NULL) +
    theme_bw(base_size = 9)
  ggsave(snakemake@output[["plot"]], p, width = 9, height = 10, dpi = 200)
} else {
  png(snakemake@output[["plot"]], width = 1200, height = 600, res = 150)
  plot.new(); text(0.5, 0.5, "No gene sets returned")
  dev.off()
}

sessionInfo()
