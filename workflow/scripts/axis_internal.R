# Turn a workflow output into a normalised evidence axis.
#
# Every axis in the convergence table, whatever its origin, is a TSV with a
# gene_id and/or gene_symbol column plus metadata columns prefixed with the
# axis name, and a sibling provenance JSON declaring which key to join on, the
# species, and whether the axis is usable. Keeping that contract uniform is
# what lets convergence.R take an arbitrary number of axes without knowing
# anything about where they came from.

log <- file(snakemake@log[[1]], open = "wt")
sink(log, type = "output")
sink(log, type = "message")

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(jsonlite)
})

axis    <- snakemake@params[["axis"]]
type    <- snakemake@params[["type"]]
species <- snakemake@params[["species"]]
opts    <- snakemake@params[["opts"]]

cat("axis:", axis, " type:", type, " species:", species, "\n")

prov <- list(
  name = axis, type = type, species = species,
  generated_utc = format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC")
)

if (type == "deseq2") {
  alpha <- as.numeric(opts$alpha %||% 0.05)
  res <- read_tsv(snakemake@input[["src"]], show_col_types = FALSE)
  tab <- res %>%
    filter(!is.na(padj), padj < alpha) %>%
    transmute(gene_id, gene_symbol,
              !!paste0(axis, "_padj") := padj,
              !!paste0(axis, "_log2fc") := log2FoldChange)
  prov$key <- "gene_id"
  prov$alpha <- alpha
  prov$skipped <- FALSE
  cat("DE genes at padj <", alpha, ":", nrow(tab), "\n")

} else if (type == "wgcna") {
  hubs <- read_tsv(snakemake@input[["src"]], show_col_types = FALSE)
  wp <- fromJSON(snakemake@input[["src_provenance"]])
  if (isTRUE(wp$skipped)) {
    # WGCNA declined to run (too few samples). Emit an empty but well-formed
    # axis so the DAG stays intact and convergence can report the axis as
    # unavailable rather than silently absent.
    tab <- tibble::tibble(gene_id = character(), gene_symbol = character())
    prov$skipped <- TRUE
    prov$skip_reason <- wp$reason
    cat("axis unavailable:", wp$reason, "\n")
  } else {
    tab <- hubs %>%
      transmute(gene_id,
                !!paste0(axis, "_module") := module,
                !!paste0(axis, "_kme") := kme,
                !!paste0(axis, "_module_trait_p") := module_trait_p)
    prov$skipped <- FALSE
    cat("hub genes:", nrow(tab), "\n")
  }
  prov$key <- "gene_id"

} else {
  stop(sprintf("axis_internal.R does not handle type '%s'", type))
}

prov$n_genes <- nrow(tab)
write_tsv(tab, snakemake@output[["tsv"]])
write_json(prov, snakemake@output[["provenance"]],
           auto_unbox = TRUE, pretty = TRUE, null = "null")
sessionInfo()
