# Turn a workflow output into a normalised evidence axis.
#
# Every axis in the convergence table, whatever its origin, is a TSV with a
# gene_id and/or gene_symbol column plus metadata columns prefixed with the
# axis name, and a sibling provenance JSON declaring which key to join on, the
# species, and whether the axis is usable.

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

`%||%` <- function(a, b) if (is.null(a)) b else a

cat("axis:", axis, " type:", type, " species:", species, "\n")

prov <- list(
  name = axis, type = type, species = species, key = "gene_id",
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
  background <- res %>% filter(!is.na(padj)) %>% select(gene_id)
  prov$alpha <- alpha
  prov$skipped <- FALSE
  prov$n_eligible <- nrow(background)
  cat("DE genes at padj <", alpha, ":", nrow(tab), "\n")

} else if (type == "wgcna") {
  suppressPackageStartupMessages(library(DESeq2))
  source(file.path(snakemake@scriptdir, "lib", "module_trait.R"))

  wp <- fromJSON(snakemake@input[["src_provenance"]])
  contrast <- unlist(snakemake@params[["contrast"]])
  design <- as.formula(snakemake@params[["design"]])
  kme_cut <- as.numeric(snakemake@params[["kme_threshold"]])
  module_alpha <- as.numeric(snakemake@params[["module_alpha"]])
  empty_trait <- data.frame(module = character(), estimate = numeric(), se = numeric(),
                            t = numeric(), df = numeric(), p = numeric(), padj = numeric())

  if (isTRUE(wp$skipped)) {
    tab <- tibble::tibble(gene_id = character())
    background <- tibble::tibble(gene_id = character())
    trait <- empty_trait
    prov$skipped <- TRUE
    prov$skip_reason <- wp$reason
    cat("axis unavailable:", wp$reason, "\n")
  } else {
    modules <- as.data.frame(read_tsv(snakemake@input[["src"]], show_col_types = FALSE))
    eig <- as.data.frame(read_tsv(snakemake@input[["eigengenes"]], show_col_types = FALSE))
    rownames(eig) <- eig$sample
    eig$sample <- NULL
    coldata <- as.data.frame(colData(readRDS(snakemake@input[["dds"]])))

    background <- tibble::tibble(gene_id = modules$gene_id)
    trait <- module_trait(eig, coldata, design, contrast)
    hubs <- select_hubs(modules, trait, kme_cut, module_alpha)
    cat("module association with", paste(contrast, collapse = " / "),
        "under", paste(deparse(design), collapse = ""), ":\n")
    print(trait)

    tab <- hubs %>%
      transmute(gene_id,
                !!paste0(axis, "_module") := module,
                !!paste0(axis, "_kme") := kme,
                !!paste0(axis, "_module_estimate") := estimate,
                !!paste0(axis, "_module_padj") := padj)
    prov$skipped <- FALSE
    prov$design <- paste(deparse(design), collapse = "")
    prov$contrast <- as.list(contrast)
    prov$kme_threshold <- kme_cut
    prov$module_alpha <- module_alpha
    prov$n_modules_tested <- nrow(trait)
    prov$n_modules_associated <- sum(trait$padj < module_alpha)
    prov$n_eligible <- nrow(background)
    cat("modules associated at padj <", module_alpha, ":",
        prov$n_modules_associated, "of", nrow(trait), "\n")
    cat("hub genes at kME >=", kme_cut, ":", nrow(tab), "\n")
  }
  write_tsv(trait, snakemake@output[["module_trait"]])

} else {
  stop(sprintf("axis_internal.R does not handle type '%s'", type))
}

prov$n_genes <- nrow(tab)
write_tsv(tab, snakemake@output[["tsv"]])
write_tsv(background, snakemake@output[["background"]])
write_json(prov, snakemake@output[["provenance"]],
           auto_unbox = TRUE, pretty = TRUE, null = "null")
sessionInfo()
