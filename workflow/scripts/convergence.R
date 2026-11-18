# Multi-evidence convergence over an arbitrary set of evidence axes.
#
# This script knows nothing about what the axes are. It reads N normalised axis
# files plus their provenance, bridges identifiers across species where an axis
# declares a different one, joins on gene, counts how many axes name each gene,
# and records how each hit was obtained.
#
# Adding a fifth or sixth line of evidence — differential ATAC peaks, CUT&RUN
# targets, a druggable-genome list, cell-type markers — is a config entry, not
# a code change.
#
# WHAT THE COUNT IS NOT
# ---------------------
# n_evidence is a triage device, not a statistic.
#   - The axes are not independent: expressed genes are likelier to have
#     detected eQTLs, and GWAS mapped genes skew toward gene-dense regions.
#   - eQTL and GWAS naming the same gene is shared gene identity, not
#     colocalisation. Establishing that needs coloc, SMR, or fine-mapping.
#   - Ortholog mapping is lossy in both directions, so every hit that depended
#     on one is flagged per axis and in aggregate.

log <- file(snakemake@log[[1]], open = "wt")
sink(log, type = "output")
sink(log, type = "message")

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(purrr)
  library(ggplot2)
  library(jsonlite)
})

`%||%` <- function(a, b) if (is.null(a) || length(a) == 0 || all(is.na(a))) b else a

axis_files <- unlist(snakemake@input[["axes"]])
prov_files <- unlist(snakemake@input[["provenance"]])
universe_f <- snakemake@input[["universe"]]

de_species    <- snakemake@params[["species"]]
ortho_enabled <- isTRUE(snakemake@params[["orthologs"]])
required      <- unlist(snakemake@params[["require"]] %||% character())

provs <- lapply(prov_files, fromJSON)
axis_names <- vapply(provs, function(p) p$name, character(1))
names(provs) <- axis_names
names(axis_files) <- axis_names

cat("axes declared:", paste(axis_names, collapse = ", "), "\n")
cat("data species:", de_species, "\n\n")

# Universe = every gene DESeq2 tested, so "no evidence" stays distinguishable
# from "never tested".
universe <- read_tsv(universe_f, show_col_types = FALSE) %>%
  transmute(gene_id, gene_symbol) %>%
  distinct(gene_id, .keep_all = TRUE)
cat("universe:", nrow(universe), "genes\n\n")

# ---------------------------------------------------------------------------
# Ortholog bridge, cached per target species so several same-species axes only
# pay for the mapping once.
# ---------------------------------------------------------------------------
ORTHO_CODE <- c("Mus musculus" = "mmusculus", "Rattus norvegicus" = "rnorvegicus",
                "Homo sapiens" = "hsapiens", "Danio rerio" = "drerio",
                "Macaca mulatta" = "mmulatta")
.ortho_cache <- new.env(parent = emptyenv())

empty_bridge <- function() {
  tibble::tibble(gene_id = character(), mapped_id = character(),
                 mapped_symbol = character(), via_ortholog = logical())
}

map_to <- function(target) {
  if (identical(target, de_species)) {
    return(universe %>% mutate(mapped_id = gene_id,
                               mapped_symbol = gene_symbol,
                               via_ortholog = FALSE))
  }
  key <- paste0(de_species, " -> ", target)
  cached <- .ortho_cache[[key]]
  if (!is.null(cached)) return(cached)

  if (!ortho_enabled) {
    cat("  ortholog mapping disabled; cannot reach", target, "\n")
    return(empty_bridge())
  }
  if (!all(c(de_species, target) %in% names(ORTHO_CODE))) {
    cat("  no g:Orth organism code for", de_species, "or", target, "\n")
    return(empty_bridge())
  }
  if (!requireNamespace("gprofiler2", quietly = TRUE)) {
    cat("  gprofiler2 unavailable; cannot map to", target, "\n")
    return(empty_bridge())
  }
  res <- try(gprofiler2::gorth(
    query = universe$gene_id,
    source_organism = ORTHO_CODE[[de_species]],
    target_organism = ORTHO_CODE[[target]],
    mthreshold = 1, filter_na = TRUE), silent = TRUE)
  if (inherits(res, "try-error") || is.null(res) || nrow(res) == 0) {
    cat("  g:Orth returned nothing for", key, "\n")
    return(empty_bridge())
  }
  out <- tibble::tibble(gene_id = res$input, mapped_id = res$ortholog_ensg,
                        mapped_symbol = res$ortholog_name, via_ortholog = TRUE) %>%
    filter(!is.na(mapped_id), mapped_id != "N/A") %>%
    distinct(gene_id, .keep_all = TRUE)
  cat("  ", key, ": mapped ", nrow(out), " of ", nrow(universe), "\n", sep = "")
  assign(key, out, envir = .ortho_cache)
  out
}

# ---------------------------------------------------------------------------
# Join one axis onto the universe
# ---------------------------------------------------------------------------
no_hits <- function(name) {
  tibble::tibble(gene_id = character()) %>%
    mutate(!!name := logical(), !!paste0(name, "_via_ortholog") := logical())
}

join_axis <- function(name) {
  p <- provs[[name]]

  if (isTRUE(p$skipped)) {
    cat(name, ": unavailable -", p$skip_reason %||% "skipped", "\n")
    return(no_hits(name))
  }

  tab <- read_tsv(axis_files[[name]], show_col_types = FALSE)
  if (nrow(tab) == 0) {
    cat(name, ": no genes in resource\n")
    return(no_hits(name))
  }

  key <- p$key %||% "gene_id"
  bridge <- map_to(p$species %||% de_species)
  if (nrow(bridge) == 0) {
    cat(name, ": no identifier bridge to", p$species %||% de_species, "\n")
    return(no_hits(name))
  }

  meta <- setdiff(names(tab), c("gene_id", "gene_symbol"))

  if (key == "gene_symbol") {
    joined <- bridge %>%
      filter(!is.na(mapped_symbol)) %>%
      mutate(.k = toupper(mapped_symbol)) %>%
      inner_join(tab %>% mutate(.k = toupper(gene_symbol)) %>%
                   select(-any_of(c("gene_symbol", "gene_id"))), by = ".k")
  } else {
    joined <- bridge %>%
      inner_join(tab %>% select(-any_of("gene_symbol")),
                 by = c("mapped_id" = "gene_id"))
  }

  if (nrow(joined) == 0) {
    cat(name, ": 0 genes matched\n")
    return(no_hits(name))
  }

  out <- joined %>%
    mutate(!!name := TRUE, !!paste0(name, "_via_ortholog") := via_ortholog) %>%
    select(gene_id, all_of(name), all_of(paste0(name, "_via_ortholog")),
           any_of(meta)) %>%
    distinct(gene_id, .keep_all = TRUE)

  cat(name, ": ", nrow(out), " genes",
      if (any(joined$via_ortholog)) " (via orthologs)" else "", "\n", sep = "")
  out
}

joined <- reduce(map(axis_names, join_axis),
                 function(acc, x) left_join(acc, x, by = "gene_id"),
                 .init = universe)

hit_cols <- axis_names[axis_names %in% names(joined)]
flag_cols <- paste0(axis_names, "_via_ortholog")
flag_cols <- flag_cols[flag_cols %in% names(joined)]

if (length(hit_cols) == 0) {
  stop("No axis contributed any column; nothing to converge on.")
}

conv <- joined %>%
  mutate(across(all_of(hit_cols), ~ replace_na(.x, FALSE)),
         across(all_of(flag_cols), ~ replace_na(.x, FALSE))) %>%
  mutate(n_evidence = rowSums(across(all_of(hit_cols))),
         evidence = apply(across(all_of(hit_cols)), 1,
                          function(r) paste(hit_cols[as.logical(r)], collapse = "|")),
         any_ortholog = if (length(flag_cols))
           rowSums(across(all_of(flag_cols))) > 0 else FALSE)

# Optional hard filter: a gene must hit these axes to appear at all.
if (length(required)) {
  absent <- setdiff(required, hit_cols)
  if (length(absent)) {
    cat("\nWARNING: required axes contributed nothing:",
        paste(absent, collapse = ", "), "\n")
  }
  keep <- intersect(required, hit_cols)
  if (length(keep)) {
    before <- sum(conv$n_evidence > 0)
    conv <- conv %>% filter(if_all(all_of(keep), ~ .x))
    cat("\nrequired [", paste(keep, collapse = ", "), "]: ",
        nrow(conv), " retained of ", before, " with any evidence\n", sep = "")
  }
}

conv <- conv %>% filter(n_evidence > 0) %>% arrange(desc(n_evidence), gene_symbol)
write_tsv(conv, snakemake@output[["table"]])

cat("\nevidence combinations:\n")
print(as.data.frame(count(conv, n_evidence, evidence) %>%
                      arrange(desc(n_evidence), desc(n))))

available <- axis_names[vapply(axis_names,
                               function(n) !isTRUE(provs[[n]]$skipped), logical(1))]
maxn <- if (nrow(conv)) max(conv$n_evidence) else 0
top <- filter(conv, n_evidence == maxn)
cat("\naxes with data: ", length(available), " of ", length(axis_names),
    " (", paste(available, collapse = ", "), ")\n", sep = "")
cat("genes at maximum evidence (", maxn, "): ", nrow(top), "\n", sep = "")
if (nrow(top)) {
  print(as.data.frame(head(select(top, gene_id, gene_symbol, evidence,
                                  any_of(c("DEG_log2fc", "DEG_padj"))), 25)))
}

# ---------------------------------------------------------------------------
# Plot: how the evidence sets intersect
# ---------------------------------------------------------------------------
png(snakemake@output[["plot"]], width = 1500, height = 900, res = 150)
mat <- conv %>% select(all_of(hit_cols)) %>% as.data.frame()
usable <- hit_cols[colSums(mat) > 0]
if (length(usable) >= 2 && requireNamespace("ComplexUpset", quietly = TRUE)) {
  print(ComplexUpset::upset(mat, usable, name = "evidence axis",
                            width_ratio = 0.25, min_size = 1) +
          ggtitle("Multi-evidence convergence"))
} else if (length(usable) >= 1) {
  barplot(colSums(mat)[usable], col = "#2c7fb8", ylab = "genes",
          main = "Genes per evidence axis")
} else {
  plot.new(); text(0.5, 0.5, "no axis returned any genes")
}
dev.off()

write_json(
  list(
    generated_utc = format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"),
    de_species = de_species,
    ortholog_mapping = ortho_enabled,
    required_axes = as.list(required),
    axes = lapply(axis_names, function(n) list(
      name = n, type = provs[[n]]$type, species = provs[[n]]$species,
      key = provs[[n]]$key, skipped = isTRUE(provs[[n]]$skipped),
      skip_reason = provs[[n]]$skip_reason,
      n_in_resource = provs[[n]]$n_genes,
      n_matched = if (n %in% names(conv)) sum(conv[[n]]) else 0)),
    n_axes_available = length(available),
    n_genes_any_evidence = nrow(conv),
    max_evidence = maxn, n_at_max = nrow(top),
    caveat = paste("n_evidence is a prioritisation heuristic, not a statistic;",
                   "axes are correlated, and eQTL/GWAS overlap is shared gene",
                   "identity, not colocalisation")
  ),
  snakemake@output[["provenance"]], auto_unbox = TRUE, pretty = TRUE, null = "null"
)

sessionInfo()
