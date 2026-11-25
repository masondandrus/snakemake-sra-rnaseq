# Join evidence axes into one table per contrast.
#
# For every gene in the results table and every axis, the script records a
# state: positive, negative, not_evaluable, or axis_skipped. It then tests,
# for each axis, whether genes positive on the reference axis are positive on
# that axis more often than other genes (one-sided Fisher's exact test over
# genes evaluable on both).
#
# n_evidence counts the axes on which a gene is positive. It is a way to rank
# candidates. The axes are not independent, and a gene shared by an eQTL axis
# and a GWAS axis is shared gene identity, not colocalisation.
#
# Ortholog tables are read from disk; this script makes no network request.

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

source(file.path(snakemake@scriptdir, "lib", "convergence_lib.R"))

`%||%` <- function(a, b) if (is.null(a)) b else a

res        <- read_tsv(snakemake@input[["results"]], show_col_types = FALSE)
axis_files <- unlist(snakemake@input[["axes"]])
prov_files <- unlist(snakemake@input[["provenance"]])
bg_files   <- unlist(snakemake@input[["backgrounds"]])
bg_axes    <- unlist(snakemake@params[["background_axes"]])
orth_files <- unlist(snakemake@input[["orthologs"]])
orth_species <- unlist(snakemake@params[["ortholog_species"]])
de_species <- snakemake@params[["species"]]
do_ortho   <- isTRUE(snakemake@params[["orthologs"]])
required   <- unlist(snakemake@params[["require"]])
reference  <- snakemake@params[["reference"]]

universe <- res %>% transmute(gene_id, gene_symbol)
cat("universe (genes in the results table):", nrow(universe), "\n")

orthologs <- list()
for (i in seq_along(orth_files)) {
  orthologs[[orth_species[i]]] <- read_tsv(orth_files[i], show_col_types = FALSE,
                                           col_types = cols(.default = col_character()))
}
backgrounds <- setNames(as.list(bg_files), bg_axes)

# The universe expressed in the identifier space of `species`.
bridge_for <- function(species) {
  if (identical(species, de_species)) {
    return(list(bridge = universe %>%
                  transmute(gene_id, mapped_id = gene_id,
                            mapped_symbol = toupper(gene_symbol), via_ortholog = FALSE),
                reason = NULL))
  }
  if (!do_ortho) return(list(bridge = NULL, reason = "ortholog mapping is disabled"))
  if (is.null(orthologs[[species]])) {
    return(list(bridge = NULL, reason = sprintf("no ortholog table for %s", species)))
  }
  list(bridge = orthologs[[species]] %>%
         filter(gene_id %in% universe$gene_id) %>%
         transmute(gene_id, mapped_id = ortholog_gene_id,
                   mapped_symbol = toupper(ortholog_symbol), via_ortholog = TRUE),
       reason = NULL)
}

evaluate_axis <- function(i) {
  p <- fromJSON(prov_files[i])
  name <- p$name
  key <- p$key %||% "gene_id"
  entry <- list(name = name, type = p$type, species = p$species, key = key)
  none <- tibble(gene_id = universe$gene_id, state = NA_character_, via_ortholog = FALSE)

  if (isTRUE(p$skipped)) {
    entry$status <- "skipped"
    entry$reason <- p$skip_reason %||% "skipped"
    return(list(entry = entry, genes = mutate(none, state = "axis_skipped"), meta = NULL))
  }
  b <- bridge_for(p$species)
  if (is.null(b$bridge)) {
    entry$status <- "not_evaluable"
    entry$reason <- b$reason
    return(list(entry = entry, genes = mutate(none, state = "not_evaluable"), meta = NULL))
  }

  tab <- read_tsv(axis_files[i], show_col_types = FALSE)
  bridge <- b$bridge
  if (key == "gene_id") {
    bridge$key <- bridge$mapped_id
    tab$key <- tab$gene_id
  } else {
    bridge$key <- bridge$mapped_symbol
    tab$key <- toupper(tab$gene_symbol)
  }
  background <- NULL
  if (!is.null(backgrounds[[name]])) {
    bg <- read_tsv(backgrounds[[name]], show_col_types = FALSE)
    background <- if (key == "gene_id") bg$gene_id else toupper(bg$gene_symbol)
  }

  state <- gene_states(universe$gene_id, bridge, unique(tab$key), background)
  via <- bridge %>% group_by(gene_id) %>%
    summarise(via_ortholog = any(via_ortholog), n_identifiers = n(), .groups = "drop")
  genes <- tibble(gene_id = universe$gene_id, state = unname(state)) %>%
    left_join(via, by = "gene_id") %>%
    mutate(via_ortholog = coalesce(via_ortholog, FALSE) & state != "not_evaluable")

  meta_cols <- setdiff(names(tab), c("gene_id", "gene_symbol", "key"))
  meta <- bridge %>% select(gene_id, key) %>%
    inner_join(tab %>% select(key, all_of(meta_cols)) %>% distinct(key, .keep_all = TRUE),
               by = "key") %>%
    distinct(gene_id, .keep_all = TRUE) %>% select(-key)

  entry$status <- "ok"
  entry$n_in_resource <- p$n_genes %||% nrow(tab)
  entry$background <- if (is.null(background)) "all genes with an identifier in the axis species" else "supplied by the axis"
  entry$n_evaluable <- sum(state != "not_evaluable")
  entry$n_positive <- sum(state == "positive")
  entry$n_one_to_many <- sum(via$n_identifiers > 1)
  list(entry = entry, genes = genes, meta = meta)
}

evaluated <- lapply(seq_along(axis_files), evaluate_axis)
axis_names <- map_chr(evaluated, ~ .x$entry$name)
names(evaluated) <- axis_names

conv <- universe
for (nm in axis_names) {
  g <- evaluated[[nm]]$genes
  conv[[nm]] <- g$state[match(conv$gene_id, g$gene_id)] == "positive"
  conv[[paste0(nm, "_state")]] <- g$state[match(conv$gene_id, g$gene_id)]
  conv[[paste0(nm, "_via_ortholog")]] <- g$via_ortholog[match(conv$gene_id, g$gene_id)]
  if (!is.null(evaluated[[nm]]$meta) && ncol(evaluated[[nm]]$meta) > 1) {
    conv <- left_join(conv, evaluated[[nm]]$meta, by = "gene_id")
  }
}

state_cols <- paste0(axis_names, "_state")
conv$n_evidence <- rowSums(as.matrix(conv[, axis_names, drop = FALSE]))
conv$n_evaluable <- rowSums(as.matrix(conv[, state_cols, drop = FALSE]) %in% c("positive", "negative") %>%
                              matrix(nrow = nrow(conv)))
conv$evidence <- apply(as.matrix(conv[, axis_names, drop = FALSE]), 1,
                       function(r) paste(axis_names[r], collapse = "|"))
conv$any_ortholog <- rowSums(as.matrix(conv[, paste0(axis_names, "_via_ortholog"), drop = FALSE])) > 0

# Overlap of every usable axis with the reference axis, over the whole universe.
overlap <- tibble()
ok_axes <- axis_names[map_chr(evaluated, ~ .x$entry$status) == "ok"]
if (!is.null(reference) && reference %in% ok_axes) {
  ref_state <- conv[[paste0(reference, "_state")]]
  overlap <- map_dfr(setdiff(ok_axes, reference), function(nm) {
    cbind(tibble(reference = reference, axis = nm),
          overlap_test(ref_state, conv[[paste0(nm, "_state")]]))
  })
  if (nrow(overlap) > 0) overlap$padj <- p.adjust(overlap$p_greater, method = "BH")
  cat("\noverlap with", reference, ":\n")
  print(as.data.frame(overlap))
} else {
  cat("\nno overlap test: reference axis", reference %||% "(none)", "is not usable\n")
}
if (nrow(overlap) == 0) {
  overlap <- tibble(reference = character(), axis = character(), n_universe = integer(),
                    n_reference = integer(), n_axis = integer(), n_both = integer(),
                    expected = numeric(), odds_ratio = numeric(), p_greater = numeric(),
                    padj = numeric())
}
write_tsv(overlap, snakemake@output[["overlap"]])

out <- conv
for (ax in intersect(required, axis_names)) out <- out[out[[ax]], , drop = FALSE]
out <- out %>% filter(n_evidence > 0) %>% arrange(desc(n_evidence), gene_symbol)
write_tsv(out, snakemake@output[["tsv"]])

cat("\ngenes with any evidence:", nrow(out), "\n")
if (nrow(out) > 0) print(count(out, evidence, sort = TRUE), n = 30)

sets <- axis_names[colSums(as.matrix(out[, axis_names, drop = FALSE])) > 0]
if (length(sets) >= 2 && requireNamespace("ComplexUpset", quietly = TRUE)) {
  plt <- ComplexUpset::upset(as.data.frame(out), sets, min_size = 1, width_ratio = 0.25,
                             name = "evidence axis") +
    labs(title = "Convergence of evidence across axes")
} else {
  bars <- tibble(axis = factor(axis_names, levels = axis_names),
                 genes = colSums(as.matrix(out[, axis_names, drop = FALSE])))
  plt <- ggplot(bars, aes(axis, genes)) +
    geom_col(fill = "#4C72B0", width = 0.6) +
    geom_text(aes(label = genes), vjust = -0.4) +
    labs(title = "Genes per evidence axis", x = NULL, y = "genes") +
    theme_minimal(base_size = 12)
}
ggsave(snakemake@output[["plot"]], plt, width = 9, height = 5, dpi = 150)

axes_prov <- lapply(evaluated, function(e) {
  entry <- e$entry
  entry$n_reported <- if (nrow(out) > 0) sum(out[[entry$name]]) else 0L
  entry
})
orth_prov <- lapply(seq_along(orth_files), function(i) {
  f <- sub("\\.tsv$", ".provenance.json", orth_files[i])
  c(list(species = orth_species[i], file = orth_files[i]),
    if (file.exists(f)) fromJSON(f) else list())
})
write_json(
  list(
    generated_utc = format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"),
    de_species = de_species,
    ortholog_mapping = do_ortho,
    ortholog_tables = orth_prov,
    universe = "genes in the results table (every modeled gene)",
    n_universe = nrow(universe),
    required_axes = as.list(required),
    reference_axis = reference,
    axes = unname(axes_prov),
    n_axes_usable = length(ok_axes),
    n_genes_any_evidence = nrow(out),
    max_evidence = if (nrow(out) > 0) max(out$n_evidence) else 0L,
    n_at_max = if (nrow(out) > 0) sum(out$n_evidence == max(out$n_evidence)) else 0L,
    caveat = paste(
      "n_evidence is a ranking aid, not a statistic; axes are correlated.",
      "The overlap test ignores expression level and gene length;",
      "eQTL/GWAS overlap is shared gene identity, not colocalisation."
    )
  ),
  snakemake@output[["provenance"]],
  auto_unbox = TRUE, pretty = TRUE, null = "null"
)

sessionInfo()
