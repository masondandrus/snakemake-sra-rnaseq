# Ortholog table from the reference species to one target species, fetched
# once from g:Profiler and saved with the service version. A failed request
# stops the job.

log <- file(snakemake@log[[1]], open = "wt")
sink(log, type = "output")
sink(log, type = "message")

suppressPackageStartupMessages({
  library(gprofiler2)
  library(dplyr)
  library(readr)
  library(jsonlite)
})

src <- snakemake@wildcards[["src"]]
dst <- snakemake@wildcards[["dst"]]
t2g <- read.delim(snakemake@input[["tx2gene"]], stringsAsFactors = FALSE)
genes <- unique(t2g$gene_id)
cat("genes to map:", length(genes), " ", src, "->", dst, "\n")

chunks <- split(genes, ceiling(seq_along(genes) / 5000))
mapped <- lapply(chunks, function(g) {
  r <- gorth(g, source_organism = src, target_organism = dst,
             mthreshold = Inf, filter_na = TRUE)
  if (is.null(r)) return(NULL)
  r[, c("input", "ortholog_ensg", "ortholog_name")]
})
mapped <- bind_rows(mapped)
if (nrow(mapped) == 0) stop("g:Profiler returned no orthologs for ", src, " -> ", dst)

out <- mapped %>%
  transmute(gene_id = input, ortholog_gene_id = ortholog_ensg, ortholog_symbol = ortholog_name) %>%
  filter(!is.na(ortholog_gene_id), ortholog_gene_id != "", ortholog_gene_id != "N/A") %>%
  distinct()
write_tsv(out, snakemake@output[["tsv"]])

version <- tryCatch(get_version_info(organism = dst),
                    error = function(e) list(gprofiler_version = NA))
per_gene <- count(out, gene_id)
write_json(
  list(source_organism = src, target_organism = dst,
       service = "g:Profiler g:Orth",
       gprofiler_version = version$gprofiler_version,
       retrieved_utc = format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"),
       n_genes_queried = length(genes),
       n_genes_mapped = nrow(per_gene),
       n_genes_one_to_many = sum(per_gene$n > 1),
       n_pairs = nrow(out)),
  snakemake@output[["provenance"]], auto_unbox = TRUE, pretty = TRUE, null = "null"
)
cat("genes mapped:", nrow(per_gene), " one-to-many:", sum(per_gene$n > 1), "\n")
sessionInfo()
