# Ranked gene list for enrichment, keyed on gene symbol and ordered by the
# Wald statistic. Returns the ranking and a count of every gene left out,
# by reason.
#
# A gene that DESeq2 flagged as a count outlier keeps a statistic in the
# results table but has no p-value, because the statistic rests on one
# extreme count. It is left out of the ranking.

rank_genes <- function(res) {
  stopifnot(all(c("gene_id", "gene_symbol", "stat", "pvalue") %in% names(res)))
  no_stat <- is.na(res$stat) | !is.finite(res$stat)
  outlier <- !no_stat & is.na(res$pvalue)
  no_symbol <- !no_stat & !outlier & (is.na(res$gene_symbol) | res$gene_symbol == "")
  r <- res[!no_stat & !outlier & !no_symbol, c("gene_id", "gene_symbol", "stat")]
  r <- r[order(-abs(r$stat), r$gene_id), ]
  dup <- duplicated(r$gene_symbol)
  r <- r[!dup, ]
  r <- r[order(-r$stat, r$gene_symbol), ]
  rownames(r) <- NULL
  list(
    ranks = r,
    excluded = c(
      total = nrow(res),
      no_statistic = sum(no_stat),
      count_outlier = sum(outlier),
      no_symbol = sum(no_symbol),
      duplicate_symbol = sum(dup),
      ranked = nrow(r)
    )
  )
}
