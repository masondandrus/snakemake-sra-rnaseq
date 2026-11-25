# Assemble the per-contrast results table. Every modeled gene keeps its row.
#
# wald:  DESeq2 results() output before shrinkage (carries stat, pvalue, padj)
# final: the same contrast after shrinkage, or `wald` itself when shrinkage is off
# sym:   data frame with gene_id and gene_symbol

padj_status <- function(base_mean, pvalue, padj) {
  ifelse(!is.na(padj), "tested",
    ifelse(base_mean == 0, "no_counts",
      ifelse(is.na(pvalue), "count_outlier", "independent_filtering")))
}

assemble_results <- function(wald, final, sym) {
  stopifnot(identical(rownames(wald), rownames(final)))
  df <- data.frame(
    gene_id = rownames(wald),
    baseMean = wald$baseMean,
    log2FoldChange = final$log2FoldChange,
    lfcSE = final$lfcSE,
    stat = wald$stat,
    pvalue = wald$pvalue,
    padj = wald$padj,
    log2FoldChange_unshrunk = wald$log2FoldChange,
    stringsAsFactors = FALSE
  )
  df$padj_status <- padj_status(df$baseMean, df$pvalue, df$padj)
  df$gene_symbol <- sym$gene_symbol[match(df$gene_id, sym$gene_id)]
  df$gene_symbol[is.na(df$gene_symbol)] <- ""
  df <- df[order(is.na(df$padj), df$padj, -abs(df$stat)), ]
  rownames(df) <- NULL
  df[, c("gene_id", "gene_symbol", "baseMean", "log2FoldChange", "lfcSE",
         "stat", "pvalue", "padj", "log2FoldChange_unshrunk", "padj_status")]
}
