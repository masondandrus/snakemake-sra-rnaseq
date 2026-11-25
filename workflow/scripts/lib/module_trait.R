# Association between each module eigengene and one contrast, fitted with the
# same design formula and the same numerator and denominator that DESeq2 uses.
#
# eigengenes: data frame, one row per sample (row names), one ME* column per module
# coldata:    sample table with the design variables, row names = samples
# design:     formula with no left-hand side, e.g. ~ sex + condition
# contrast:   c(factor, numerator, denominator)

module_trait <- function(eigengenes, coldata, design, contrast) {
  stopifnot(length(contrast) == 3)
  f <- contrast[1]; num <- contrast[2]; den <- contrast[3]
  if (!f %in% names(coldata)) {
    stop(sprintf("contrast factor '%s' is not a column of the sample table", f))
  }
  coldata <- as.data.frame(coldata)[rownames(eigengenes), , drop = FALSE]
  levels_present <- unique(as.character(coldata[[f]]))
  for (lv in c(num, den)) {
    if (!lv %in% levels_present) {
      stop(sprintf("level '%s' of factor '%s' is not in the sample table", lv, f))
    }
  }
  coldata[[f]] <- relevel(factor(coldata[[f]]), ref = den)
  coef_name <- paste0(f, num)
  modules <- setdiff(grep("^ME", colnames(eigengenes), value = TRUE), "MEgrey")
  if (length(modules) == 0) {
    return(data.frame(module = character(), estimate = numeric(), se = numeric(),
                      t = numeric(), df = numeric(), p = numeric(), padj = numeric()))
  }
  rows <- lapply(modules, function(m) {
    d <- coldata
    d$.me <- eigengenes[[m]]
    fit <- lm(update(design, .me ~ .), data = d)
    co <- coef(summary(fit))
    if (!coef_name %in% rownames(co)) {
      stop(sprintf("coefficient '%s' is not estimable under design %s",
                   coef_name, paste(deparse(design), collapse = "")))
    }
    data.frame(module = sub("^ME", "", m),
               estimate = co[coef_name, "Estimate"],
               se = co[coef_name, "Std. Error"],
               t = co[coef_name, "t value"],
               df = fit$df.residual,
               p = co[coef_name, "Pr(>|t|)"])
  })
  out <- do.call(rbind, rows)
  out$padj <- p.adjust(out$p, method = "BH")
  out <- out[order(out$p), ]
  rownames(out) <- NULL
  out
}

select_hubs <- function(modules, trait, kme_threshold, module_alpha) {
  keep_modules <- trait$module[!is.na(trait$padj) & trait$padj < module_alpha]
  hubs <- modules[modules$module %in% keep_modules &
                    !is.na(modules$kme) & modules$kme >= kme_threshold, ]
  hubs <- merge(hubs, trait[, c("module", "estimate", "p", "padj")], by = "module")
  hubs <- hubs[order(hubs$padj, -hubs$kme),
               c("gene_id", "module", "kme", "estimate", "p", "padj")]
  rownames(hubs) <- NULL
  hubs
}
