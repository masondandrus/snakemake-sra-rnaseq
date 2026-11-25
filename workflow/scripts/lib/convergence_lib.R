# Gene-level state on one evidence axis, and the overlap test between two axes.
#
# A gene is in one of three states on an axis:
#   positive       evaluable and named by the axis
#   negative       evaluable and not named by the axis
#   not_evaluable  the axis could not have named it (no identifier in the
#                  axis's species, or outside the axis's background)
# When the axis itself was skipped, every gene gets axis_skipped.

# bridge:     data frame gene_id, key  (key = the gene's identifier in the
#             axis's identifier space; several rows per gene when a gene has
#             several orthologs; genes without an identifier are absent)
# axis_keys:  identifiers the axis names
# background: identifiers the axis could have named, or NULL for "any"
gene_states <- function(universe_ids, bridge, axis_keys, background = NULL) {
  bridge <- bridge[!is.na(bridge$key) & bridge$key != "", , drop = FALSE]
  in_bg <- if (is.null(background)) rep(TRUE, nrow(bridge)) else bridge$key %in% background
  evaluable <- unique(bridge$gene_id[in_bg])
  positive <- unique(bridge$gene_id[in_bg & bridge$key %in% axis_keys])
  state <- ifelse(universe_ids %in% positive, "positive",
             ifelse(universe_ids %in% evaluable, "negative", "not_evaluable"))
  setNames(state, universe_ids)
}

# One-sided Fisher's exact test: are genes positive on `ref` positive on
# `axis` more often than other genes? Only genes evaluable on both count.
overlap_test <- function(ref_state, axis_state) {
  both <- ref_state %in% c("positive", "negative") & axis_state %in% c("positive", "negative")
  r <- ref_state[both] == "positive"
  a <- axis_state[both] == "positive"
  n <- sum(both)
  tab <- matrix(c(sum(r & a), sum(r & !a), sum(!r & a), sum(!r & !a)), nrow = 2)
  out <- data.frame(n_universe = n, n_reference = sum(r), n_axis = sum(a),
                    n_both = sum(r & a),
                    expected = if (n > 0) sum(r) * sum(a) / n else NA_real_,
                    odds_ratio = NA_real_, p_greater = NA_real_)
  if (n > 0 && sum(r) > 0 && sum(a) > 0 && sum(!r) > 0 && sum(!a) > 0) {
    ft <- fisher.test(tab, alternative = "greater")
    out$odds_ratio <- unname(ft$estimate)
    out$p_greater <- ft$p.value
  }
  out
}
