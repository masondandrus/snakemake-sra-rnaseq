source(file.path(lib_dir, "ranking.R"))

res <- data.frame(
  gene_id = c("g1", "g2", "g3", "g4", "g5", "g6", "g7"),
  gene_symbol = c("A", "B", "C", "", "D", "D", "E"),
  stat = c(98.4, 0.5, -8, 3, 2, -6, NA),
  pvalue = c(0, 0.6, NA, 0.003, 0.05, 1e-9, NA),
  padj = c(0, NA, NA, 0.01, 0.1, 1e-7, NA),
  stringsAsFactors = FALSE
)

test_that("a gene with a p-value of zero stays in the ranking", {
  r <- rank_genes(res)
  expect_true("A" %in% r$ranks$gene_symbol)
  expect_equal(r$ranks$gene_symbol[1], "A")
})

test_that("a gene without an adjusted p-value is ranked when it has a p-value", {
  r <- rank_genes(res)
  expect_true("B" %in% r$ranks$gene_symbol)
})

test_that("count outliers, genes without a statistic or a symbol are counted and left out", {
  r <- rank_genes(res)
  expect_false("C" %in% r$ranks$gene_symbol)
  expect_equal(unname(r$excluded["count_outlier"]), 1)
  expect_equal(unname(r$excluded["no_statistic"]), 1)
  expect_equal(unname(r$excluded["no_symbol"]), 1)
})

test_that("a duplicated symbol keeps the row with the larger absolute statistic", {
  r <- rank_genes(res)
  expect_equal(sum(r$ranks$gene_symbol == "D"), 1)
  expect_equal(r$ranks$stat[r$ranks$gene_symbol == "D"], -6)
  expect_equal(unname(r$excluded["duplicate_symbol"]), 1)
})

test_that("the counts add up and the order is decreasing", {
  r <- rank_genes(res)
  e <- r$excluded
  expect_equal(unname(e["total"]),
               unname(e["no_statistic"] + e["count_outlier"] + e["no_symbol"] +
                        e["duplicate_symbol"] + e["ranked"]))
  expect_false(is.unsorted(rev(r$ranks$stat)))
})
