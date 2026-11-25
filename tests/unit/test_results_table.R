source(file.path(lib_dir, "results_table.R"))

make <- function() {
  wald <- data.frame(
    baseMean = c(100, 5, 300, 0),
    log2FoldChange = c(2, 0.5, -4, NA),
    lfcSE = c(0.2, 1, 0.5, NA),
    stat = c(10, 0.5, -8, NA),
    pvalue = c(1e-20, 0.6, NA, NA),
    padj = c(1e-18, NA, NA, NA),
    row.names = c("g1", "g2", "g3", "g4")
  )
  final <- wald
  final$log2FoldChange <- c(1.8, 0.1, -3, NA)
  sym <- data.frame(gene_id = c("g1", "g2", "g3"), gene_symbol = c("A", "B", "C"))
  list(wald = wald, final = final, sym = sym)
}

test_that("every modeled gene keeps its row", {
  x <- make()
  out <- assemble_results(x$wald, x$final, x$sym)
  expect_equal(nrow(out), 4)
  expect_setequal(out$gene_id, c("g1", "g2", "g3", "g4"))
})

test_that("the reason for a missing adjusted p-value is recorded", {
  x <- make()
  out <- assemble_results(x$wald, x$final, x$sym)
  status <- setNames(out$padj_status, out$gene_id)
  expect_equal(unname(status["g1"]), "tested")
  expect_equal(unname(status["g2"]), "independent_filtering")
  expect_equal(unname(status["g3"]), "count_outlier")
  expect_equal(unname(status["g4"]), "no_counts")
})

test_that("the Wald statistic and both fold changes are kept", {
  x <- make()
  out <- assemble_results(x$wald, x$final, x$sym)
  g1 <- out[out$gene_id == "g1", ]
  expect_equal(g1$stat, 10)
  expect_equal(g1$log2FoldChange, 1.8)
  expect_equal(g1$log2FoldChange_unshrunk, 2)
})

test_that("tested genes sort first and a gene without a symbol gets an empty one", {
  x <- make()
  out <- assemble_results(x$wald, x$final, x$sym)
  expect_equal(out$gene_id[1], "g1")
  expect_equal(out$gene_symbol[out$gene_id == "g4"], "")
})
