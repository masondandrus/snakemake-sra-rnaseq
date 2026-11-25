source(file.path(lib_dir, "convergence_lib.R"))

universe <- paste0("g", 1:10)

test_that("genes without an identifier in the axis's species are not evaluable", {
  bridge <- data.frame(gene_id = c("g1", "g2", "g3"), key = c("H1", "H2", NA))
  s <- gene_states(universe, bridge, axis_keys = "H1")
  expect_equal(unname(s["g1"]), "positive")
  expect_equal(unname(s["g2"]), "negative")
  expect_equal(unname(s["g3"]), "not_evaluable")
  expect_equal(unname(s["g9"]), "not_evaluable")
})

test_that("a gene with several orthologs is positive when any of them is on the axis", {
  bridge <- data.frame(gene_id = c("g1", "g1", "g2"), key = c("H1", "H7", "H2"))
  s <- gene_states(universe, bridge, axis_keys = "H7")
  expect_equal(unname(s["g1"]), "positive")
  expect_equal(unname(s["g2"]), "negative")
})

test_that("a background limits which genes are evaluable", {
  bridge <- data.frame(gene_id = c("g1", "g2", "g3"), key = c("H1", "H2", "H3"))
  s <- gene_states(universe, bridge, axis_keys = c("H1", "H3"), background = c("H1", "H2"))
  expect_equal(unname(s[c("g1", "g2", "g3")]), c("positive", "negative", "not_evaluable"))
})

test_that("the overlap test counts only genes evaluable on both axes", {
  ref  <- c(rep("positive", 10), rep("negative", 90), rep("not_evaluable", 50))
  axis <- c(rep("positive", 8), rep("negative", 2), rep("positive", 9), rep("negative", 81),
            rep("positive", 50))
  o <- overlap_test(ref, axis)
  expect_equal(o$n_universe, 100)
  expect_equal(o$n_reference, 10)
  expect_equal(o$n_axis, 17)
  expect_equal(o$n_both, 8)
  expect_equal(o$expected, 1.7)
  expected_p <- fisher.test(matrix(c(8, 2, 9, 81), 2), alternative = "greater")$p.value
  expect_equal(o$p_greater, expected_p)
  expect_true(o$p_greater < 1e-5)
})

test_that("an overlap at the chance rate is not significant, and empty margins give no test", {
  ref  <- c(rep("positive", 20), rep("negative", 180))
  axis <- c(rep("positive", 2), rep("negative", 18), rep("positive", 18), rep("negative", 162))
  o <- overlap_test(ref, axis)
  expect_equal(o$odds_ratio, 1, tolerance = 0.05)
  expect_true(o$p_greater > 0.5)
  none <- overlap_test(rep("negative", 5), c("positive", rep("negative", 4)))
  expect_true(is.na(none$p_greater))
})
