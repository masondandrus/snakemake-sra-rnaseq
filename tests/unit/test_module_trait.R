source(file.path(lib_dir, "module_trait.R"))

set.seed(7)
n <- 24
coldata <- data.frame(
  dose = rep(c("vehicle", "low", "high"), each = 8),
  sex = rep(c("F", "M"), times = 12),
  row.names = paste0("s", seq_len(n))
)
eig <- data.frame(
  MEhigh = ifelse(coldata$dose == "high", 1, 0) + rnorm(n, 0, 0.2),
  MElow = ifelse(coldata$dose == "low", 1, 0) + rnorm(n, 0, 0.2),
  MEsex = ifelse(coldata$sex == "F", 1, 0) + rnorm(n, 0, 0.2),
  MEnoise = rnorm(n),
  MEgrey = rnorm(n),
  row.names = rownames(coldata)
)

gap <- function(me, f, a, b) mean(eig[[me]][coldata[[f]] == a]) - mean(eig[[me]][coldata[[f]] == b])

test_that("a three-level factor is tested for the requested pair of levels", {
  hi <- module_trait(eig, coldata, ~ sex + dose, c("dose", "high", "vehicle"))
  lo <- module_trait(eig, coldata, ~ sex + dose, c("dose", "low", "vehicle"))
  expect_true(hi$padj[hi$module == "high"] < 0.001)
  expect_true(lo$padj[lo$module == "low"] < 0.001)
  # The design is balanced, so each estimate equals the difference in group means.
  expect_equal(hi$estimate[hi$module == "high"], gap("MEhigh", "dose", "high", "vehicle"))
  expect_equal(hi$estimate[hi$module == "low"], gap("MElow", "dose", "high", "vehicle"))
  expect_equal(lo$estimate[lo$module == "low"], gap("MElow", "dose", "low", "vehicle"))
  expect_equal(lo$estimate[lo$module == "high"], gap("MEhigh", "dose", "low", "vehicle"))
})

test_that("swapping numerator and denominator flips the sign", {
  a <- module_trait(eig, coldata, ~ sex + dose, c("dose", "high", "vehicle"))
  b <- module_trait(eig, coldata, ~ sex + dose, c("dose", "vehicle", "high"))
  expect_equal(a$estimate[a$module == "high"], -b$estimate[b$module == "high"])
})

test_that("a module driven by another factor is tested on its own factor", {
  hi <- module_trait(eig, coldata, ~ sex + dose, c("dose", "high", "vehicle"))
  sx <- module_trait(eig, coldata, ~ sex + dose, c("sex", "F", "M"))
  expect_equal(hi$estimate[hi$module == "sex"], gap("MEsex", "dose", "high", "vehicle"))
  expect_true(abs(hi$estimate[hi$module == "sex"]) < 0.3)
  expect_true(sx$padj[sx$module == "sex"] < 0.001)
  expect_equal(sx$estimate[sx$module == "sex"], gap("MEsex", "sex", "F", "M"))
})

test_that("grey is never tested and p-values are corrected across modules", {
  hi <- module_trait(eig, coldata, ~ sex + dose, c("dose", "high", "vehicle"))
  expect_false("grey" %in% hi$module)
  expect_equal(hi$padj, p.adjust(hi$p, "BH"))
})

test_that("a missing factor or level stops with a message", {
  expect_error(module_trait(eig, coldata, ~ dose, c("batch", "a", "b")), "not a column")
  expect_error(module_trait(eig, coldata, ~ dose, c("dose", "medium", "vehicle")),
               "not in the sample table")
})

test_that("hubs come only from associated modules and only above the kME threshold", {
  modules <- data.frame(
    gene_id = paste0("g", 1:8),
    module = c("high", "high", "high", "noise", "noise", "sex", "grey", "low"),
    kme = c(0.95, 0.85, 0.5, 0.99, 0.9, 0.9, 0.9, 0.9)
  )
  trait <- data.frame(module = c("high", "low", "sex", "noise"),
                      estimate = c(1, 0.1, 0, 0.2), p = c(1e-9, 0.2, 0.5, 0.6),
                      padj = c(4e-9, 0.4, 0.6, 0.6))
  hubs <- select_hubs(modules, trait, kme_threshold = 0.8, module_alpha = 0.05)
  expect_setequal(hubs$gene_id, c("g1", "g2"))
})
