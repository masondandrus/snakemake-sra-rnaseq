suppressPackageStartupMessages(library(testthat))
args <- commandArgs(trailingOnly = TRUE)
repo <- if (length(args)) args[1] else "."
lib_dir <- normalizePath(file.path(repo, "workflow", "scripts", "lib"))
res <- test_dir(file.path(repo, "tests", "unit"), env = environment(),
                reporter = "summary", stop_on_failure = FALSE)
df <- as.data.frame(res)
cat(sprintf("\ntests: %d  passed expectations: %d  failed: %d  errors: %d\n",
            nrow(df), sum(df$passed), sum(df$failed), sum(df$error)))
quit(status = if (sum(df$failed) + sum(df$error) > 0) 1 else 0)
