# Extracted from test-objFunSpread-pruneAbove.R:35

# setup ------------------------------------------------------------------------
library(testthat)
test_env <- simulate_test_env(package = "fireSenseUtils", path = "..")
attach(test_env, warn.conflicts = FALSE)

# test -------------------------------------------------------------------------
src <- paste(deparse(fireSenseUtils:::.objfunSpreadFit), collapse = "\n")
expect_match(src, "min\\(thresh \\* numYrsDone, pruneAbove\\)")
