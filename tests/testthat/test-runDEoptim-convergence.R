## runDEoptim() passes the early-stop settings (`convergence`) to clusters::DEoptimIterative(). They are in
## the fit's cache key, so a stopped fit rerun with a stricter rule is a new call (which resumes from
## clusters' per-generation cache); NULL leaves the key what it was before the argument existed.

test_that("convergence reaches clusters::DEoptimIterative(), NULL by default", {
  seen <- new.env(); seen$calls <- 0L
  testthat::local_mocked_bindings(
    clusterSetup = function(...) list(itermax = 5, trace = FALSE, strategy = 2L, NP = 40L),
    DEoptimIterative = function(fn, lower, upper, control, ..., convergence = "unset") {
      seen$convergence <- convergence
      list()
    },
    .package = "clusters")
  testthat::local_mocked_bindings(termsInDEoptim = function(...) invisible(NULL))
  withr::local_options(reproducible.useCache = FALSE)
  rule <- list(window = 100L, minGenerations = 600L, se = 0.5, sdChange = 0.05)
  callRunDEoptimPlotEvery(withr::local_tempdir(), convergence = rule)
  expect_identical(seen$convergence, rule)
  callRunDEoptimPlotEvery(withr::local_tempdir())
  expect_null(seen$convergence)
})

test_that("convergence is in the fit's cache key; NULL is not", {
  seen <- new.env(); seen$calls <- 0L
  mockDEoptimPlotEvery(seen)
  withr::local_options(reproducible.useCache = TRUE)
  cp <- withr::local_tempdir()
  callRunDEoptimPlotEvery(cp)
  callRunDEoptimPlotEvery(cp, convergence = NULL)        # a cache hit: same key as before the argument
  expect_identical(seen$calls, 1L)
  callRunDEoptimPlotEvery(cp, convergence = list(minGenerations = 600L))   # a stricter rule is a new fit
  expect_identical(seen$calls, 2L)
  callRunDEoptimPlotEvery(cp, convergence = list(minGenerations = 600L))   # ... once
  expect_identical(seen$calls, 2L)
})
