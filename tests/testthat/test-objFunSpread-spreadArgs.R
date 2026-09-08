## The objective function calls SpaDES.tools::spread() Nreps times per fire year,
## so an argument that is silently ignored there is paid hundreds of times per
## evaluation. It passed `skipChecks = TRUE`, which spread() does not have --
## that belongs to spread3() -- so it went into `...` and did nothing, and every
## call re-validated the entire per-cell spreadProb vector: na.omit() copies it,
## inRange() scans it. Measured on synthetic landscapes with identical output:
## 1.9x the whole call at 1M cells, 3.4x at 4M, 5.1x at 9M.

test_that("objFunInner asks spread() to skip its checks, by the name spread() uses", {
  ## objFunInner() is where the per-year, per-replicate spread() call lives.
  ## `:::`: objFunInner() is internal, called from .objfunSpreadFit().
  src <- paste(deparse(fireSenseUtils:::objFunInner), collapse = "\n")
  expect_match(src, "quick = TRUE")
  expect_false(grepl("skipChecks", src))
})

test_that("`quick` is the name spread() acts on", {
  skip_if_not_installed("SpaDES.tools")
  args <- names(formals(SpaDES.tools::spread))
  expect_true("quick" %in% args)
})

test_that("spread() honours whichever of the two names is supplied", {
  ## This test replaces one that asserted `skipChecks` was NOT a formal of
  ## spread(). That was true when this file was written and is why the objective
  ## function's `skipChecks = TRUE` was silently ignored -- but SpaDES.tools then
  ## added `skipChecks` as an alias for `quick` (spread(), R/spread.R: the formal
  ## `skipChecks = quick`, then `if (isTRUE(skipChecks)) quick <- TRUE`), so the
  ## assertion became false and reddened this package's CI on development.
  ## Asserting the behaviour rather than the absence of an argument is what should
  ## have been tested: whichever name a caller uses, the checks get skipped.
  skip_if_not_installed("SpaDES.tools")
  skip_if_not_installed("terra")

  args <- names(formals(SpaDES.tools::spread))
  skip_if_not("skipChecks" %in% args,
              "this SpaDES.tools predates the skipChecks alias")

  n <- 100
  r <- terra::rast(nrows = n, ncols = n, xmin = 0, xmax = n, ymin = 0, ymax = n, vals = 1)
  loci <- c(2550L, 5050L)

  set.seed(42)
  viaQuick <- SpaDES.tools::spread(landscape = r, loci = loci, spreadProb = 0.23,
                                   returnIndices = TRUE, quick = TRUE)
  set.seed(42)
  viaAlias <- SpaDES.tools::spread(landscape = r, loci = loci, spreadProb = 0.23,
                                   returnIndices = TRUE, skipChecks = TRUE)

  expect_equal(viaQuick, viaAlias)
})

test_that("skipping the checks does not change what spread() returns", {
  skip_if_not_installed("SpaDES.tools")
  skip_if_not_installed("terra")
  n <- 300
  set.seed(1)
  r <- terra::rast(nrows = n, ncols = n, xmin = 0, xmax = n, ymin = 0, ymax = n, vals = 1)
  loci <- sample(terra::ncell(r), 20)
  sp <- runif(terra::ncell(r), 0.1, 0.3)
  ms <- pmax(2, round(rlnorm(20, log(200), 1)))
  set.seed(42); checked <- SpaDES.tools::spread(landscape = r, maxSize = ms, loci = loci,
                                                spreadProb = sp, returnIndices = TRUE,
                                                allowOverlap = FALSE, quick = FALSE)
  set.seed(42); quick <- SpaDES.tools::spread(landscape = r, maxSize = ms, loci = loci,
                                              spreadProb = sp, returnIndices = TRUE,
                                              allowOverlap = FALSE, quick = TRUE)
  expect_identical(checked, quick)
})
