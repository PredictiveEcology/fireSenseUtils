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

test_that("`quick` is spread()'s argument and `skipChecks` is not, which is why this matters", {
  skip_if_not_installed("SpaDES.tools")
  args <- names(formals(SpaDES.tools::spread))
  expect_true("quick" %in% args)
  expect_false("skipChecks" %in% args)   # would land in `...` and be ignored
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
