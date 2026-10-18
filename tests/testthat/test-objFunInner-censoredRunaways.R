## censoredRunaways: a simulated fire that reaches its own buffer's edge ring is right-censored at the size it
## burned. It keeps burning inside its buffer (no early stop), so the size terms score a buffer-bounded size as
## the observed fires are; the per-fire t likelihood fits the replicates' location and scale by censored-normal
## maximum likelihood, a runaway contributing P(size >= burned size). Before, a runaway was dropped from the
## likelihood (with a share penalty) and scored as runawayBufferMultiple x its buffer in the size terms, which made
## the objective prefer yearSpreadSD = 0 whatever the true year effect (known-sigma test, ELF 14.3, 2026-10-17).
## Fixture: helper-runawayFixture.R.

test_that("censoredNormalFit(): no censoring gives the sample mean and sd; censoring points are 'at least this'", {
  x <- c(3, 4, 5, 6, 7)
  expect_equal(censoredNormalFit(x, rep(FALSE, 5)), list(mu = 5, s = sd(x)))
  cens <- censoredNormalFit(x, c(FALSE, FALSE, FALSE, TRUE, TRUE))
  expect_gt(cens$mu, 5)
  ## known parameters, 30% censored at 11: the fit recovers them, the naive mean does not
  set.seed(1)
  z <- rnorm(4000, 10, 2)
  isCens <- z > 11
  xc <- ifelse(isCens, 11, z)
  fit <- censoredNormalFit(xc, isCens)
  expect_equal(fit$mu, 10, tolerance = 0.02)
  expect_equal(fit$s, 2, tolerance = 0.05)
  expect_lt(mean(xc), 9.8)
})

test_that("sizeLikTCensored() is sizeLikT() without censoring; censored runaways give a big observed fire more likelihood", {
  N <- c(20, 25, 30, 22, 28, 24, 26, 27, 23, 29)
  expect_equal(sizeLikTCensored(26, N, rep(FALSE, 10), 5), sizeLikT(26, N, 5))
  ran <- N >= 28 # the three largest replicates reached the ring: the fire is at least that big in them
  dropped <- mean(!ran) * sizeLikT(60, N[!ran], 5) # the share-and-drop rule
  expect_gt(sizeLikTCensored(60, N, ran, 5), dropped)
})

test_that("with censoredRunaways a runaway keeps burning, the size terms score what it burned, the likelihood changes", {
  skip_if_not_installed("SpaDES.tools")
  seen <- new.env()
  orig <- SpaDES.tools::spreadCpp
  local_mocked_bindings(
    spreadCpp = function(...) {
      seen$stop <- c(seen$stop, !is.null(list(...)$stopCells))
      orig(...)
    },
    .package = "SpaDES.tools")
  sp <- c(0.22, 0.3) # 10 x 10 buffers, fires 5 cells from the edge: some replicates reach the ring, not all
  cens <- runawayObjective(sp, seed = 3, sizeLik = "t", censorRunaways = TRUE, censoredRunaways = TRUE)
  expect_false(any(seen$stop))                      # no early stop: the burned size is the censoring point
  expect_gt(cens$runaways, 0)
  expect_lt(cens$runaways, cens$nSims)
  ## penaliseRunaways = FALSE draws the same fires (no stop either); only the likelihood differs
  off <- runawayObjective(sp, seed = 3, sizeLik = "t", censorRunaways = FALSE)
  expect_identical(cens$runaways, off$runaways)
  expect_identical(cens$allFireSizes, off$allFireSizes)       # burned sizes, bounded by the 100-pixel buffers
  expect_identical(cens$annualAreaByRep, off$annualAreaByRep)
  expect_true(all(cens$allFireSizes <= 100L))
  expect_false(isTRUE(all.equal(cens$SNLL_FS, off$SNLL_FS)))
  ## the default rule stops a runaway early
  seen$stop <- NULL
  runawayObjective(sp, seed = 3, sizeLik = "t", censorRunaways = TRUE)
  expect_true(any(seen$stop))
})

test_that(".objfunSpreadFit() hands censoredRunaways on, drops the runaway size substitutes with it, and refuses what it cannot do", {
  seen <- list()
  local_mocked_bindings(
    objFunInner = function(runawaySize, runawayBufferMultiple, censorRunaways, censoredRunaways, ...) {
      seen[[length(seen) + 1L]] <<- list(size = runawaySize, multiple = runawayBufferMultiple,
                                         censor = censorRunaways, censored = censoredRunaways)
      list(SNLL_FS = 0)
    },
    .package = "fireSenseUtils"
  )
  land <- terra::rast(nrows = 20, ncols = 20, xmin = 0, xmax = 20, ymin = 0, ymax = 20, vals = 1)
  run <- function(...) {
    seen <<- list()
    fireSenseUtils::.objfunSpreadFit(
      par = c(0.27, 1, 1, 1), landscape = land,
      annualDTx1000 = list(year2004 = data.table::data.table(pixelID = 1:400, cov = 0L)),
      nonAnnualDTx1000 = list(`year2004` = data.table::data.table(pixelID = 1:400)),
      formulaToFit = "~ 0 + cov",
      historicalFires = list(year2004 = data.frame(size = 40L, cells = 1L, ids = 1L)),
      fireBufferedListDT = list(year2004 = data.table::data.table(ids = 1L, pixelID = 1:400)),
      doAssertions = FALSE, sizeLik = "t", ...
    )
    seen[[1]]
  }
  expect_equal(run(runawayBufferMultiple = 2), list(size = NULL, multiple = 2, censor = TRUE, censored = FALSE))
  expect_equal(run(runawayBufferMultiple = 2, censoredRunaways = TRUE),
               list(size = NULL, multiple = NULL, censor = TRUE, censored = TRUE))
  expect_error(run(censoredRunaways = TRUE, penaliseRunaways = FALSE), "penaliseRunaways")
  expect_error(run(censoredRunaways = TRUE, sizeLik = "kde"), "sizeLik")
})
