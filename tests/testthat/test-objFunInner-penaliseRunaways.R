## A simulated fire that burns any pixel of the edge ring of its OWN buffer has run away: it is
## censored ("at least this big"). Before, a fire was capped at multiplier() of its observed size and
## "capHit" meant reaching that cap; the cap is gone, and spread is bounded only by the year's buffers.
## `objFunInner()` runs the real spreadCpp() here; only the spreadProb surface is faked.

r60 <- terra::rast(nrows = 60, ncols = 60, xmin = 0, xmax = 60, ymin = 0, ymax = 60)
cellOf60 <- function(row, col) (row - 1L) * 60L + col

## five fires, each in its own n x n buffer block, ignited at the block's centre
runawayFixture <- function(n = 10L) {
  tops <- c(2L, 2L, 22L, 22L, 42L); lefts <- c(2L, 22L, 2L, 22L, 2L)
  buf <- data.table::rbindlist(lapply(1:5, function(i)
    data.table::data.table(ids = i, pixelID = as.vector(outer(tops[i] + 0:(n - 1L), lefts[i] + 0:(n - 1L),
                                                              cellOf60)), buffer = 0L)))
  centre <- cellOf60(tops + n %/% 2L, lefts + n %/% 2L)
  list(buf = buf, fires = data.table::data.table(cells = centre, size = c(25L, 30L, 20L, 28L, 22L), ids = 1:5))
}

## spreadProb is `sp` inside the buffers, 0 elsewhere: the mocked surface is over the buffers' pixels only
runawayObjective <- function(sp, seed = 1, fx = runawayFixture(), Nreps = 6L, ...) {
  pix <- fx$buf$pixelID
  local_mocked_bindings(
    paramsSeparate = function(...) list(logisticPars = c(0.27, 1, 1), covPars = 1),
    spreadProbFromIntegerCovs = function(...) data.table::data.table(pixelID = pix, cov = 0),
    logisticAll = function(...) seq(sp[1], sp[length(sp)], length.out = length(pix)), # a flat surface is refused as not spread out
    .package = "fireSenseUtils"
  )
  set.seed(seed)
  fireSenseUtils:::objFunInner(
    yr = "2004", annDTx1000 = NULL, par = 1, parsModel = NULL,
    annualFires = fx$fires, nonAnnualDTx1000 = NULL, annualFireBufferedDT = data.table::copy(fx$buf),
    indexNonAnnual = NULL, colsToUse = "cov", covMinMax = NULL,
    mutuallyExclusive = NULL, doAssertions = FALSE, maxFireSpread = 1,
    lowerSpreadProb = 0.05, cells = numeric(3600), lanscape1stQuantileThresh = 1,
    weighted = FALSE, r = r60, Nreps = Nreps,
    doSNLL_FSTest = TRUE, doMADTest = TRUE, doADTest = TRUE, doYearArea = TRUE,
    plot.it = FALSE, verbose = 0, ...
  )
}

test_that("a saturated parameter set runs away in every replicate and scores far worse with penaliseRunaways", {
  skip_if_not_installed("SpaDES.tools")
  ## buffers of 36 pixels around fires of 20-30: a fire that fills its buffer looks about right to the
  ## likelihood unless it is recognised as a runaway
  fx <- runawayFixture(6L)
  off <- runawayObjective(c(0.9, 0.99), fx = fx)
  on <- runawayObjective(c(0.9, 0.99), fx = fx, runawaySize = 3600)
  expect_equal(on$runaways, on$nSims)             # every replicate reached its buffer's edge
  expect_equal(off$runaways, off$nSims)           # counted whether or not they are penalised
  expect_gt(on$SNLL_FS, off$SNLL_FS + 50)
  expect_true(all(on$allFireSizes == 3600))       # censored sizes reach the AD term
  expect_true(all(off$allFireSizes <= 36))        # uncensored: what the fire burned, bounded by its buffer
  expect_gt(min(on$annualAreaByRep), max(off$annualAreaByRep))
})

test_that("by default a runaway is censored in the per-fire likelihood only; the area terms see what it burned", {
  skip_if_not_installed("SpaDES.tools")
  fx <- runawayFixture(6L)
  off <- runawayObjective(c(0.9, 0.99), fx = fx)
  on <- runawayObjective(c(0.9, 0.99), fx = fx, censorRunaways = TRUE)  # runawaySize NULL, as .objfunSpreadFit() now passes
  expect_equal(on$runaways, on$nSims)
  expect_gt(on$SNLL_FS, off$SNLL_FS + 50)          # the likelihood still penalises them
  expect_true(all(on$allFireSizes <= 36))          # not 3600, nor any landscape-wide size
  expect_lte(max(on$annualAreaByRep), 5 * 36)      # five fires, each bounded by its 36-pixel buffer
})

test_that("with no fire reaching the edge of its buffer, the objective is identical with and without the penalty", {
  skip_if_not_installed("SpaDES.tools")
  off <- runawayObjective(c(0.06, 0.08))
  on <- runawayObjective(c(0.06, 0.08), runawaySize = 3600)
  expect_equal(on$runaways, 0)
  expect_identical(on[setdiff(names(on), c("runaways", "nSims"))],
                   off[setdiff(names(off), c("runaways", "nSims"))])
})

test_that("the runaway count rises with spread, and is neither all nor none in between", {
  skip_if_not_installed("SpaDES.tools")
  ## 10 x 10 buffers, fires start at the centre, 5 cells from the edge
  mid <- vapply(c(0.1, 0.14, 0.2, 0.4), function(sp) runawayObjective(sp * c(0.9, 1.1), seed = 3)$runaways,
                numeric(1))
  expect_true(all(diff(mid) >= 0))
  expect_true(mid[1] < 30 && mid[4] > 0)
})

test_that("no size cap: spreadCpp gets no maxSize, and a fire can exceed multiplier() of its observed size", {
  skip_if_not_installed("SpaDES.tools")
  seen <- new.env()
  orig <- SpaDES.tools::spreadCpp
  local_mocked_bindings(
    spreadCpp = function(...) {
      seen$args <- names(list(...))
      orig(...)
    },
    .package = "SpaDES.tools")
  ## 18 x 18 = 324 pixel buffers; the old cap for these fires (20-30 pixels) was 220-318
  out <- runawayObjective(c(0.9, 0.99), returnSims = TRUE, fx = runawayFixture(18L))
  expect_false("maxSize" %in% seen$args)
  expect_true(all(out$sims$sim <= 324L) && stats::median(out$sims$sim) >= 300)
  expect_true(stats::median(out$sims$sim) > max(multiplier(out$sims$size, minSize = 100)))
})

test_that("returnSims keeps the simulated sizes, penalty or not, with the same columns as before", {
  skip_if_not_installed("SpaDES.tools")
  off <- runawayObjective(c(0.9, 0.99), returnSims = TRUE)
  on <- runawayObjective(c(0.9, 0.99), returnSims = TRUE, runawaySize = 3600)
  expect_identical(on$sims, off$sims)
  expect_named(on$sims, c("yr", "rep", "initialLocus", "sim", "ids", "size"))
  expect_equal(nrow(on$sims), 6L * 5L)
  expect_true(all(on$sims$sim <= 100L) && stats::median(on$sims$sim) >= 95)  # real sizes, not 3600: the 10 x 10 buffer, filled
  withBurned <- runawayObjective(c(0.9, 0.99), returnSims = TRUE, returnBurned = TRUE)
  expect_named(withBurned$burned, c("rep", "pixelID"))
})

test_that(".objfunSpreadFit() censors runaways by default without a runawaySize, and switches off", {
  seen <- list()
  local_mocked_bindings(
    objFunInner = function(runawaySize, censorRunaways, ...) {
      seen[[length(seen) + 1L]] <<- list(size = runawaySize, censor = censorRunaways)
      list(SNLL_FS = 0)
    },
    .package = "fireSenseUtils"
  )
  land <- terra::rast(nrows = 20, ncols = 20, xmin = 0, xmax = 20, ymin = 0, ymax = 20, vals = 1)
  land[1:40] <- NA
  run <- function(...) {
    seen <<- list()
    fireSenseUtils::.objfunSpreadFit(
      par = c(0.27, 1, 1, 1), landscape = land,
      annualDTx1000 = list(year2004 = data.table::data.table(pixelID = 1:400, cov = 0L)),
      nonAnnualDTx1000 = list(`year2004` = data.table::data.table(pixelID = 1:400)),
      formulaToFit = "~ 0 + cov",
      historicalFires = list(year2004 = data.frame(size = 40L, cells = 1L, ids = 1L)),
      fireBufferedListDT = list(year2004 = data.table::data.table(ids = 1L, pixelID = 1:400)),
      indexNonAnnual = data.table::data.table(ind = 1L, date = "2004"), doAssertions = FALSE,
      lrgSmallFireYears = list(1L, integer(0)), ...
    )
    seen
  }
  ## the landscape's 360 non-NA pixels are no longer a runaway's size: one runaway replicate made a whole
  ## year look like the landscape burned, and fits chose parameters that under-burned (2026-10-01)
  expect_equal(run()[[1]], list(size = NULL, censor = TRUE))
  expect_equal(run(runawaySize = 99)[[1]], list(size = 99, censor = TRUE))
  expect_equal(run(penaliseRunaways = FALSE)[[1]], list(size = NULL, censor = FALSE))
})

test_that("capSizes and penaliseCapHits are deprecated: one warning, then ignored", {
  local_mocked_bindings(objFunInner = function(...) list(SNLL_FS = 0), .package = "fireSenseUtils")
  rm(list = ls(fireSenseUtils:::.capArgsWarned), envir = fireSenseUtils:::.capArgsWarned)
  land <- terra::rast(nrows = 20, ncols = 20, xmin = 0, xmax = 20, ymin = 0, ymax = 20, vals = 1)
  call1 <- function(...) fireSenseUtils::.objfunSpreadFit(
    par = c(0.27, 1, 1, 1), landscape = land,
    annualDTx1000 = list(year2004 = data.table::data.table(pixelID = 1:400, cov = 0L)),
    nonAnnualDTx1000 = list(`year2004` = data.table::data.table(pixelID = 1:400)),
    formulaToFit = "~ 0 + cov",
    historicalFires = list(year2004 = data.frame(size = 40L, cells = 1L, ids = 1L)),
    fireBufferedListDT = list(year2004 = data.table::data.table(ids = 1L, pixelID = 1:400)),
    indexNonAnnual = data.table::data.table(ind = 1L, date = "2004"), doAssertions = FALSE, ...)
  expect_warning(call1(capSizes = FALSE, penaliseCapHits = TRUE), "penaliseRunaways")
  expect_no_warning(call1(capSizes = FALSE))
})

test_that("returnTerms reports the first block's SNLL per year and whether it bailed", {
  ## the early stop's own numbers, so a caller calibrating `thresh` need not read the printed log
  mk <- function(snll, gated = FALSE, thresh = Inf, ...) {
    local_mocked_bindings(
      objFunInner = function(yr, ...) list(SNLL_FS = snll[[yr]], runaways = 0L, nSims = 0L, gated = gated[[yr]]),
      .package = "fireSenseUtils")
    land <- terra::rast(nrows = 20, ncols = 20, xmin = 0, xmax = 20, ymin = 0, ymax = 20, vals = 1)
    yrs <- c("2001", "2002", "2003")
    dt <- function(...) data.table::data.table(...)
    fireSenseUtils::.objfunSpreadFit(
      par = c(0.27, 1, 1, 1), landscape = land,
      annualDTx1000 = stats::setNames(lapply(yrs, function(y) dt(pixelID = 1:2, cov = 100L)), yrs),
      nonAnnualDTx1000 = list(`2001_2003` = dt(pixelID = 1:2)),
      formulaToFit = "~ 0 + cov",
      historicalFires = list(`2001` = data.frame(size = c(100, 200), cells = 1:2, ids = 1:2),
                             `2002` = data.frame(size = c(150, 250), cells = 3:4, ids = 1:2),
                             `2003` = data.frame(size = c(5, 6), cells = 5:6, ids = 1:2)),
      fireBufferedListDT = stats::setNames(lapply(yrs, function(y) dt(pixelID = 1:2, ids = 1L)), yrs),
      tests = "snll_fs", Nreps = 1L, doAssertions = FALSE, verbose = 0, thresh = thresh,
      returnTerms = TRUE, ...)
  }
  snll <- c("2001" = 100, "2002" = 300, "2003" = 40)  # 2001 and 2002 are the first block (the largest areas)
  ok <- mk(snll, gated = c("2001" = FALSE, "2002" = FALSE, "2003" = FALSE))
  expect_equal(ok[["firstBlockSNLL"]], (100 + 300) / 2)
  expect_equal(ok[["bailed"]], 0)
  gated <- mk(snll, gated = c("2001" = FALSE, "2002" = TRUE, "2003" = FALSE))
  expect_equal(gated[["bailed"]], 1)               # a year refused as too burny / not spread out
  expect_equal(gated[["firstBlockSNLL"]], 200)     # still reported, with thresh = Inf
  over <- mk(snll, gated = c("2001" = FALSE, "2002" = FALSE, "2003" = FALSE), thresh = 150)
  expect_equal(over[["bailed"]], 1)                # 200 per year against a threshold of 150
  expect_equal(over[["firstBlockSNLL"]], 200)
})

## one fire whose buffer is pixels 1..(nRing + 10); the first nRing of them are its edge ring. spreadCpp is
## mocked to burn the ignition plus `touch` ring pixels, so the rule is tested, not the spread
runawayEdgeRuns <- function(nRing, touch, ...) {
  pix <- seq_len(nRing + 10L)
  fx <- list(buf = data.table::data.table(ids = 1L, pixelID = pix, buffer = 0L, edge = pix <= nRing),
             fires = data.table::data.table(cells = nRing + 5L, size = 20L, ids = 1L))
  local_mocked_bindings(
    cropToCells = function(r, cells, ...) list(r = r, ncell = terra::ncell(r), toCrop = identity, toFull = identity),
    .package = "fireSenseUtils")
  local_mocked_bindings(
    spreadCpp = function(landscape, loci, ...)
      data.table::data.table(initialLocus = loci, indices = c(loci, seq_len(touch))),
    .package = "SpaDES.tools")
  runawayObjective(c(0.06, 0.08), fx = fx, ...)$runaways
}

test_that("a replicate is a runaway only when it burns at least max(3, 1% of ring) distinct ring pixels", {
  skip_if_not_installed("SpaDES.tools")
  expect_equal(runawayEdgeRuns(100L, 0L), 0)
  expect_equal(runawayEdgeRuns(100L, 1L), 0)      # one touched pixel is luck
  expect_equal(runawayEdgeRuns(100L, 2L), 0)
  expect_equal(runawayEdgeRuns(100L, 3L), 6)      # all 6 replicates
  expect_equal(runawayEdgeRuns(1000L, 9L), 0)     # ceiling(0.01 * 1000) = 10
  expect_equal(runawayEdgeRuns(1000L, 10L), 6)
  expect_equal(runawayEdgeRuns(2L, 1L), 0)        # ring smaller than runawayEdgeMin: k is the ring size
  expect_equal(runawayEdgeRuns(2L, 2L), 6)
  expect_equal(runawayEdgeRuns(100L, 1L, runawayEdgeMin = 1, runawayEdgeFrac = 0), 6)
})

## The early stop: with the penalty, spreadCpp() is told to stop a fire once it has burned k of its edge-ring
## pixels. `stopRuns(stop)` runs the objective with spreadCpp() wrapped to keep or drop those arguments, and
## records how many cells each call burned.
stopRuns <- function(sp, stop, fx, seed = 1, Nreps = 6L, ...) {
  orig <- SpaDES.tools::spreadCpp
  rec <- new.env(); rec$burned <- integer(0); rec$stopArgs <- NULL
  local_mocked_bindings(
    spreadCpp = function(...) {
      a <- list(...)
      rec$stopArgs <- names(Filter(Negate(is.null), a[c("stopCells", "stopEvent", "stopAt")]))
      if (!stop) a[c("stopCells", "stopEvent", "stopAt")] <- NULL
      out <- do.call(orig, a)
      rec$burned <- c(rec$burned, NROW(out))
      out
    },
    .package = "SpaDES.tools")
  res <- runawayObjective(sp, seed = seed, fx = fx, Nreps = Nreps, ...)
  list(res = res, burned = rec$burned, stopArgs = rec$stopArgs)
}

oneFire <- function(n = 10L) {
  fx <- runawayFixture(n)
  list(buf = fx$buf[fx$buf$ids == 1L, , drop = FALSE], fires = fx$fires[1L, , drop = FALSE])
}

test_that("early stop: one fire, one replicate, the score and the classification are those of the full burn", {
  skip_if_not_installed("SpaDES.tools")
  skip_if_not("stopCells" %in% names(formals(SpaDES.tools::spreadCpp)))
  fx <- oneFire()
  sawRunaway <- sawNone <- FALSE
  for (sp in list(c(0.06, 0.08), c(0.12, 0.16), c(0.2, 0.3), c(0.9, 0.99)))
    for (seed in 1:5) {
      on <- stopRuns(sp, TRUE, fx, seed = seed, Nreps = 1L, runawaySize = 3600)
      off <- stopRuns(sp, FALSE, fx, seed = seed, Nreps = 1L, runawaySize = 3600)
      expect_setequal(on$stopArgs, c("stopCells", "stopEvent", "stopAt"))
      expect_equal(on$res$runaways, off$res$runaways)
      expect_identical(on$res, off$res)
      sawRunaway <- sawRunaway || on$res$runaways > 0
      sawNone <- sawNone || on$res$runaways == 0
    }
  expect_true(sawRunaway && sawNone)       # the comparison covers both kinds of replicate
  ## returnTerms-style fields are in the compared result: SNLL_FS, the AD and area terms
  expect_true(all(c("SNLL_FS", "allFireSizes", "annualAreaByRep") %in% names(on$res)))
})

test_that("early stop: a runaway burns fewer cells, and only with the penalty and without returnSims", {
  skip_if_not_installed("SpaDES.tools")
  skip_if_not("stopCells" %in% names(formals(SpaDES.tools::spreadCpp)))
  fx <- oneFire(18L)
  on <- stopRuns(c(0.9, 0.99), TRUE, fx, runawaySize = 3600)
  off <- stopRuns(c(0.9, 0.99), FALSE, fx, runawaySize = 3600)
  expect_equal(on$res$runaways, on$res$nSims)
  expect_lt(sum(on$burned), sum(off$burned))
  expect_equal(off$res$runaways, off$res$nSims)
  ## no penalty: a runaway is not censored, so it keeps burning and spreadCpp() is not asked to stop
  noPen <- stopRuns(c(0.9, 0.99), TRUE, fx, runawaySize = NULL)
  expect_length(noPen$stopArgs, 0L)
  expect_equal(sum(noPen$burned), sum(off$burned))
  sims <- stopRuns(c(0.9, 0.99), TRUE, fx, runawaySize = 3600, returnSims = TRUE)
  expect_length(sims$stopArgs, 0L)
})

test_that("early stop: timing on the five-fire fixture (reported, not asserted)", {
  skip_if_not_installed("SpaDES.tools")
  skip_if_not("stopCells" %in% names(formals(SpaDES.tools::spreadCpp)))
  fx <- runawayFixture(18L)
  tm <- function(stop) system.time(for (i in 1:20) stopRuns(c(0.9, 0.99), stop, fx, runawaySize = 3600))[["elapsed"]]
  tOn <- tm(TRUE); tOff <- tm(FALSE)
  message("early stop timing, 20 evaluations of 5 saturated fires x 6 replicates: with stop ",
          round(tOn, 2), " s, without ", round(tOff, 2), " s")
  expect_true(is.finite(tOn) && is.finite(tOff))
})
