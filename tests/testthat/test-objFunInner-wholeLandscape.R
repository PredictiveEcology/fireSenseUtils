## A replicate fire cut at an open buffer edge is scored at the size it reached inside its buffer. With
## `wholeLandscape = TRUE` every fire spreads on all the pixels of the tables (the whole landscape) instead, with
## spreadProb computed by spreadCpp() only on the cells a fire reaches (link = "powerLaw"). See
## helper-wholeLandscape.R.

## the buffer of wlFixture() with its edge columns saying no pixel is open: the buffered objective on the whole
## landscape's tables then spreads every fire on the whole landscape with a full spreadProb vector (the reference)
closedBuffer <- function(side = 11L) {
  fx <- wlFixture(side)
  data.table::set(fx$buf, NULL, c("edge", "open"), list(FALSE, FALSE))
  fx
}

test_that("on the whole landscape the fires burn the cells of the full spreadProb vector, to the bit", {
  skip_if_not_installed("SpaDES.tools")
  for (args in list(list(), list(yearSpreadSD = 0.5), list(yearSpreadSD = 0.5, escapeMinPx = 9, jumpTries = 20,
                                                          jumpMeanDist = 3))) {
    full <- do.call(wlObjective, c(list(fx = closedBuffer(), annPix = seq_len(3600L), returnSims = TRUE,
                                        returnBurned = TRUE, seed = 3), args))
    lazy <- do.call(wlObjective, c(list(whole = TRUE, returnSims = TRUE, returnBurned = TRUE, seed = 3), args))
    expect_identical(lazy$burned, full$burned)
    expect_identical(lazy$sims, full$sims)
    expect_gt(max(lazy$sims$sim), 121)                       # past what the 11 x 11 buffer holds
    ## and every term of the objective
    full <- do.call(wlObjective, c(list(fx = closedBuffer(), annPix = seq_len(3600L), seed = 4), args))
    lazy <- do.call(wlObjective, c(list(whole = TRUE, seed = 4), args))
    expect_identical(lazy, full)
  }
})

test_that("buffers over the whole landscape: wholeLandscape changes nothing", {
  skip_if_not_installed("SpaDES.tools")
  for (args in list(list(), list(yearSpreadSD = 0.5, escapeMinPx = 9, jumpTries = 20, jumpMeanDist = 3))) {
    off <- do.call(wlObjective, c(list(side = 60L, seed = 2), args))
    on <- do.call(wlObjective, c(list(side = 60L, whole = TRUE, seed = 2), args))
    expect_true(on$runaways == 0)
    expect_equal(off$runaways, 0)
    expect_identical(on, off)
  }
})

test_that("on the whole landscape no fire is cut, and fires burn past the buffer", {
  skip_if_not_installed("SpaDES.tools")
  off <- wlObjective(yearSpreadSD = 0.5, Nreps = 40L)
  on <- wlObjective(yearSpreadSD = 0.5, Nreps = 40L, whole = TRUE)
  expect_gt(off$runaways, 0)
  expect_equal(on$runaways, 0)
  expect_identical(on$nSims, off$nSims)
  expect_lte(max(off$allFireSizes), 121)
  expect_gt(max(on$allFireSizes), 121)
})

test_that("with wholeLandscape FALSE the objective is what it was before it existed", {
  skip_if_not_installed("SpaDES.tools")
  ## equal to the last bit (tolerance = 0); identical() would also compare data.table internals
  ## fixtures/objFunInnerBuffered.rds was made by the same calls on fireSenseUtils 0.2.3.9155 (c86d0a4), on Linux
  ## x86_64. Other platforms' maths libraries differ in the last bit (macOS arm64: 1.8e-16 in crpsFire), so there the
  ## comparison allows that.
  tol <- if (identical(Sys.info()[["sysname"]], "Linux") && identical(R.version$arch, "x86_64")) 0 else 1e-12
  before <- readRDS(test_path("fixtures", "objFunInnerBuffered.rds"))
  expect_equal(wlObjective(p = 0.2, yearSpreadSD = 0.4, Nreps = 30L, seed = 5), before$crps, tolerance = tol)
  expect_equal(wlObjective(p = 0.2, Nreps = 30L, seed = 6), before$noEffect, tolerance = tol)
  sims <- wlObjective(p = 0.2, yearSpreadSD = 0.4, Nreps = 30L, seed = 5, returnSims = TRUE, returnBurned = TRUE)
  expect_equal(sims, before$sims, tolerance = tol)
})

test_that("wholeLandscape needs the powerLaw link, and not penaliseRunaways", {
  expect_error(checkWholeLandscape(TRUE, TRUE, "powerLaw"), "penaliseRunaways = FALSE")
  expect_error(checkWholeLandscape(TRUE, FALSE, NULL), "powerLaw")
  expect_false(checkWholeLandscape(FALSE, TRUE, NULL))
  expect_false(checkWholeLandscape(NULL, TRUE, NULL))
  expect_true(checkWholeLandscape(TRUE, FALSE, "powerLaw"))
})

test_that("fitted ignitions outside every buffer of their year are warned about", {
  fx <- wlFixture()
  fires <- data.table::data.table(cells = c(wlCell(30L, 30L), wlCell(2L, 2L)), size = 20L, ids = 1:2)
  expect_warning(n <- warnIgnitionsOutsideBuffers(list(y1 = fires, y2 = fx$fires), list(y1 = fx$buf, y2 = fx$buf)),
                 "1 of 3 fitted ignitions are outside every buffer of their year \\(y1: 1\\)")
  expect_identical(unname(n), c(1L, 0L))
  expect_no_warning(warnIgnitionsOutsideBuffers(list(y2 = fx$fires), list(y2 = fx$buf)))
})

test_that("under gateRefusal = 'simulate' the gates do not run, and the objective is the same", {
  skip_if_not_installed("SpaDES.tools")
  calls <- 0L
  realGate <- spreadProbGateTest
  local_mocked_bindings(spreadProbGateTest = function(...) { calls <<- calls + 1L; realGate(...) })
  for (whole in c(FALSE, TRUE)) {
    calls <- 0L
    skipped <- wlObjective(whole = whole, yearSpreadSD = 0.5, seed = 7)
    expect_identical(calls, 0L)
    ## not for verbose, returnSims or plotting either
    utils::capture.output(wlObjective(whole = whole, seed = 7, verbose = 2))
    grDevices::pdf(nullFile <- tempfile(fileext = ".pdf"))
    ## the plots run without the gates' panel (the toy buffers are small, which the plotting path warns about)
    expect_no_error(suppressWarnings(wlObjective(whole = whole, seed = 7, plot.it = "screen")))
    grDevices::dev.off(); unlink(nullFile)
    sims <- wlObjective(whole = whole, seed = 7, returnSims = TRUE)
    expect_identical(calls, 0L)
    expect_false(is.null(sims$pSummary))
    ## "minLik" decides with them; gates that pass leave the objective as it is without them
    ran <- wlObjective(whole = whole, yearSpreadSD = 0.5, seed = 7, gateRefusal = "minLik",
                       lanscape1stQuantileThresh = 1, maxFireSpread = 1)
    expect_identical(calls, 1L)
    expect_identical(ran, skipped)
  }
})

test_that("a fit evaluation (returnSims FALSE) under 'simulate' does no gate or pGate work", {
  skip_if_not_installed("SpaDES.tools")
  n <- new.env(); n$gate <- n$summ <- n$p <- 0L
  realGate <- spreadProbGateTest; realSumm <- spreadProbSummary; realP <- spreadProbFromCovs
  local_mocked_bindings(spreadProbGateTest = function(...) { n$gate <- n$gate + 1L; realGate(...) },
                        spreadProbSummary = function(...) { n$summ <- n$summ + 1L; realSumm(...) },
                        spreadProbFromCovs = function(...) { n$p <- n$p + 1L; realP(...) })
  wlObjective(whole = TRUE, yearSpreadSD = 0.5, seed = 8)
  expect_identical(c(n$gate, n$summ, n$p), c(0L, 0L, 0L))  # on the whole landscape not even the spreadProb column
  wlObjective(yearSpreadSD = 0.5, seed = 8)
  expect_identical(c(n$gate, n$summ), c(0L, 0L))           # buffered: the column feeds the spread, nothing else
  wlObjective(whole = TRUE, seed = 8, returnSims = TRUE)
  expect_identical(c(n$gate, n$summ), c(0L, 1L))           # the saturation summary, only when asked for
})
