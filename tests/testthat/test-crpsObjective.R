## tests = "crps": per-fire CRPS on log size, adTest, areaDist, CRPS of the log total area and of the 10 largest
## fires, all proper scores. Its inner pieces on the toy landscape are tested in test-bufferIndependentObjective.R.

test_that("crpsEnsemble() is E|X - y| - 0.5 E|X - X'|, worked by hand", {
  expect_equal(crpsEnsemble(c(0, 1), 0), 0.5 - 0.5 * 0.5)      # E|X-y| = 1/2, E|X-X'| = (0+1+1+0)/4
  expect_equal(crpsEnsemble(c(1, 2, 3), 2), 2/3 - 0.5 * 8/9)   # pairs: 0,1,2,1,0,1,2,1,0
  expect_equal(crpsEnsemble(rep(5, 7), 5), 0)                    # a point mass at the observation
  expect_equal(crpsEnsemble(3, 1), 2)                            # one member: the absolute error
  x <- c(2.5, 0.1, 7, 3, 3, 1)                                   # unsorted, ties
  allPairs <- mean(abs(outer(x, x, "-")))
  expect_equal(crpsEnsemble(x, 2), mean(abs(x - 2)) - 0.5 * allPairs)
  skip_if_not_installed("scoringRules")
  expect_equal(crpsEnsemble(x, 2), scoringRules::crps_sample(2, x))
})

test_that("crpsFireSizes() works on log size, and a cut replicate is at least the observed size", {
  expect_equal(crpsFireSizes(10, c(5, 20)), crpsEnsemble(log(c(5, 20)), log(10)))
  ## cut replicates: 5 stopped short of the observed 10 counts as 10; 20 stays 20
  expect_equal(crpsFireSizes(10, c(5, 20), cut = c(TRUE, TRUE)), crpsEnsemble(log(c(10, 20)), log(10)))
  expect_equal(crpsFireSizes(10, c(5, 20), cut = c(FALSE, TRUE)), crpsFireSizes(10, c(5, 20)))
  expect_equal(crpsFireSizes(10, numeric(0)), log(10))           # no replicate: 1 pixel
})

test_that("crpsTotalArea() scores the log of the total over all years; a replicate with no fire is 1 pixel", {
  rs <- data.table::data.table(rep = c(1, 1, 2, 2, 2), N = c(10, 20, 5, 5, 20))   # totals 30 and 30; rep 3 has none
  expect_equal(crpsTotalArea(rs, c(20, 10), Nreps = 3L), crpsEnsemble(log(c(30, 30, 1)), log(30)))
})

test_that("crpsTopFires() sums the per-rank CRPS of the sorted largest fires, and copes with fewer than 10", {
  rs <- data.table::data.table(rep = c(1, 1, 1, 2, 2), N = c(50, 100, 3, 90, 40))
  obs <- c(60, 100)                       # 2 fitted fires: ranks 1 and 2 only
  expect_equal(crpsTopFires(rs, obs, Nreps = 2L),
               crpsEnsemble(log(c(100, 90)), log(100)) + crpsEnsemble(log(c(50, 40)), log(60)))
  ## a replicate with fewer fires than ranks has 1 pixel there
  expect_equal(crpsTopFires(rs[rep == 1], c(100, 60, 30), Nreps = 2L),
               crpsEnsemble(log(c(100, 1)), log(100)) + crpsEnsemble(log(c(50, 1)), log(60)) +
                 crpsEnsemble(log(c(3, 1)), log(30)))
  ## more than 10 observed fires: the 10 largest
  big <- data.table::data.table(rep = 1L, N = 1:12)
  expect_equal(crpsTopFires(big, 12:1, Nreps = 1L), 0)
  expect_equal(crpsTopFires(big, c(12:1, 100), Nreps = 1L),
               crpsEnsemble(log(12), log(100)) + sum(vapply(2:10, function(i) abs(log(13 - i) - log(14 - i)), 1)))
})

test_that("the default weights are the 2026-10-28 profile's, and crpsWeights overrides by name", {
  expect_identical(crpsWeightsDefault(), list(fire = 1, adTest = 0.34, areaDist = 2.6, totalArea = 164, top10 = 16.5))
})

## .objfunSpreadFit() with the per-year work mocked: one year, 2 fitted fires (25, 60 pixels), Nreps = 2
crpsFit <- function(..., sizes = data.table::data.table(rep = c(1, 1, 2, 2), N = c(20, 70, 30, 50))) {
  local_mocked_bindings(
    objFunInner = function(yr, ...) list(SNLL_FS = 5, runaways = 0L, nSims = 4L, gated = FALSE,
                                         allFireSizes = sizes$N, crpsFire = 1.5, crpsSizes = sizes),
    .package = "fireSenseUtils")
  land <- terra::rast(nrows = 20, ncols = 20, xmin = 0, xmax = 20, ymin = 0, ymax = 20, vals = 1)
  fireSenseUtils::.objfunSpreadFit(
    par = c(0.27, 1, 1, 1), landscape = land,
    annualDTx1000 = list(`2004` = data.table::data.table(pixelID = 1:400, cov = 0L)),
    nonAnnualDTx1000 = list(`2004` = data.table::data.table(pixelID = 1:400)), formulaToFit = "~ 0 + cov",
    historicalFires = list(`2004` = data.frame(size = c(25, 60), cells = 1:2, ids = 1:2)),
    fireBufferedListDT = list(`2004` = data.table::data.table(ids = 1L, pixelID = 1:400)),
    Nreps = 2L, doAssertions = FALSE, verbose = 0, returnTerms = TRUE, ...)
}

test_that("tests = \"crps\" is the weighted sum of its five terms, with no annual-area term and no early stop", {
  w <- crpsWeightsDefault()
  fit <- crpsFit(tests = "crps", yearAreaWeight = "auto", thresh = 0.001)
  obs <- c(25, 60); sizes <- c(20, 70, 30, 50)
  ta <- crpsEnsemble(log(c(90, 80)), log(85))
  top <- crpsEnsemble(log(c(70, 50)), log(60)) + crpsEnsemble(log(c(20, 30)), log(25))
  ad <- fireSenseUtils:::adStatistic(sizes, obs)
  aw <- fireSenseUtils:::areaWeightedCvM(sizes, obs)
  expect_equal(fit[["crpsFire"]], 1.5)
  expect_equal(fit[["crpsTotalArea"]], ta)
  expect_equal(fit[["crpsTop10"]], top)
  expect_equal(fit[["adTest"]], w$adTest * ad)
  expect_equal(fit[["areaDist"]], w$areaDist * aw)
  expect_true(is.na(fit[["yearArea"]]))
  expect_equal(fit[["objective"]], w$fire * 1.5 + w$adTest * ad + w$areaDist * aw + w$totalArea * ta + w$top10 * top)
  ## weights by name; the rest stay default
  fit2 <- crpsFit(tests = "crps", crpsWeights = list(totalArea = 0, top10 = 1))
  expect_equal(fit2[["objective"]], 1.5 + w$adTest * ad + w$areaDist * aw + top)
})

test_that("a refused year scores as every fire left at one pixel (gateRefusal = \"minLik\")", {
  ## objFunInner's refused-year output, through the whole objective
  sizes <- data.table::data.table(rep = rep(1:2, each = 2), N = 1L)
  fit <- crpsFit(tests = "crps", sizes = sizes)
  expect_equal(fit[["crpsTotalArea"]], crpsEnsemble(log(c(2, 2)), log(85)))
  expect_equal(fit[["crpsTop10"]], crpsEnsemble(log(c(1, 1)), log(60)) + crpsEnsemble(log(c(1, 1)), log(25)))
})

test_that("the defaults leave the objective's value as it was before the CRPS options", {
  ## A real objFunInner() and spreadCpp() on a 20 x 20 landscape, 3 years, 6 fires, seed 42. The values were
  ## computed on development before this change (2026-10-28).
  skip_if_not_installed("SpaDES.tools")
  dt <- function(...) data.table::data.table(...)
  yrs <- c("2001", "2002", "2003")
  land <- terra::rast(nrows = 20, ncols = 20, xmin = 0, xmax = 20, ymin = 0, ymax = 20, vals = 1)
  local_mocked_bindings(
    paramsSeparate = function(...) list(logisticPars = c(0.27, 1, 1), covPars = 1),
    spreadProbFromIntegerCovs = function(...) dt(pixelID = 1:400, cov = 0),
    logisticAll = function(...) seq(0.18, 0.27, length.out = 400),
    .package = "fireSenseUtils")
  ob <- function(...) {
    set.seed(42)
    fireSenseUtils::.objfunSpreadFit(
      par = c(0.27, 1, 1, 1), landscape = land,
      annualDTx1000 = stats::setNames(lapply(yrs, function(y) dt(pixelID = 1:400, cov = 0L)), yrs),
      nonAnnualDTx1000 = list(`2001_2003` = dt(pixelID = 1:400)), formulaToFit = "~ 0 + cov",
      historicalFires = list(`2001` = data.frame(size = c(30, 60), cells = c(45L, 250L), ids = 1:2),
                             `2002` = data.frame(size = c(12, 90), cells = c(100L, 300L), ids = 1:2),
                             `2003` = data.frame(size = c(8, 25), cells = c(150L, 360L), ids = 1:2)),
      fireBufferedListDT = stats::setNames(lapply(yrs, function(y) dt(pixelID = 1:400, ids = rep(1:2, each = 200))), yrs),
      Nreps = 8L, doAssertions = FALSE, verbose = 0, thresh = Inf, returnTerms = TRUE, ...)
  }
  a <- ob(tests = c("snll_fs", "adTest"), yearAreaWeight = "auto", areaDistWeight = "auto")
  expect_equal(unname(a[c("objective", "SNLL_FS", "adTest", "yearArea", "areaDist")]),
               c(41.7991463621852, 15.2772756261106, 2.93199853484711, 17.1236740922091, 6.46619810901835),
               tolerance = 1e-10)
  b <- ob(tests = c("snll_fs", "adTest"), sizeLik = "t", weighted = FALSE)
  expect_equal(unname(b[c("objective", "SNLL_FS", "adTest")]),
               c(15.5170892083938, 15.0812307018957, 0.435858506498095), tolerance = 1e-10)
  expect_equal(ob(tests = c("snll_fs", "adTest"), yearAreaWeight = "auto", areaDistWeight = "auto",
                  gateRefusal = "minLik")[["objective"]], a[["objective"]])
})
