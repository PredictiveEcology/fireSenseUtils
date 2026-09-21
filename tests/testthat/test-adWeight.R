## `adWeight` multiplies the Anderson-Darling statistic before it is added to the fire-size SNLL.
## The fixed 50 balanced the two terms only for the kernel-density likelihood with no size weight:
## measured across six ELFs, the AD term's share of the influence on the objective ranged from 0.16
## to 0.94 over the four sizeLik x weighted combinations. `adWeight = "auto"` scales with the number
## of fitted fires, which is what both terms grow with, and holds that share to 0.36-0.68.

test_that("adWeightAuto() is c * sqrt(nFires), with c from sizeLik and weighted", {
  expect_equal(adWeightAuto(100, "kde", FALSE), 1.897 * 10)
  expect_equal(adWeightAuto(100, "kde", "sqrt"), 5.098 * 10)
  expect_equal(adWeightAuto(100, "t", FALSE), 0.282 * 10)
  expect_equal(adWeightAuto(100, "t", "sqrt"), 0.804 * 10)
  ## TRUE means the log(size) weight, which was measured with the unweighted constants
  expect_equal(adWeightAuto(400, "kde", TRUE), adWeightAuto(400, "kde", FALSE))
  expect_gt(adWeightAuto(1624, "kde", FALSE), adWeightAuto(138, "kde", FALSE))
})

## the fixture of test-objFunSpread-pruneAboveBehaviour.R: 3 fire years, 6 fires, objFunInner mocked
mkFixture <- function() {
  dt <- function(...) data.table::data.table(...)
  yrs <- c("2001", "2002", "2003")
  list(
    par = c(1, 1),
    landscape = terra::rast(nrows = 10, ncols = 10, xmin = 0, xmax = 10, ymin = 0, ymax = 10, vals = 1),
    annualDTx1000 = stats::setNames(lapply(yrs, function(y) dt(pixelID = 1:2, cov1 = c(100L, 200L))), yrs),
    nonAnnualDTx1000 = list(`2001_2003` = dt(pixelID = 1:2, cov1 = c(100L, 200L))),
    historicalFires = list(`2001` = data.frame(size = c(100, 200), cells = 1:2),
                           `2002` = data.frame(size = c(150, 250), cells = 3:4),
                           `2003` = data.frame(size = c(5, 6), cells = 5:6)),
    fireBufferedListDT = stats::setNames(lapply(yrs, function(y) dt(pixelID = 1:2, buffer = c(1L, 0L), ids = 1L)), yrs),
    formulaToFit = "~ cov1"
  )
}
runObjFun <- function(...) {
  fx <- mkFixture()
  local_mocked_bindings(
    objFunInner = function(...) list(SNLL_FS = 0, allFireSizes = c(2, 10, 300)),
    adStatistic = function(x, y) 7
  )
  do.call(fireSenseUtils:::.objfunSpreadFit, c(
    fx[c("par", "landscape", "annualDTx1000", "nonAnnualDTx1000", "formulaToFit",
         "historicalFires", "fireBufferedListDT")],
    list(tests = c("adTest", "snll_fs"), Nreps = 1L, doAssertions = FALSE, plot.it = FALSE,
         verbose = 0, thresh = Inf, ...)
  ))
}

test_that("the default scales with the fitted fires, and follows sizeLik and weighted", {
  skip_if_not_installed("terra")
  ## 6 fires in the fixture; minFireSize = 2 drops none of them
  expect_equal(runObjFun(), adWeightAuto(6, "kde", TRUE) * 7)
  expect_equal(runObjFun(), runObjFun(adWeight = "auto"))
  expect_equal(runObjFun(sizeLik = "t", weighted = "sqrt"), adWeightAuto(6, "t", "sqrt") * 7)
})

test_that("a number is used as given, including the former fixed default", {
  skip_if_not_installed("terra")
  expect_equal(runObjFun(adWeight = 3), 3 * 7)
  expect_equal(runObjFun(adWeight = 50), 50 * 7)
})
