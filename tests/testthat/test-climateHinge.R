## climateHinge*(), addClimateHinge(), ledgerClimateHingeKnots(): the climate covariate with a hinge term at the ELF's median

refTable <- function() {
  data.table::data.table(year = rep(c("year1990", "year1991"), each = 4), pixelID = rep(1:4, 2),
                         CMD = c(10, 20, 30, 40, 100, 110, 120, 130), Tmaxsm = c(1, 2, 3, 4, 5, 6, 7, NA))
}

test_that("the term is named <var>_hinge, per variable", {
  expect_identical(climateHingeCols(c("CMD", "Tmaxsm")), c("CMD_hinge", "Tmaxsm_hinge"))
  expect_identical(climateHingeBase(c("CMD_hinge", "CMD_yr", "CMD")), c("CMD", "CMD_yr", "CMD"))
  expect_identical(isClimateHingeCol(c("CMD_hinge", "CMD_yr", "CMD", "class1")), c(TRUE, FALSE, FALSE, FALSE))
})

test_that("the knot is the median over the reference rows of all years, NA left out", {
  knots <- climateHingeKnots(refTable(), c("CMD", "Tmaxsm"))
  expect_identical(names(knots), c("CMD", "Tmaxsm"))
  expect_identical(knots[["CMD"]], 70)
  expect_identical(knots[["Tmaxsm"]], 4)
})

test_that("the hinge is 0 below the knot and C - k above it; C stays; NA stays NA", {
  dt <- refTable()
  addClimateHinge(dt, c(CMD = 70, Tmaxsm = 4))
  expect_identical(dt$CMD_hinge, c(0, 0, 0, 0, 30, 40, 50, 60))
  expect_identical(dt$Tmaxsm_hinge, c(0, 0, 0, 0, 1, 2, 3, NA))
  expect_identical(dt$CMD, refTable()$CMD)
})

test_that("fit and prediction agree: the stored knot rebuilds the fit's hinge on a one-year table", {
  fit <- refTable()
  knots <- climateHingeKnots(fit, "CMD")
  addClimateHinge(fit, knots)
  ## the ledger row stores the knots as a list; the prediction's table is one year's layer
  sa <- data.frame(polygonID = c("A", "B"))
  sa[[spreadFitClimateHingeKnotsTxt]] <- I(list(as.list(knots), NULL))
  pred <- refTable()[refTable()$year == "year1991", c("pixelID", "CMD")]
  addClimateHinge(pred, ledgerClimateHingeKnots(sa, 1L))
  expect_identical(pred$CMD_hinge, fit$CMD_hinge[fit$year == "year1991"])
})

test_that("a row without knots (no column, or an empty entry) gives NULL", {
  expect_null(ledgerClimateHingeKnots(data.frame(polygonID = "A"), 1L))
  sa <- data.frame(polygonID = "A")
  sa[[spreadFitClimateHingeKnotsTxt]] <- I(list(NULL))
  expect_null(ledgerClimateHingeKnots(sa, 1L))
})

test_that("the fixed range is 0 to half the variable's width; given ranges win", {
  r <- climateHingeRanges(list(CMD = c(0, 100), Tmaxsm = c(0, 10)))
  expect_identical(r$CMD_hinge, c(0, 50))
  expect_identical(r$Tmaxsm_hinge, c(0, 5))
  expect_setequal(names(climateHingeRanges(climateCovRanges)),
                  c(names(climateCovRanges), climateHingeCols(names(climateCovRanges))))
  expect_identical(climateHingeRanges(list()), list())
  expect_identical(climateHingeRanges(list(CMD = c(0, 100), CMD_hinge = c(0, 1)))$CMD_hinge, c(0, 1))
})
