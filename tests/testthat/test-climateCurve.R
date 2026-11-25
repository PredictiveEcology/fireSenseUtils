## softplus(), climateCurve*(), addClimateCurve(), ledgerClimateCurvePars(): the climate covariate with a smooth curve term
## at the ELF's median, with the ELF's sd as the scale of the bend

refTable <- function() {
  data.table::data.table(year = rep(c("year1990", "year1991"), each = 4), pixelID = rep(1:4, 2),
                         CMD = c(10, 20, 30, 40, 100, 110, 120, 130), Tmaxsm = c(1, 2, 3, 4, 5, 6, 7, NA))
}

test_that("softplus is s * log(1 + exp(x / s)), smooth, near 0 below and near x above, and does not overflow", {
  x <- c(-50, -2, 0, 2, 50)
  expect_equal(softplus(x, 3), 3 * log(1 + exp(x / 3)))
  expect_equal(softplus(0, 3), 3 * log(2))
  expect_lt(softplus(-50, 3), 1e-6)
  expect_equal(softplus(50, 3), 50, tolerance = 1e-6)
  expect_identical(softplus(1e5, 1), 1e5)                     # exp(1e5) would be Inf
  expect_true(all(diff(softplus(seq(-10, 10, 0.5), 2)) > 0))   # strictly increasing, no corner
  expect_identical(softplus(NA_real_), NA_real_)
  expect_identical(softplus(1:3, s = 1e-9), as.numeric(1:3))   # a tiny scale is pmax(0, x)
})

test_that("the term is named <var>_curve, per variable", {
  expect_identical(climateCurveCols(c("CMD", "Tmaxsm")), c("CMD_curve", "Tmaxsm_curve"))
  expect_identical(climateCurveBase(c("CMD_curve", "CMD_yr", "CMD")), c("CMD", "CMD_yr", "CMD"))
  expect_identical(isClimateCurveCol(c("CMD_curve", "CMD_yr", "CMD", "class1")), c(TRUE, FALSE, FALSE, FALSE))
})

test_that("the centre is the median and the scale the sd over the reference rows of all years, NA left out", {
  pars <- climateCurvePars(refTable(), c("CMD", "Tmaxsm"))
  expect_identical(names(pars), c("CMD", "Tmaxsm"))
  expect_identical(pars$CMD, c(k = 70, s = stats::sd(refTable()$CMD)))
  expect_identical(pars$Tmaxsm, c(k = 4, s = stats::sd(1:7)))
})

test_that("a constant climate has no curve", {
  ref <- data.table::data.table(CMD = rep(5, 4))
  expect_error(climateCurvePars(ref, "CMD"), "scale `s` above 0 for CMD")
})

test_that("the term is softplus(C - k, s); C stays; NA stays NA", {
  dt <- refTable()
  addClimateCurve(dt, list(CMD = c(k = 70, s = 30), Tmaxsm = c(k = 4, s = 2)))
  expect_equal(dt$CMD_curve, 30 * log1p(exp((refTable()$CMD - 70) / 30)))
  expect_equal(dt$Tmaxsm_curve, 2 * log1p(exp((refTable()$Tmaxsm - 4) / 2)))
  expect_true(is.na(dt$Tmaxsm_curve[8]))
  expect_identical(dt$CMD, refTable()$CMD)
  expect_true(all(dt$CMD_curve > 0))                         # smooth: never exactly 0
  expect_gt(dt$CMD_curve[8] - dt$CMD_curve[7], dt$CMD_curve[2] - dt$CMD_curve[1])   # steeper above the centre
})

test_that("fit and prediction agree: the stored pars rebuild the fit's term on a one-year table", {
  fit <- refTable()
  pars <- climateCurvePars(fit, "CMD")
  addClimateCurve(fit, pars)
  sa <- data.frame(polygonID = c("A", "B"))
  sa[[spreadFitClimateCurveParsTxt]] <- I(list(pars, NULL))
  pred <- refTable()[refTable()$year == "year1991", c("pixelID", "CMD")]
  addClimateCurve(pred, ledgerClimateCurvePars(sa, 1L))
  expect_identical(pred$CMD_curve, fit$CMD_curve[fit$year == "year1991"])
})

test_that("a row without pars (no column, or an empty entry) gives NULL", {
  expect_null(ledgerClimateCurvePars(data.frame(polygonID = "A"), 1L))
  sa <- data.frame(polygonID = "A")
  sa[[spreadFitClimateCurveParsTxt]] <- I(list(NULL))
  expect_null(ledgerClimateCurvePars(sa, 1L))
})

test_that("the fixed range is 0 to half the variable's width; given ranges win", {
  r <- climateCurveRanges(list(CMD = c(0, 100), Tmaxsm = c(0, 10)))
  expect_identical(r$CMD_curve, c(0, 50))
  expect_identical(r$Tmaxsm_curve, c(0, 5))
  expect_setequal(names(climateCurveRanges(climateCovRanges)),
                  c(names(climateCovRanges), climateCurveCols(names(climateCovRanges))))
  expect_identical(climateCurveRanges(list()), list())
  expect_identical(climateCurveRanges(list(CMD = c(0, 100), CMD_curve = c(0, 1)))$CMD_curve, c(0, 1))
})
