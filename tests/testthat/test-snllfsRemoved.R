## The "SNLL_FS" test (the per-fire likelihood objective) and the first-block early stop that went with it
## (`thresh`, `pruneAbove`, `firstBlockPerFire`, `bailed`) were removed; "crps" is the default objective.

test_that("the spread objective and runDEoptim() default to tests = \"crps\"", {
  expect_identical(eval(formals(fireSenseUtils:::.objfunSpreadFit)$tests), "crps")
  expect_identical(eval(formals(runDEoptim)$tests), "crps")
})

test_that("asking for SNLL_FS (or SNLL) stops, and says the test was removed", {
  for (tst in list("SNLL_FS", "snll_fs", c("adTest", "SNLL_FS"), "SNLL")) {
    expect_error(fireSenseUtils:::.objfunSpreadFit(tests = tst), "were removed")
    expect_error(runDEoptim(tests = tst), "were removed")   # on the master, before any cluster is built
  }
})

test_that("thresh and pruneAbove are no longer arguments; thresh warns once, then is ignored", {
  expect_false(any(c("thresh", "pruneAbove") %in% names(formals(fireSenseUtils:::.objfunSpreadFit))))
  expect_false("thresh" %in% names(formals(runDEoptim)))
  expect_true("thresh" %in% names(fireSenseUtils:::retiredObjectiveArgs))
  rm(list = ls(fireSenseUtils:::.deprecatedArgsWarned), envir = fireSenseUtils:::.deprecatedArgsWarned)
  expect_warning(fireSenseUtils:::deprecatedArgs("thresh"), "deprecated and ignored")
})

test_that("a crps fit scores every year and reports no first-block terms", {
  dt <- function(...) data.table::data.table(...)
  yrs <- c("2001", "2002", "2003")
  seen <- character()
  local_mocked_bindings(
    objFunInner = function(yr, ...) {
      seen <<- c(seen, yr)
      list(crpsFire = 1, crpsSizes = dt(year = yr, rep = 1L, N = 5L))
    },
    .package = "fireSenseUtils")
  out <- fireSenseUtils:::.objfunSpreadFit(
    par = c(1, 1),
    landscape = terra::rast(nrows = 10, ncols = 10, xmin = 0, xmax = 10, ymin = 0, ymax = 10, vals = 1),
    annualDTx1000 = stats::setNames(lapply(yrs, function(y) dt(pixelID = 1:2, cov1 = c(100L, 200L))), yrs),
    nonAnnualDTx1000 = list(`2001_2003` = dt(pixelID = 1:2, cov1 = c(100L, 200L))),
    historicalFires = list(`2001` = data.frame(size = c(100, 200), cells = 1:2),
                           `2002` = data.frame(size = c(150, 250), cells = 3:4),
                           `2003` = data.frame(size = c(5, 6), cells = 5:6)),
    fireBufferedListDT = stats::setNames(lapply(yrs, function(y) dt(pixelID = 1:2, buffer = c(1L, 0L), ids = 1L)), yrs),
    formulaToFit = "~ 0 + cov1", tests = "crps", Nreps = 1L, doAssertions = FALSE, verbose = 0, returnTerms = TRUE)
  expect_setequal(seen, yrs)
  expect_false(any(c("SNLL_FS", "firstBlockPerFire", "bailed") %in% names(out)))
  expect_lt(out[["objective"]], 1e6)
})
