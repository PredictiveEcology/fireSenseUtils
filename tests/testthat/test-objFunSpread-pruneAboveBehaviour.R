## Behaviour of `pruneAbove`, exercising the real code path rather than inspecting source.
##
## The other pruneAbove file asserts the argument and the bound by source inspection, which is this
## package's usual idiom -- but it never EXECUTES
## `threshold <- min(thresh * numYrsDone, pruneAbove)`, so the semantics went untested (and patch
## coverage for the change was 0%).
##
## Reaching that line normally needs a full spread fit: landscape, covariates, buffered fires and
## Nreps spread simulations per year. It does not have to. `.objfunSpreadFit()` dispatches the years
## of a batch through `purrr::pmap(.f = objFunInner)`, and `objFunInner` is where all the simulation
## lives; everything the early-bail arithmetic needs is in what it RETURNS. Mocking it -- the idiom
## test-objFunInner-noTests.R already uses for this chain -- makes the whole bail path testable with
## a synthetic fixture and no spread() calls.
##
## The fixture gives three fire years. `lrgSmallFireYears` splits them into the two largest (batch 1,
## the only batch with a checkpoint) and the rest (batch 2), so `numYrsDone` is 2 and the static
## bound is `thresh * 2`. Choosing a mocked SNLL_FS that sits BELOW the static bound but ABOVE
## `pruneAbove` isolates the new behaviour: with the default `Inf` the evaluation continues, and only
## the caller-supplied bound stops it.
##
## Note `.objfunSpreadFit()` reads `results$SNLL` while `objFunInner` returns `SNLL_FS`; that resolves
## by `$` partial matching on a list. The mock therefore returns `SNLL_FS`, as the real function does.

mkFixture <- function() {
  skip_if_not_installed("terra")
  skip_if_not_installed("purrr")
  dt <- function(...) data.table::data.table(...)
  yrs <- c("2001", "2002", "2003")
  list(
    par = c(1, 1),
    landscape = terra::rast(nrows = 10, ncols = 10, xmin = 0, xmax = 10, ymin = 0, ymax = 10, vals = 1),
    annualDTx1000 = stats::setNames(lapply(yrs, function(y) dt(pixelID = 1:2, cov1 = c(100L, 200L))), yrs),
    nonAnnualDTx1000 = list(`2001_2003` = dt(pixelID = 1:2, cov1 = c(100L, 200L))),
    ## 2002 (400) and 2001 (300) are the two largest, so they are batch 1; 2003 (11) is batch 2
    historicalFires = list(`2001` = data.frame(size = c(100, 200), cells = 1:2),
                           `2002` = data.frame(size = c(150, 250), cells = 3:4),
                           `2003` = data.frame(size = c(5, 6),     cells = 5:6)),
    fireBufferedListDT = stats::setNames(lapply(yrs, function(y) dt(pixelID = 1:2, buffer = c(1L, 0L), ids = 1L)), yrs),
    formulaToFit = "~ cov1"
  )
}

runObjFun <- function(fx, pruneAbove, snllPerYear, thresh = 550) {
  do.call(fireSenseUtils:::.objfunSpreadFit, c(
    fx[c("par", "landscape", "annualDTx1000", "nonAnnualDTx1000", "formulaToFit",
         "historicalFires", "fireBufferedListDT")],
    list(tests = "snll_fs", Nreps = 1L, doAssertions = FALSE, plot.it = FALSE, verbose = 0,
         thresh = thresh, pruneAbove = pruneAbove)
  ))
}

test_that("with the default Inf, a batch under the static threshold is NOT pruned", {
  fx <- mkFixture()
  ## 2 years x 450 = 900, below the static bound of 550 * 2 = 1100, so the fit continues and the
  ## objective is a real value rather than the 1e6 fail sentinel.
  local_mocked_bindings(objFunInner = function(...) list(SNLL_FS = 450))
  out <- runObjFun(fx, pruneAbove = Inf, snllPerYear = 450)
  expect_true(is.numeric(out))
  expect_lt(out, 1e6)
})

test_that("a caller-supplied pruneAbove below the batch's value DOES prune it", {
  fx <- mkFixture()
  ## Same 900, still below the static 1100 -- so only pruneAbove can stop it. 800 < 900 does.
  local_mocked_bindings(objFunInner = function(...) list(SNLL_FS = 450))
  out <- runObjFun(fx, pruneAbove = 800, snllPerYear = 450)
  expect_equal(out, 1e6)
})

test_that("pruneAbove above the batch's value leaves it alone", {
  ## Guard on the direction of the comparison: a bound the trial is already better than must not
  ## prune it, or every trial would be discarded and the search would collapse.
  fx <- mkFixture()
  local_mocked_bindings(objFunInner = function(...) list(SNLL_FS = 450))
  out <- runObjFun(fx, pruneAbove = 5000, snllPerYear = 450)
  expect_lt(out, 1e6)
})

test_that("the static threshold still prunes on its own, as it did before pruneAbove existed", {
  ## Regression: 2 x 700 = 1400 exceeds 550 * 2 = 1100 with pruneAbove left at its Inf default,
  ## so behaviour without a caller-supplied bound is unchanged.
  fx <- mkFixture()
  local_mocked_bindings(objFunInner = function(...) list(SNLL_FS = 700))
  out <- runObjFun(fx, pruneAbove = Inf, snllPerYear = 700)
  expect_equal(out, 1e6)
})
