## termsInDEoptim() names the fitted parameters; visualizeDE() draws one histogram per parameter
## of the final DEoptim population.

test_that("termsInDEoptim names the logistic terms first, then the formula's covariates", {
  ## deprecated in 0.2.3.9066; it returns the same until it is removed
  expect_warning(expect_message(
    terms <- termsInDEoptim(~ 0 + youngAge + class3, thresh = 550, numParams = 4),
    "Using a 2 parameter logistic equation"), "deprecated")
  expect_identical(terms, c("logit1", "logit2", "youngAge", "class3"))
})

fakeDE <- function(npar = 3, np = 20) {
  set.seed(1)
  structure(list(member = list(pop = matrix(stats::runif(npar * np), np, npar))),
            class = "DEoptim")
}

test_that("visualizeDE draws one histogram per parameter of the last population", {
  titles <- c("logit1", "logit2", "youngAge")
  lims <- stats::setNames(rep(0, 3), titles)
  p <- visualizeDE(list(fakeDE(), fakeDE()), titles = titles, lower = lims, upper = lims + 1)
  expect_s3_class(p, "ggplot")
  ## one panel per parameter
  expect_length(p$layers, 3L)
})

test_that("visualizeDE needs DE or a cachePath", {
  expect_error(visualizeDE(titles = "a"), "Must provide either DE or cachePath")
})

test_that("visualizeDE loads the most recent DEoptim result from the cache when DE is missing", {
  titles <- c("a", "b", "c")
  lims <- stats::setNames(rep(0, 3), titles)
  loaded <- NULL
  local_mocked_bindings(showCache = function(...) data.frame(cacheId = c("old", "new")))
  local_mocked_bindings(
    loadFromCache = function(cachePath, cacheId, ...) {
      loaded <<- cacheId
      list(fakeDE())   # as cached: one DEoptim result per block of iterations
    }, .package = "reproducible")

  expect_message(visualizeDE(cachePath = "cp", titles = titles, lower = lims, upper = lims + 1),
                 "visualizing the most recent")
  expect_identical(loaded, "new")
})

test_that("runDEoptim stops unless yearSpreadSD is the last parameter", {
  lower <- c(p1 = 0, yearSpreadSD = 0, p2 = 0)
  local_mocked_bindings(clusterSetup = function(...) list(NP = 30L), .package = "clusters")
  local_mocked_bindings(termsInDEoptim = function(...) invisible(NULL))
  expect_error(
    suppressMessages(runDEoptim(
      landscape = NULL, annualDTx1000 = NULL, nonAnnualDTx1000 = NULL,
      fireBufferedListDT = NULL, historicalFires = NULL, itermax = 5, trace = FALSE,
      strategy = 2L, cores = NA, paths = list(cachePath = withr::local_tempdir()),
      lower = lower, upper = lower + 1, mutuallyExclusive = NULL, formulaToFit = NULL,
      objFunCoresInternal = 1L, covMinMax = NULL, maxFireSpread = 0.3, Nreps = 1L,
      .verbose = FALSE)),
    "`yearSpreadSD` must be the last element")
})
