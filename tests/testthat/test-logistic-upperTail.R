## logistic3pUpper(): the spread link with Stukel's (1988) upper tail. In logistic3p the approach to the
## ceiling is set by the slope alone, so in fitted models most spreadable pixels sat pressed against the
## ceiling (72% of pixel-years in ELF 5.2.1, 2026-09-22) and no parameter could change the upper end
## without moving everything else. `upperTail1` does only that, and 0 is the old link.

test_that("upperTail() matches sirt::pgenlogis(), the reference implementation of Stukel's link", {
  ## plogis(upperTail(x, a)) is pgenlogis(x, alpha1 = a, alpha2 = 0). Values computed with
  ## sirt 4.2-133's pgenlogis() (R/pgenlogis.R, sourced verbatim); sirt is not a dependency.
  x <- c(-3, -0.5, 0, 0.5, 2, 6, 15)
  ref <- list(
    "-0.8" = c(0.0474258731775668, 0.377540668798145, 0.5, 0.603624493827386, 0.767525168871424,
               0.900008803419917, 0.961066430928329),
    "-0.3" = c(0.0474258731775668, 0.377540668798145, 0.5, 0.614406520325504, 0.827309624272457,
               0.96869174742818, 0.996606498207396),
    "0"    = c(0.0474258731775668, 0.377540668798145, 0.5, 0.622459331201855, 0.880797077977882,
               0.997527376843365, 0.999999694097773),
    "0.4"  = c(0.0474258731775668, 0.377540668798145, 0.5, 0.634948834165064, 0.955376817340328,
               0.999999999986894, 1))
  for (a in names(ref))
    expect_equal(stats::plogis(fireSenseUtils:::upperTail(x, as.numeric(a))), ref[[a]],
                 tolerance = 1e-12, info = paste("alpha1 =", a))
})

test_that("upperTail1 = 0 is logistic3p, so the new link nests the old one", {
  x <- seq(-10, 30, by = 0.25)
  for (par in list(c(0.26, 0.58, 0.17), c(0.252, 0.98, 3.9), c(0.27, 1.5, 1)))
    expect_equal(logistic3pUpper(x, c(par, 0), par1 = 0.13), logistic3p(x, par, par1 = 0.13),
                 tolerance = 1e-12)
})

test_that("a negative upper tail lowers the top of the curve and leaves the lower half alone", {
  par <- c(0.26, 0.58, 0.17)
  low <- seq(-10, 0, by = 0.5); high <- seq(5, 30, by = 1)
  expect_equal(logistic3pUpper(low, c(par, -0.5), 0.13), logistic3pUpper(low, c(par, 0), 0.13))
  expect_true(all(logistic3pUpper(high, c(par, -0.5), 0.13) < logistic3pUpper(high, c(par, 0), 0.13)))
  ## ...and a positive one raises it
  expect_true(all(logistic3pUpper(high, c(par, 0.5), 0.13) > logistic3pUpper(high, c(par, 0), 0.13)))
})

test_that("the curve stays increasing and inside [floor, ceiling] across the parameter box", {
  x <- seq(-20, 60, by = 0.1)
  for (a in c(-1, -0.3, 0, 0.3, 1)) for (h in c(0.2, 1, 2)) for (cc in c(0.1, 1, 4)) {
    p <- logistic3pUpper(x, c(0.26, h, cc, a), par1 = 0.13)
    lab <- sprintf("a %g h %g c %g", a, h, cc)
    expect_true(all(diff(p) >= -1e-15), info = lab)
    expect_true(all(p >= 0.13 - 1e-12 & p <= 0.26 + 1e-12), info = lab)
  }
})

test_that("logisticAll() dispatches the upper-tail link by name or by `link`, and nothing else changes", {
  mat <- matrix(c(-2, 0, 1, 4, 12), ncol = 1); cp <- 1
  upper <- c(maxAsymptote = 0.26, hillSlope1 = 0.58, inflectionPoint1 = 0.17, upperTail1 = -0.4)
  want <- logistic3pUpper(mat %*% cp, upper, par1 = 0.13)
  expect_equal(logisticAll(upper, mat, cp, 0.13), want)                                   # by name
  expect_equal(logisticAll(unname(upper), mat, cp, 0.13, link = "logistic3pUpper"), want) # explicit
  ## an unnamed 3-vector is logistic3p, as before
  expect_equal(logisticAll(unname(upper[1:3]), mat, cp, 0.13),
               logistic3p(mat %*% cp, unname(upper[1:3]), par1 = 0.13))
  expect_error(logisticAll(unname(upper), mat, cp, 0.13, link = "cloglog"), "should be")
})

test_that("objFunInner() hands `link` to logisticAll()", {
  seen <- new.env()
  local_mocked_bindings(
    paramsSeparate = function(...) list(logisticPars = c(0.27, 1, 1, -0.3), covPars = 1),
    spreadProbFromIntegerCovs = function(...) data.table::data.table(pixelID = 1:400, cov = 0),
    logisticAll = function(...) { seen$link <- list(...)$link; rep(0.2, 400) },
    .package = "fireSenseUtils")
  local_mocked_bindings(
    spreadCpp = function(...) data.table::data.table(initialLocus = list(...)$loci, indices = 1:10,
                                                    id = 1L, active = FALSE),
    .package = "SpaDES.tools")
  r20 <- terra::rast(nrows = 20, ncols = 20, xmin = 0, xmax = 20, ymin = 0, ymax = 20)
  fireSenseUtils:::objFunInner(
    yr = "2004", annDTx1000 = NULL, par = 1, parsModel = NULL,
    annualFires = data.table::data.table(cells = 210L, size = 10, ids = 1L),
    nonAnnualDTx1000 = NULL,
    annualFireBufferedDT = data.table::data.table(ids = 1L, pixelID = 1:400, buffer = 1L),
    indexNonAnnual = NULL, colsToUse = "cov", covMinMax = NULL,
    mutuallyExclusive = NULL, doAssertions = FALSE, maxFireSpread = 0.28,
    lowerSpreadProb = 0.13, cells = numeric(400), lanscape1stQuantileThresh = 0.25,
    weighted = FALSE, r = r20, Nreps = 2L,
    doSNLL_FSTest = TRUE, doMADTest = FALSE, doADTest = FALSE,
    plot.it = FALSE, verbose = 0, link = "logistic3pUpper")
  expect_identical(seen$link, "logistic3pUpper")
})

test_that("runDEoptim() hands `link` to the fit and to the re-score", {
  seen <- new.env()
  finalPop <- matrix(c(0.26, 1, 1, -0.3, 0.5), nrow = 1)
  lower <- stats::setNames(c(0.25, 0.2, 0.1, -1, 0),
                           c("maxAsymptote", "hillSlope1", "inflectionPoint1", "upperTail1", "x"))
  testthat::local_mocked_bindings(
    clusterSetup = function(...) list(itermax = 5, trace = FALSE, strategy = 2L, NP = 40L, cluster = NULL),
    DEoptimIterative2 = function(fn, lower, upper, control, ...) {
      seen$fitLink <- list(...)$link
      list(list(member = list(pop = finalPop)))
    },
    .package = "clusters")
  testthat::local_mocked_bindings(
    termsInDEoptim = function(...) invisible(NULL),
    rescorePopulation = function(pop, fn, reps, cl, seed = 1L, fnArgs = list()) {
      seen$rescoreLink <- fnArgs$link
      data.table::data.table(member = 1L, rep = 1L, value = 1)
    })
  withr::local_options(reproducible.useCache = FALSE)
  suppressMessages(runDEoptim(
    landscape = NULL, annualDTx1000 = NULL, nonAnnualDTx1000 = NULL,
    fireBufferedListDT = NULL, historicalFires = NULL,
    itermax = 5, trace = FALSE, strategy = 2L, cores = c("hostA", "hostB"),
    paths = list(cachePath = withr::local_tempdir()),
    lower = lower, upper = lower + 1, mutuallyExclusive = NULL, formulaToFit = NULL,
    objFunCoresInternal = 1L, covMinMax = NULL, maxFireSpread = 0.3, Nreps = 1L,
    .verbose = FALSE, link = "logistic3pUpper"))
  expect_identical(seen$fitLink, "logistic3pUpper")
  expect_identical(seen$rescoreLink, "logistic3pUpper")
})
