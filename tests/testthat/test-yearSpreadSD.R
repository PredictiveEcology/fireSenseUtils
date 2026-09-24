## Per-year random effect (fitYearSpreadSD), a seasonal departure. Without it the simulated fire sizes were too
## alike: medians too large and the largest fires too small at once, and weighting big fires only shifted the whole
## distribution (2026-09-23). Each replicate draws one eps ~ N(0, sd^2) for the year, applied on the logit scale to
## every pixel of that year, inside the objective; spreadCpp() is unchanged.

## objFunInner with the spread mocked: 400 pixels, two fires with separate buffers (pixels 1:200, 201:400)
innerNoise <- function(seen, yearSpreadSD, Nreps = 2L) {
  local_mocked_bindings(
    paramsSeparate = function(...) list(logisticPars = c(0.27, 1, 1), covPars = 1),
    spreadProbFromIntegerCovs = function(...) data.table::data.table(pixelID = 1:400, cov = 0),
    logisticAll = function(...) rep(c(0.2, 0.22), each = 200),
    .package = "fireSenseUtils", .env = parent.frame())
  seen$sp <- list()
  local_mocked_bindings(
    spreadCpp = function(...) {
      a <- list(...); seen$sp[[length(seen$sp) + 1L]] <- a$spreadProb
      data.table::data.table(initialLocus = a$loci, indices = a$loci, id = seq_along(a$loci), active = FALSE)
    },
    .package = "SpaDES.tools", .env = parent.frame())
  r20 <- terra::rast(nrows = 20, ncols = 20, xmin = 0, xmax = 20, ymin = 0, ymax = 20)
  fireSenseUtils:::objFunInner(
    yr = "2004", annDTx1000 = NULL, par = 1, parsModel = NULL,
    annualFires = data.table::data.table(cells = c(50L, 350L), size = c(10, 10), ids = c(1L, 2L)),
    nonAnnualDTx1000 = NULL,
    annualFireBufferedDT = data.table::data.table(ids = rep(1:2, each = 200), pixelID = 1:400, buffer = 1L),
    indexNonAnnual = NULL, colsToUse = "cov", covMinMax = NULL,
    mutuallyExclusive = NULL, doAssertions = FALSE, maxFireSpread = 0.28,
    lowerSpreadProb = 0.13, cells = numeric(400), lanscape1stQuantileThresh = Inf,
    weighted = FALSE, r = r20, Nreps = Nreps,
    doSNLL_FSTest = TRUE, doMADTest = FALSE, doADTest = FALSE,
    plot.it = FALSE, verbose = 0, returnSims = TRUE, yearSpreadSD = yearSpreadSD)
}

test_that("sd = 0 leaves the spread probabilities as they were and draws no random numbers", {
  seen <- new.env()
  set.seed(1); before <- .Random.seed
  innerNoise(seen, yearSpreadSD = 0)
  expect_identical(.Random.seed, before)
  sp <- seen$sp[[1]]
  expect_equal(sp[c(1, 49, 51, 200)], rep(0.2, 4))
  expect_equal(sp[c(50, 350)], c(1, 1))          # ignition cells
})

test_that("with sd > 0 every pixel of the year shifts by one eps, shared by all its fires", {
  seen <- new.env()
  set.seed(2)
  innerNoise(seen, yearSpreadSD = 1, Nreps = 3L)
  for (sp in seen$sp) {
    e1 <- qlogis(sp[setdiff(1:200, 50)]) - qlogis(0.2)
    e2 <- qlogis(sp[setdiff(201:400, 350)]) - qlogis(0.22)
    expect_lt(diff(range(c(e1, e2))), 1e-9)          # one eps for the whole year: both fires share it
    expect_gt(abs(e1[1]), 1e-6)                      # and it is not zero
    expect_equal(sp[c(50, 350)], c(1, 1))            # ignition cells
  }
  ## and replicates redraw
  expect_gt(abs(qlogis(seen$sp[[1]][1]) - qlogis(seen$sp[[2]][1])), 1e-6)
})

test_that(".objfunSpreadFit takes the sd from the end of a named par, and only when it is last", {
  skip_if_not_installed("purrr")
  seen <- new.env()
  local_mocked_bindings(objFunInner = function(par, yearSpreadSD, ...) {
    seen$par <- par; seen$sd <- yearSpreadSD; list(SNLL_FS = 1)
  })
  dt <- function(...) data.table::data.table(...)
  call <- function(par, ...) fireSenseUtils:::.objfunSpreadFit(
    par = par,
    landscape = terra::rast(nrows = 10, ncols = 10, xmin = 0, xmax = 10, ymin = 0, ymax = 10, vals = 1),
    annualDTx1000 = list(year2001 = dt(pixelID = 1:2, cov1 = c(100L, 200L))),
    nonAnnualDTx1000 = list(`year2001_year2001` = dt(pixelID = 1:2, cov1 = c(100L, 200L))),
    historicalFires = list(year2001 = data.frame(size = c(100, 200), cells = 1:2)),
    fireBufferedListDT = list(year2001 = dt(pixelID = 1:2, buffer = c(1L, 0L), ids = 1L)),
    formulaToFit = "~ cov1", tests = "snll_fs", Nreps = 1L, doAssertions = FALSE, verbose = 0, ...)
  call(c(a = 0.26, b = 1, c = 1, cov1 = 2, yearSpreadSD = 0.7))
  expect_equal(seen$sd, 0.7)
  expect_identical(names(seen$par), c("a", "b", "c", "cov1"))
  call(c(0.26, 1, 1, 2))                                       # unnamed: off unless asked
  expect_equal(seen$sd, 0)
  call(c(0.26, 1, 1, 2, 0.4), fitYearSpreadSD = TRUE)           # DEoptim's unnamed par, told explicitly
  expect_equal(seen$sd, 0.4)
  expect_error(call(c(a = 0.26, yearSpreadSD = 0.7, b = 1, c = 1, cov1 = 2)), "must be the last")
})

test_that("runDEoptim fits the sd when the bounds name it, in the fit and the re-score", {
  seen <- new.env()
  lower <- stats::setNames(c(0.25, 0.2, 0.1, 0, 0), c("maxAsymptote", "hillSlope1", "inflectionPoint1", "x", "yearSpreadSD"))
  testthat::local_mocked_bindings(
    clusterSetup = function(...) list(itermax = 5, trace = FALSE, strategy = 2L, NP = 40L, cluster = NULL),
    DEoptimIterative2 = function(fn, lower, ...) { seen$fit <- list(...)$fitYearSpreadSD; list(list(member = list(pop = matrix(lower + 0.1, 1)))) },
    .package = "clusters")
  testthat::local_mocked_bindings(
    termsInDEoptim = function(...) invisible(NULL),
    rescorePopulation = function(pop, fn, reps, cl, seed = 1L, fnArgs = list()) {
      seen$rescore <- fnArgs$fitYearSpreadSD; data.table::data.table(member = 1L, rep = 1L, value = 1) })
  run <- function(lower) suppressMessages(runDEoptim(
    landscape = NULL, annualDTx1000 = NULL, nonAnnualDTx1000 = NULL, fireBufferedListDT = NULL,
    historicalFires = NULL, itermax = 5, trace = FALSE, strategy = 2L, cores = c("hostA", "hostB"),
    paths = list(cachePath = withr::local_tempdir()), lower = lower, upper = lower + 1,
    mutuallyExclusive = NULL, formulaToFit = NULL, objFunCoresInternal = 1L, covMinMax = NULL,
    maxFireSpread = 0.3, Nreps = 1L, .verbose = FALSE))
  withr::local_options(reproducible.useCache = FALSE)
  run(lower)
  expect_true(seen$fit); expect_true(seen$rescore)
  run(lower[-5])
  expect_false(seen$fit); expect_false(seen$rescore)
  expect_error(run(lower[c(1:3, 5, 4)]), "must be the last")
})
