## spreadClimate = "activation" (spreadLinkInputs(activationIntercept =)): fuel sets each pixel's own ceiling,
## pc_i = f + (p50 - f) * plogis(eta_i), eta_i the linear predictor of the non-climate columns, and climate is the
## linear predictor, x_i = a + sum(beta_clim * (clim - climateRef)); p_i = powerLawSmooth(x_i; pc_i), lower = floor.
## Fixture: helper-spreadLinkFixture.R (f = 0.13, lp3, matOf(), climFix, spreadProbSeenC()).

p50 <- 0.3
a <- 1.2
b <- c(clim = 1.5, fuel = 3, youngAge = -1)
ref <- list(clim = 0.3)
gam <- 1.9
shp <- 0.5
actInputs <- function(mat, ...) spreadLinkInputs(mat, lp3, b, f, "powerLawSmooth", ref, p50, powerLawGamma = gam,
                                                 powerLawShape = shp, lowerAsymptote = f, activationIntercept = a, ...)

test_that("p of each pixel is the formula: pc_i from the fuel terms, x_i from climate and the activation intercept", {
  mat <- matOf("2005")
  inputs <- actInputs(mat)
  eta <- unname(mat[, "fuel"] * b[["fuel"]] + mat[, "youngAge"] * b[["youngAge"]])
  pc <- f + (p50 - f) * plogis(eta)
  x <- unname(a + b[["clim"]] * (mat[, "clim"] - ref$clim))
  expect_equal(unname(inputs$pars[, "pc"]), pc)
  expect_equal(inputs$x, x)
  expect_true(all(inputs$pars[, "floor"] == f) && all(inputs$pars[, "lower"] == f))
  expect_true(all(inputs$pars[, "gamma"] == gam) && all(inputs$pars[, "shape"] == shp))
  ## by hand: the power law above the floor, the floor below it (lower = floor)
  byHand <- ifelse(x >= 0, pc - (pc - f) * exp(-(pmax(x, 0) / gam)^shp), f)
  p <- spreadProbFromCovs(mat, lp3, b, f, link = "powerLawSmooth", climateRef = ref, powerLawPc50 = p50,
                          powerLawGamma = gam, powerLawShape = shp, lowerAsymptote = f, activationIntercept = a)
  expect_equal(as.numeric(p), byHand)
  expect_identical(as.numeric(p), linearLinkProb(inputs$x, "powerLawSmooth", inputs$pars))
  ## pc0 (logisticPars[1]) plays no part
  expect_identical(spreadLinkInputs(mat, c(0.2, 1, 1), b, f, "powerLawSmooth", ref, p50, powerLawGamma = gam,
                                    powerLawShape = shp, lowerAsymptote = f, activationIntercept = a)$x, inputs$x)
  ## at the reference climate x is the activation intercept
  expect_equal(actInputs(mat, climateAtRef = TRUE)$x, rep(a, nrow(mat)))
})

test_that("the ratio of p - floor between two pixels does not change with climate", {
  mat <- matOf("2004")[1:6, ]
  mat[, "fuel"] <- c(0.05, 0.1, 0.2, 0.4, 0.6, 0.9)
  mat[, "youngAge"] <- 0
  ratios <- sapply(c(0.35, 0.6, 1, 2, 4), function(clim) {
    m <- mat; m[, "clim"] <- clim
    p <- as.numeric(linearLinkProb(actInputs(m)$x, "powerLawSmooth", actInputs(m)$pars))
    expect_true(all(p < actInputs(m)$pars[, "pc"]))   # p is still below every pixel's ceiling
    (p - f) / (p[1] - f)
  })
  for (i in seq_len(nrow(ratios))) expect_equal(ratios[i, ], rep(ratios[i, 1], ncol(ratios)))
  ## and it is the ratio of the fuel potentials
  expect_equal(ratios[, 1], plogis(3 * mat[, "fuel"]) / plogis(3 * mat[1, "fuel"]), ignore_attr = TRUE)
})

test_that("p never exceeds p50, whatever the fuel and the climate", {
  mat <- matOf("2005")
  for (clim in c(-50, -1, 0.3, 2, 50, 1e3)) for (fuelBy in c(-100, 0, 1, 100)) {
    m <- mat; m[, "clim"] <- clim; m[, "fuel"] <- m[, "fuel"] * fuelBy
    inputs <- actInputs(m)
    expect_true(all(inputs$pars[, "pc"] <= p50))
    p <- linearLinkProb(inputs$x, "powerLawSmooth", inputs$pars)
    expect_lte(max(p), p50)
    expect_gte(min(p), f)
  }
})

test_that("rows with no climate or no fuel cannot burn", {
  mat <- matOf("2005")
  mat[3, "clim"] <- NA
  mat[7, "fuel"] <- NA
  inputs <- actInputs(mat)
  expect_false(anyNA(inputs$pars))
  expect_true(all(is.na(inputs$x[c(3, 7)])))
  expect_false(anyNA(inputs$x[-c(3, 7)]))
  expect_identical(inputs$pars[c(3, 7), , drop = FALSE], unburnableLinkPars("powerLawSmooth", 2L))
  expect_true(all(is.na(spreadProbFromCovs(mat, lp3, b, f, link = "powerLawSmooth", climateRef = ref, powerLawPc50 = p50,
                                           lowerAsymptote = f, activationIntercept = a)[c(3, 7)])))
})

test_that("the activation stops without what it needs", {
  mat <- matOf("2004")
  call <- function(...) {
    args <- utils::modifyList(list(mat = mat, logisticPars = lp3, covPars = b, lowerSpreadProb = f, link = "powerLawSmooth",
                                   climateRef = ref, powerLawPc50 = p50, lowerAsymptote = f, activationIntercept = a),
                              list(...))
    do.call(spreadLinkInputs, args)
  }
  expect_error(call(climateRef = NULL), "needs `climateRef`")
  expect_error(call(powerLawPc50 = NA_real_), "powerLawPc50")
  expect_error(call(link = "powerLaw"), "powerLawSmooth")
  expect_error(call(link = "richards"), "power-law link")
  expect_error(call(lowerAsymptote = 0.05), "lowerAsymptote = lowerSpreadProb")
  expect_error(call(activationIntercept = NA_real_), "one finite number")
  ## a logistic link cannot take it either
  expect_error(spreadProbFromCovs(mat, lp3, b, f, link = NULL, climateRef = ref, powerLawPc50 = p50,
                                  activationIntercept = a), "must be one of")
  ## the same check the module makes of its parameters and of a ledger row
  expect_true(checkLinkParams("powerLawSmooth", f, activation = TRUE, floor = f))
  expect_error(checkLinkParams("powerLawSmooth", 0.05, activation = TRUE), "lowerAsymptote = lowerSpreadProb")
  expect_error(checkLinkParams("powerLaw", f, activation = TRUE, floor = f), "powerLawSmooth")
  expect_true(checkLinkParams("powerLaw", 5))   # without the activation, nothing new is checked
  ## pc0 above p50 is no error here: it is not used
  expect_silent(call(logisticPars = c(0.4, 1, 1)))
})

test_that("splitSpreadPar() takes the activation intercept out and puts p50 first as maxAsymptote", {
  split <- fireSenseUtils:::splitSpreadPar
  ## as DEoptim gives it: unnamed, covariates, the intercept, the year sd
  s <- split(c(1.5, 3, -1, a, 0.7), fitYearSpreadSD = TRUE, fitActivationIntercept = TRUE, powerLawPc50 = p50)
  expect_identical(unname(s$par), c(p50, 1, 1, 1.5, 3, -1))
  expect_identical(s$activationIntercept, a)
  expect_identical(s$yearSpreadSD, 0.7)
  ## named, as drawActivePars() and the validation give it: found by name
  nm <- c(clim = 1.5, fuel = 3, youngAge = -1, activationIntercept = a, yearSpreadSD = 0.7)
  s2 <- split(nm, powerLawPc50 = p50)
  expect_identical(unname(s2$par), unname(s$par))
  expect_identical(names(s2$par)[1:3], c("maxAsymptote", "hillSlope1", "inflectionPoint1"))
  ## a ledger row's set already has maxAsymptote (p50, with the fixed parameters): kept, not doubled
  row <- c(maxAsymptote = p50, hillSlope1 = 1, inflectionPoint1 = 1, nm)
  expect_identical(split(row)$par, s2$par)
  expect_identical(split(row)$activationIntercept, a)
  expect_error(split(c(nm[c(1, 2, 4, 3)], yearSpreadSD = 0.7), powerLawPc50 = p50), "right before the random-effect")
  expect_error(split(nm), "powerLawPc50")
  ## without it nothing changes, and the element is NULL
  old <- split(c(0.27, 1.5, 3, -1, 0.7), fitYearSpreadSD = TRUE)
  expect_identical(unname(old$par), c(0.27, 1, 1, 1.5, 3, -1))
  expect_null(old$activationIntercept)
})

test_that("spreadRandomEffectArgs() sets fitActivationIntercept only for bounds that have it", {
  lo <- c(clim = 0, fuel = 0, activationIntercept = -10, yearSpreadSD = 0)
  expect_identical(spreadRandomEffectArgs(lo), list(fitYearSpreadSD = TRUE, fitActivationIntercept = TRUE))
  expect_identical(spreadRandomEffectArgs(lo[-4]), list(fitYearSpreadSD = FALSE, fitActivationIntercept = TRUE))
  ## the other fits get what they got
  expect_identical(spreadRandomEffectArgs(c(maxAsymptote = 0.14, clim = 0, yearSpreadSD = 0)), list(fitYearSpreadSD = TRUE))
  expect_error(spreadRandomEffectArgs(lo[c(3, 1, 2, 4)]), "right before the random-effect sds")
})

test_that("the objective uses the formula, and its gates and saturation summary look at p with p50 as the ceiling", {
  par <- c(clim = 1.5, fuel = 3, youngAge = -1, activationIntercept = a)
  args <- list(link = "powerLawSmooth", powerLawPc50 = p50, powerLawGamma = gam, powerLawShape = shp, lowerAsymptote = f)
  for (yr in c("2004", "2005")) {
    ## objFunInner() gets what .objfunSpreadFit() hands it: p50 first, the intercept on its own
    got <- do.call(spreadProbSeenC, c(list(c(maxAsymptote = p50, par), yr, climateRef = ref), args,
                                      list(activationIntercept = a)))
    ## the year's p as fireSense_spreadPredict computes it (spreadProbFromCovs() on the same design matrix)
    p <- as.numeric(spreadProbFromCovs(matOf(yr), lp3, b, f, link = "powerLawSmooth", climateRef = ref, powerLawPc50 = p50,
                                       powerLawGamma = gam, powerLawShape = shp, lowerAsymptote = f, activationIntercept = a))
    ign <- climFix$historicalFires[[yr]]$cells   # the ignition cell is 1
    expect_true(got$reached)
    expect_equal(got$sp[-ign], p[-ign])
    expect_identical(got$pSummary, spreadProbSummary(p, ceiling = p50))
  }
  ## spreadProbGates() agrees: the same p, the same ceiling
  g <- do.call(spreadProbGates, c(list(par = par, annualDTx1000 = climFix$annualDTx1000,
                                       nonAnnualDTx1000 = climFix$nonAnnualDTx1000, historicalFires = climFix$historicalFires,
                                       formulaToFit = formC, climateRef = ref, years = "2005", mutuallyExclusive = NULL,
                                       lanscape1stQuantileThresh = 0.2, maxFireSpread = 0.28), args))
  p5 <- as.numeric(spreadProbFromCovs(matOf("2005"), lp3, b, f, link = "powerLawSmooth", climateRef = ref, powerLawPc50 = p50,
                                      powerLawGamma = gam, powerLawShape = shp, lowerAsymptote = f, activationIntercept = a))
  expect_identical(g$pass, spreadProbGateTest(p5, p50, f, 0.28, 0.2)$pass)
})

test_that("on the whole landscape spreadCpp() gets per-cell pc from fuel, and burns as the full probability vector", {
  skip_if_not_installed("SpaDES.tools")
  args <- list(clim = function(pix) 0.5 + 0.4 * cos(pix), link = "powerLawSmooth", powerLawGamma = 1.9, powerLawShape = 0.5,
               lowerAsymptote = 0.05, climateRef = c(clim = 0.5), powerLawPc50 = 0.5, activationIntercept = 1,
               returnSims = TRUE, returnBurned = TRUE, seed = 3)
  for (more in list(list(), list(yearSpreadSD = 0.5))) {
    full <- do.call(wlObjective, c(list(fx = closedBuffer(), annPix = seq_len(3600L)), args, more))
    lazy <- do.call(wlObjective, c(list(whole = TRUE), args, more))
    expect_identical(lazy$burned, full$burned)
    expect_identical(lazy$sims, full$sims)
  }
  ## a drier year (climate above its reference) gives bigger fires
  hot <- do.call(wlObjective, utils::modifyList(c(list(whole = TRUE), args), list(clim = function(pix) rep(1.5, length(pix)))))
  ref0 <- do.call(wlObjective, utils::modifyList(c(list(whole = TRUE), args), list(clim = function(pix) rep(0.5, length(pix)))))
  expect_gt(mean(hot$sims$sim), mean(ref0$sims$sim))
})

test_that("spreadFitValidationData() predicts p with the activation and records its intercept", {
  par <- c(clim = 1.5, fuel = 3, youngAge = -1, activationIntercept = a)
  r <- terra::rast(nrows = 20, ncols = 20, xmin = 0, xmax = 20, ymin = 0, ymax = 20, vals = 1)
  fb <- lapply(climFix$historicalFires, function(h) dtC(pixelID = 1:nPixC, buffer = as.integer(1:nPixC <= 30), ids = 1L))
  d <- spreadFitValidationData(par, r, climFix$annualDTx1000, climFix$nonAnnualDTx1000, formC,
                               climFix$historicalFires, fb, mutuallyExclusive = NULL, Nreps = 2L, link = "powerLawSmooth",
                               powerLawGamma = gam, powerLawShape = shp, lowerAsymptote = f, climateRef = ref,
                               powerLawPc50 = p50, seed = 1)
  info <- attr(d, "spreadFitValidation")
  expect_identical(info$activationIntercept, a)
  expect_identical(info$logisticPars[[1]], p50)
  expect_lte(max(d$p), p50)
  g <- plotSpreadFitResponse(d)
  hl <- Filter(function(l) inherits(l$geom, "GeomHline"), g$layers)[[1]]
  expect_equal(hl$data$yintercept[2], qlogis(p50))
})

test_that("runDEoptim() tells the fit and the re-score of the activation, and leaves the other fits' calls as they were", {
  seen <- new.env()
  testthat::local_mocked_bindings(
    clusterSetup = function(...) list(itermax = 5, trace = FALSE, strategy = 2L, NP = 40L, cluster = NULL),
    DEoptimIterative = function(fn, lower, ...) { seen$fit <- list(...)$fitActivationIntercept; list(list(member = list(pop = matrix(lower + 0.1, 1)))) },
    .package = "clusters")
  testthat::local_mocked_bindings(
    rescorePopulation = function(pop, fn, reps, cl, seed = 1L, fnArgs = list()) {
      seen$rescore <- fnArgs; data.table::data.table(member = 1L, rep = 1L, value = 1) })
  seen$omit <- list()
  testthat::local_mocked_bindings(Cache = function(FUN, ..., omitArgs = NULL, cachePath = NULL, .functionName = NULL,
                                                   .cacheExtra = NULL, useCache = NULL) {
    seen$omit[[length(seen$omit) + 1L]] <- omitArgs
    if (is.function(FUN)) FUN(...) else FUN
  }, .package = "fireSenseUtils")
  run <- function(lower) suppressMessages(runDEoptim(
    landscape = NULL, annualDTx1000 = NULL, nonAnnualDTx1000 = NULL, fireBufferedListDT = NULL,
    historicalFires = NULL, itermax = 5, trace = FALSE, strategy = 2L, cores = c("hostA", "hostB"),
    paths = list(cachePath = withr::local_tempdir()), lower = lower, upper = lower + 1,
    mutuallyExclusive = NULL, formulaToFit = "~ 0 + clim + fuel", objFunCoresInternal = 1L, covMinMax = NULL,
    maxFireSpread = 0.3, Nreps = 1L, .verbose = FALSE, climateRef = list(clim = 0.4), powerLawPc50 = p50))
  run(c(clim = 0, fuel = 0, activationIntercept = -10))
  expect_true(seen$fit)
  expect_true(seen$rescore$fitActivationIntercept)
  expect_false("fitActivationIntercept" %in% seen$omit[[1]])
  seen$omit <- list()
  run(c(maxAsymptote = 0.25, clim = 0, fuel = 0))
  expect_null(seen$fit)
  expect_false("fitActivationIntercept" %in% names(seen$rescore))
  expect_true("fitActivationIntercept" %in% seen$omit[[1]])
})

test_that("'predictor' and 'ceiling' give the values they gave before the activation existed", {
  ## computed with fireSenseUtils 0.2.3.9181 (origin/development 970e7b5), the last version without it
  m <- matOf("2005")
  pred <- spreadLinkInputs(m, lp3, b, f, "powerLawSmooth", lowerAsymptote = 0.05, powerLawShape = 0.5)
  expect_equal(pred$x[1:4], c(2.3085, 3.0975, 1.589, 2.7405), tolerance = 1e-12)
  expect_equal(unname(pred$pars), c(0.27, 1.9, 0.13, 0.5, 0.05))
  ceil <- spreadLinkInputs(m, lp3, b, f, "powerLawSmooth", ref, p50, lowerAsymptote = 0.05, powerLawShape = 0.5)
  expect_equal(unname(ceil$pars[1:4, "pc"]), c(0.285198960990328, 0.286472868874645, 0.289501867526665, 0.286594066845817),
               tolerance = 1e-12)
  expect_equal(ceil$x[1:4], c(1.152, 1.851, 0.089, 1.485), tolerance = 1e-12)
})
