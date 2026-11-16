## Fuel x climate interactions in the spread model: a formula term "clim:fuel" is a column of the design, the
## product of its two covariates on the 0-1 scale (after rescaling and the mutually exclusive step, before
## centring), with its own coefficient. Fit and prediction form it in spreadProbFromIntegerCovs(). Without
## one, nothing changes. Fixtures (climFix, spreadProbSeenC, ...) are in helper-spreadLinkFixture.R.

formI <- "~ 0 + clim + fuel + clim:fuel"

test_that("the interaction is a design column after the main effects; helpers find it and its sources", {
  cols <- spreadDesignCols("~ 1 + clim + fuel + clim:fuel")
  expect_identical(cols, c(spreadInterceptTxt, "clim", "fuel", "clim:fuel"))
  expect_identical(spreadInteractionCols(cols), "clim:fuel")
  expect_identical(spreadCovCols(cols), c("clim", "fuel", "clim:fuel"))
  expect_identical(spreadSourceCols(cols), c("clim", "fuel"))
  ## the covariates of an interaction are sources even without a main effect
  expect_identical(spreadSourceCols(spreadDesignCols("~ 0 + clim + clim:fuel")), c("clim", "fuel"))
  ## without interactions nothing is found and the sources are the covariates
  expect_identical(spreadInteractionCols(c("a", "b")), character(0))
  expect_identical(spreadSourceCols(c(spreadInterceptTxt, "a", "b")), c("a", "b"))
  expect_identical(spreadFuelCols(c("clim", "dom_agb_Pice_mar", "nfLCC_1_2", "treedWetland", "treedWetland_agb")),
                   c("clim", "dom_agb_Pice_mar", "treedWetland_agb"))
})

test_that("the objective's spreadProb with an interaction is the hand-built linear predictor", {
  skip_if_not_installed("SpaDES.tools")
  par <- c(maxAsymptote = 0.27, clim = 2, fuel = 3, `clim:fuel` = -4)
  got <- spreadProbSeenC(par, "2005", form = formI)$sp
  m <- matOf("2005", c("clim", "fuel"))
  expected <- f + (0.27 - f) * plogis(2 * m[, "clim"] + 3 * m[, "fuel"] - 4 * m[, "clim"] * m[, "fuel"])
  expected[2L] <- 1 # the ignition cell of 2005
  expect_equal(got[1:nPixC], unname(expected))
  ## the interaction coefficient is used: without it the probability is another one
  par0 <- par; par0[["clim:fuel"]] <- 0
  expect_false(isTRUE(all.equal(got, spreadProbSeenC(par0, "2005", form = formI)$sp)))
})

test_that("a formula without interactions gives the columns and the probabilities it gave", {
  skip_if_not_installed("SpaDES.tools")
  par <- c(maxAsymptote = 0.27, clim = 2, fuel = 3, youngAge = -1)
  got <- spreadProbSeenC(par, "2005")$sp
  m <- matOf("2005")
  expected <- f + (0.27 - f) * plogis(2 * m[, "clim"] + 3 * m[, "fuel"] - m[, "youngAge"])
  expected[2L] <- 1
  expect_equal(got[1:nPixC], unname(expected))
})

test_that("youngAge clears the fuel, and so its interaction; centring then treats it like every covariate", {
  ann <- dtC(pixelID = 1:4, clim = c(100L, 400L, 700L, 1000L))
  nonAnn <- list(`2000` = dtC(pixelID = 1:4, fuel = c(200L, 600L, 800L, 1000L), youngAge = c(0L, 0L, 1000L, 1000L)))
  cols <- spreadDesignCols("~ 1 + clim + fuel + youngAge + clim:fuel")
  mex <- list(youngAge = "fuel")
  run <- function(covCentre = NULL)
    spreadProbFromIntegerCovs(NULL, data.table::copy(ann), data.table::copy(nonAnn), nonAnnualIndex(nonAnn), 2005,
                              covMinMax = NULL, mutuallyExclusive = mex, colsToUse = cols, doAssertions = TRUE,
                              logisticPars = c(0.2, 1, 1), covPars = rep(1, length(cols)), maxFireSpread = 0.28,
                              lowerSpreadProb = 0.1, covCentre = covCentre)
  x <- run()
  expect_equal(x$fuel, c(0.2, 0.6, 0, 0))
  expect_equal(x[["clim:fuel"]], c(0.1 * 0.2, 0.4 * 0.6, 0, 0)) # 0 where the fuel was cleared
  ## with a centre, the interaction is centred by its own mean, the way fuel is: cleared pixels are at -centre for both
  ctr <- list(clim = 0.5, fuel = 0.3, youngAge = 0.1, `clim:fuel` = 0.06)
  y <- run(ctr)
  expect_equal(y$fuel, x$fuel - 0.3)
  expect_equal(y[["clim:fuel"]], x[["clim:fuel"]] - 0.06)
  expect_equal(y$fuel[3] - ctr$fuel, y$fuel[4] - ctr$fuel)
  expect_identical(y[[spreadInterceptTxt]], rep(1, 4))
})

test_that("an interaction is rescaled from its covariates' ranges, not its own", {
  ann <- dtC(pixelID = 1:3, clim = c(0L, 500L, 1000L))
  nonAnn <- list(`2000` = dtC(pixelID = 1:3, fuel = c(0L, 1000L, 2000L)))
  mm <- data.table::data.table(clim = c(0, 2), fuel = c(0, 4)) # the data's range x 1000 is inside it
  x <- spreadProbFromIntegerCovs(NULL, ann, nonAnn, nonAnnualIndex(nonAnn), 2005, covMinMax = mm,
                                 mutuallyExclusive = NULL, colsToUse = c("clim", "fuel", "clim:fuel"),
                                 doAssertions = TRUE, logisticPars = c(0.2, 1, 1), covPars = rep(1, 3),
                                 maxFireSpread = 0.28, lowerSpreadProb = 0.1)
  expect_equal(x$clim, c(0, 0.25, 0.5))
  expect_equal(x$fuel, c(0, 0.25, 0.5))
  expect_equal(x[["clim:fuel"]], c(0, 0.0625, 0.25))
})

test_that("spreadCovCentre() centres the interaction on the mean of its product", {
  annualDTx1000 <- climFix$annualDTx1000
  ctr <- spreadCovCentre(annualDTx1000, climFix$nonAnnualDTx1000, "~ 1 + clim + fuel + clim:fuel")
  expect_named(ctr, c("clim", "fuel", "clim:fuel"))
  prod <- unlist(lapply(annualDTx1000, function(a) {
    d <- climFix$nonAnnualDTx1000[[1]][a, on = "pixelID"]
    d$clim / 1000 * d$fuel / 1000
  }))
  expect_equal(ctr[["clim:fuel"]], mean(prod))
  ## without an intercept there is no centre, with or without the interaction
  expect_null(spreadCovCentre(annualDTx1000, climFix$nonAnnualDTx1000, formI))
})

test_that("an interaction with a climate covariate is not itself a climate covariate", {
  expect_identical(spreadClimateCols("~ 0 + clim + fuel + clim:fuel", climFix$annualDTx1000), "clim")
})

## The validation data and the response curves take the interaction into account.
test_that("spreadFitValidationData() has the interaction's probability; the response curve of a covariate includes it", {
  skip_if_not_installed("SpaDES.tools")
  inp <- toyValidationInputs()
  inp$formulaToFit <- "~ 1 + clim + youngAge + agb + clim:agb"
  inp$par <- c(maxAsymptote = 0.26, `(Intercept)` = -0.4, clim = 2, youngAge = -2, agb = 3, `clim:agb` = -4)
  centre <- list(clim = 0.5, youngAge = 0.1, agb = 0.2, `clim:agb` = 0.1)
  set.seed(3)
  d <- spreadFitValidationData(
    par = inp$par, landscape = inp$landscape, annualDTx1000 = inp$annualDTx1000,
    nonAnnualDTx1000 = inp$nonAnnualDTx1000, formulaToFit = inp$formulaToFit,
    historicalFires = inp$historicalFires, fireBufferedListDT = inp$fireBufferedListDT,
    covMinMax = inp$covMinMax, mutuallyExclusive = inp$mutuallyExclusive, covCentre = centre, Nreps = 2L)
  ## by hand: the covariates are in the table in raw units; the model rescales them with covMinMax, and youngAge clears agb
  clim <- d$clim / 100; agb <- ifelse(d$youngAge > 0, 0, d$agb / 1e4)
  lp <- -0.4 + 2 * (clim - 0.5) - 2 * (d$youngAge - 0.1) + 3 * (agb - 0.2) - 4 * (clim * agb - 0.1)
  expect_equal(d$p, 0.13 + (0.26 - 0.13) * plogis(lp), tolerance = 1e-6)
  ## its curves are those of the covariates, with the interaction in them: agb alone from 0, climate at its median
  g <- plotSpreadFitResponse(d)
  cur <- g$layers[[3]]$data
  expect_identical(as.character(unique(cur$covariate)), c("clim", "agb"))
  a <- cur[cur$covariate == "agb"]
  cm <- stats::median(d$clim) / 100
  eta <- -0.4 + 2 * (cm - 0.5) - 2 * (0 - 0.1) + 3 * (a$value / 1e4 - 0.2) - 4 * (cm * a$value / 1e4 - 0.1)
  expect_equal(a$logitP, stats::qlogis(0.13 + 0.13 / (1 + exp(-eta))))
  expect_no_error(ggplot2::ggplot_build(g))
})
