## How much of the fire signal is climate, how much fuel: spreadVariancePartition() (model-based, one fit) and
## commonalityPartition() (from fits with and without each group of terms). Toy data where the answer is known.

## Five years and 60 pixels, as fireSense_dataPrepFit gives them (real scale, fuel as logMinB()). `climYear`: climate is one
## value per year (it varies only between years); `climPixel`: it also varies between the pixels of a year. Fuel and youngAge
## vary between pixels and are the same every year.
vpTables <- function(climPixel = FALSE, nPix = 60L, years = 2001:2005) {
  set.seed(3)
  fuel <- logMinB(runif(nPix, 0, 8000))
  young <- as.integer(seq_len(nPix) %% 6 == 0)
  annual <- lapply(seq_along(years), function(i)
    data.table::data.table(pixelID = seq_len(nPix), CMD = 10 * i + if (climPixel) runif(nPix, 0, 30) else 0))
  names(annual) <- years
  list(annual = annual,
       nonAnnual = list(`2000` = data.table::data.table(pixelID = seq_len(nPix), dom_agb_Pice_mar = fuel, youngAge = young)))
}
vpCovMinMax <- data.table::data.table(CMD = c(0, 100), dom_agb_Pice_mar = fuelLinearRange, youngAge = c(0, 1))
vpRow <- function(pars, centre = NULL) {
  row <- data.frame(polygonID = "4.1")
  row$params <- list(data.frame(maxAsymptote = 0.25, hillSlope1 = 1, inflectionPoint1 = 1, pars, yearSpreadSD = 1,
                                check.names = FALSE))
  row$covMinMax_spread <- list(vpCovMinMax)
  if (!is.null(centre)) row[[spreadFitCovCentreTxt]] <- list(centre)
  row
}
shareOf <- function(vp, signal, part) vp$share[vp$signal == signal & vp$part == part]

test_that("spreadTermGroups() puts each design column in its group; an interaction with climate is climate", {
  g <- spreadTermGroups(c("(Intercept)", "CMD", "dom_agb_Pice_mar", "sec_agb_Pinu_ban", "treedWetland_agb", "youngAge",
                          "nfLCC_1", "CMD:dom_agb_Pice_mar", "dom_agb_Pice_mar:youngAge", "treedWetland"),
                        climateCols = "CMD")
  expect_identical(unname(g[c("CMD", "dom_agb_Pice_mar", "sec_agb_Pinu_ban", "treedWetland_agb", "youngAge", "nfLCC_1",
                              "CMD:dom_agb_Pice_mar", "dom_agb_Pice_mar:youngAge", "treedWetland")]),
                   c("climate", "fuel", "fuel", "fuel", "youngAge", "nonForest", "climate", "other", "other"))
  expect_false("(Intercept)" %in% names(g))           # a constant: no variance to share
})

test_that("spreadDesignByYear() builds the fit's design: linear fuel, youngAge clears fuel, centred as the row says", {
  tb <- vpTables(climPixel = TRUE)
  row <- vpRow(list(`(Intercept)` = 0, CMD = 1, dom_agb_Pice_mar = 2, youngAge = -1))
  d <- spreadDesignByYear(row, tb$annual, tb$nonAnnual)
  expect_identical(attr(d, "climateCols"), "CMD")
  expect_identical(sort(unique(d$year)), 2001:2005)
  expect_equal(nrow(d), 5L * 60L)
  fuelLinear <- fuelLogToLinear(tb$nonAnnual[[1]]$dom_agb_Pice_mar) / 1e4
  young <- tb$nonAnnual[[1]]$youngAge == 1
  y1 <- d[d$year == 2001L]
  expect_equal(y1$dom_agb_Pice_mar, ifelse(young, 0, fuelLinear), tolerance = 1e-3)  # x1000 integers, as the fit
  expect_equal(y1$CMD, tb$annual[["2001"]]$CMD / 100, tolerance = 1e-3)
  ## the input tables are not changed
  expect_equal(tb$annual[["2001"]]$CMD, vpTables(climPixel = TRUE)$annual[["2001"]]$CMD)
  ## with the row's centres, each column is shifted by them
  cen <- list(CMD = 0.3, dom_agb_Pice_mar = 0.1, youngAge = 0.2)
  dc <- spreadDesignByYear(vpRow(list(`(Intercept)` = 0, CMD = 1, dom_agb_Pice_mar = 2, youngAge = -1), centre = cen),
                           tb$annual, tb$nonAnnual)
  expect_equal(dc$CMD, d$CMD - 0.3)
  expect_equal(dc$dom_agb_Pice_mar, d$dom_agb_Pice_mar - 0.1)
})

test_that("only climate varies between years: the between-year signal is all climate, the within-year all fuel and youngAge", {
  tb <- vpTables(climPixel = FALSE)
  row <- vpRow(list(CMD = 3, dom_agb_Pice_mar = 2, youngAge = -1))
  vp <- spreadVariancePartition(row, spreadDesignByYear(row, tb$annual, tb$nonAnnual))
  expect_setequal(unique(vp$signal), c("between", "within"))
  expect_equal(shareOf(vp, "between", "climate"), 1)
  expect_equal(shareOf(vp, "between", "fuel"), 0)
  expect_equal(shareOf(vp, "within", "climate"), 0)
  expect_equal(shareOf(vp, "within", "fuel") + shareOf(vp, "within", "youngAge"), 1)
  ## each signal's shares add to 1; the two signals add to the variance of the linear predictor
  expect_equal(as.vector(tapply(vp$share, vp$signal, sum)), c(1, 1))
  sig <- vp[!duplicated(vp$signal), ]
  expect_equal(sum(sig$signalShare), 1)
})

test_that("a coefficient of 0 gets no share; the shares follow the variances of independent groups", {
  tb <- vpTables(climPixel = TRUE)
  row <- vpRow(list(CMD = 2, dom_agb_Pice_mar = 0, youngAge = 0))
  vp <- spreadVariancePartition(row, spreadDesignByYear(row, tb$annual, tb$nonAnnual))
  expect_equal(shareOf(vp, "within", "climate"), 1)
  expect_equal(shareOf(vp, "within", "fuel"), 0)
  ## lmgShares(): uncorrelated parts share by their variances; two copies of one part split it evenly
  expect_equal(unname(lmgShares(diag(c(1, 3)))), c(0.25, 0.75))
  expect_equal(unname(lmgShares(matrix(1, 2, 2))), c(0.5, 0.5))
  ## a part with no variance gets 0, the others still add to 1
  expect_equal(unname(lmgShares(diag(c(2, 0, 2)))), c(0.5, 0, 0.5))
})

test_that("own variance and covariance: correlated climate and fuel share the covariance", {
  tb <- vpTables(climPixel = TRUE)
  ## fuel tracks the pixel's climate within each year: make the fuel table a function of the 2001 climate
  tb$nonAnnual[[1]]$dom_agb_Pice_mar <- logMinB(200 * tb$annual[["2001"]]$CMD)
  for (y in names(tb$annual)) tb$annual[[y]]$CMD <- tb$annual[["2001"]]$CMD + 10 * (as.integer(y) - 2001L)
  row <- vpRow(list(CMD = 1, dom_agb_Pice_mar = 1, youngAge = 0))
  vp <- spreadVariancePartition(row, spreadDesignByYear(row, tb$annual, tb$nonAnnual))
  w <- vp[vp$signal == "within", ]
  expect_gt(sum(w$covariance), 0.1)
  expect_equal(sum(w$own) + sum(w$covariance), 1)
  expect_equal(sum(w$share), 1)
})

test_that("spreadVariancePartition() with several parameter sets gives their median, and byParamSet each", {
  tb <- vpTables(climPixel = TRUE)
  row <- vpRow(list(CMD = c(1, 2, 3), dom_agb_Pice_mar = 2, youngAge = -1))
  d <- spreadDesignByYear(row, tb$annual, tb$nonAnnual)
  each <- spreadVariancePartition(row, d, byParamSet = TRUE)
  expect_setequal(unique(each$paramSet), 1:3)
  med <- spreadVariancePartition(row, d)
  expect_equal(shareOf(med, "within", "climate"),
               stats::median(each$share[each$signal == "within" & each$part == "climate"]))
  expect_error(spreadVariancePartition(row, d[, !"CMD"]), "CMD")
})

test_that("commonalityPartition() splits the explained share into unique climate, unique fuel, shared and unexplained", {
  scores <- data.frame(signal = "between", climate = c(FALSE, TRUE, FALSE, TRUE), fuel = c(FALSE, FALSE, TRUE, TRUE),
                       explained = c(0.1, 0.5, 0.3, 0.6))
  cp <- commonalityPartition(scores)
  s <- stats::setNames(cp$share, cp$part)
  expect_equal(s[c("baseline", "climate", "fuel", "shared", "unexplained")],
               c(baseline = 0.1, climate = 0.3, fuel = 0.1, shared = 0.1, unexplained = 0.4))
  expect_equal(sum(cp$share), 1)
  expect_error(commonalityPartition(scores[-4, ]), "climate and fuel")
})

test_that("commonalityScores() reads each fit's in-sample annual areas and fire sizes from its ledger row", {
  yrs <- 2001:2008
  obs <- c(10, 200, 50, 3000, 20, 800, 5, 100)
  fit <- function(sim, fireRMSE) {
    row <- data.frame(polygonID = "4.1", objFunVal = 1)
    row[[spreadFitInSampleSizesTxt]] <- list(list(sizes = data.table::data.table(fireRMSE = fireRMSE),
                                                  byYear = data.table::data.table(yr = yrs, obs = obs, simMedian = sim)))
    row
  }
  rows <- list(none = fit(rep(100, 8), 1), climate = fit(obs * 2, 0.8), fuel = fit(rev(obs), 0.5),
               climateFuel = fit(obs, 0.4))
  sc <- commonalityScores(rows)
  b <- sc[sc$signal == "between", ]
  expect_equal(b$explained[b$climate & !b$fuel], 1, tolerance = 1e-3) # twice the observed: all but the log1p offset
  expect_equal(b$explained[!b$climate & !b$fuel], 0)           # the same every year: nothing
  w <- sc[sc$signal == "within", ]
  expect_equal(w$explained[w$climate & w$fuel], 1 - 0.4^2 / 1^2)
  expect_equal(w$explained[!w$climate & !w$fuel], 0)
  cp <- commonalityPartition(sc)
  expect_equal(as.vector(tapply(cp$share, cp$signal, sum)), c(1, 1))
  expect_error(commonalityScores(rows[1:3]), "climateFuel")
})

test_that("plotSignalPartition() draws both partitions", {
  tb <- vpTables(climPixel = TRUE)
  row <- vpRow(list(CMD = 1, dom_agb_Pice_mar = 2, youngAge = -1))
  vp <- spreadVariancePartition(row, spreadDesignByYear(row, tb$annual, tb$nonAnnual))
  expect_s3_class(plotSignalPartition(vp, title = "4.1"), "ggplot")
})
