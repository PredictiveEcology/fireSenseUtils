## Validation of a spread fit against its data (?spreadFitValidation): the objective records which
## pixels burn in its simulations (returnBurned), and the observed and simulated burned shares are
## compared per covariate. A synthetic 40 x 40 landscape, two fire years, the real spreadCpp().

toyValidationInputs <- function() {
  dt <- data.table::data.table
  landscape <- terra::rast(nrows = 40, ncols = 40, xmin = 0, xmax = 40, ymin = 0, ymax = 40, vals = 1)
  cellOf <- function(row, col) (row - 1L) * 40L + col
  block <- function(r0, c0, n) as.integer(outer(seq(r0, r0 + n - 1L), seq(c0, c0 + n - 1L), cellOf))
  ## each fire: a 15 x 15 buffer, of which the central 5 x 5 burned; the ignition at the centre
  fire <- function(r0, c0, id) {
    px <- block(r0, c0, 15L)
    dt(pixelID = px, buffer = as.integer(px %in% block(r0 + 5L, c0 + 5L, 5L)), ids = id)
  }
  fb <- list(year2001 = rbind(fire(1L, 1L, 1L), fire(20L, 20L, 2L)),
             year2002 = fire(10L, 5L, 3L))
  hist <- list(year2001 = data.frame(size = c(25L, 25L), cells = c(cellOf(8L, 8L), cellOf(27L, 27L)),
                                     ids = 1:2, date = "year2001"),
               year2002 = data.frame(size = 25L, cells = cellOf(17L, 12L), ids = 3L, date = "year2002"))
  set.seed(10)
  allPx <- sort(unique(unlist(lapply(fb, `[[`, "pixelID"))))
  nonAnnual <- list(`2001` = dt(   # named by first year, as fireSense_spreadFit names them
    pixelID = allPx,
    agb = as.integer(round(ifelse(stats::runif(length(allPx)) < 0.3, 0, stats::runif(length(allPx), 0, 5000)) * 1000)),
    youngAge = as.integer(stats::runif(length(allPx)) < 0.1) * 1000L))
  annual <- lapply(fb, function(b) dt(pixelID = b$pixelID, clim = as.integer(stats::runif(nrow(b), 0, 100) * 1000)))
  list(landscape = landscape, annualDTx1000 = annual, nonAnnualDTx1000 = nonAnnual,
       historicalFires = hist, fireBufferedListDT = fb,
       formulaToFit = "~ 0 + clim + youngAge + agb",
       covMinMax = data.table::data.table(clim = c(0, 100), youngAge = c(0, 1), agb = c(0, 1e4)),
       mutuallyExclusive = list(youngAge = "agb"),
       par = c(maxAsymptote = 0.26, clim = 2, youngAge = -2, agb = 3))
}

callObjective <- function(inp, ...) {
  fireSenseUtils::.objfunSpreadFit(
    par = inp$par, landscape = inp$landscape, annualDTx1000 = inp$annualDTx1000,
    nonAnnualDTx1000 = inp$nonAnnualDTx1000, formulaToFit = inp$formulaToFit,
    historicalFires = inp$historicalFires, fireBufferedListDT = inp$fireBufferedListDT,
    covMinMax = inp$covMinMax, mutuallyExclusive = inp$mutuallyExclusive,
    tests = c("snll_fs", "adTest"), Nreps = 5L, doAssertions = FALSE, verbose = 0, weighted = FALSE, ...)
}

test_that("the objective's value is identical with recording off, and recording draws no random numbers", {
  skip_if_not_installed("SpaDES.tools")
  inp <- toyValidationInputs()
  set.seed(123); v0 <- callObjective(inp)
  set.seed(123); v1 <- callObjective(inp, returnBurned = FALSE)
  expect_true(is.finite(v0))
  expect_lt(v0, 1e6)                                     # a real value, not the fail value
  expect_identical(v1, v0)
  expect_false(formals(fireSenseUtils::.objfunSpreadFit)$returnBurned)
  ## the simulations are the same draws whether or not the pixels are recorded
  set.seed(5); s0 <- callObjective(inp, returnSims = TRUE)
  set.seed(5); s1 <- callObjective(inp, returnBurned = TRUE)
  expect_null(attr(s0, "burned"))
  expect_false(is.null(attr(s1, "burned")))
  data.table::setattr(s1, "burned", NULL)
  expect_identical(s1, s0)
})

test_that("returnBurned records, per year and replicate, pixels inside the buffers matching the sizes", {
  skip_if_not_installed("SpaDES.tools")
  inp <- toyValidationInputs()
  set.seed(7)
  s <- callObjective(inp, returnBurned = TRUE)
  b <- attr(s, "burned")
  expect_setequal(names(b), c("year2001", "year2002"))
  for (y in names(b)) {
    expect_named(b[[y]], c("rep", "pixelID"))
    expect_true(all(b[[y]]$pixelID %in% inp$fireBufferedListDT[[y]]$pixelID))
    expect_setequal(unique(b[[y]]$rep), 1:5)
    expect_false(anyDuplicated(b[[y]][, c("rep", "pixelID")]) > 0)   # one fire per cell per replicate
    ## the pixels burned in a replicate add up to that replicate's simulated fire sizes
    sizes <- s[s$yr == y, list(total = sum(sim)), by = "rep"]
    counted <- b[[y]][, list(n = .N), by = "rep"]
    expect_equal(counted$n[match(sizes$rep, counted$rep)], sizes$total)
  }
})

test_that("spreadFitValidationData() gives one row per spreadable buffer pixel-year, with shares and raw values", {
  skip_if_not_installed("SpaDES.tools")
  inp <- toyValidationInputs()
  d <- spreadFitValidationData(
    par = inp$par, landscape = inp$landscape, annualDTx1000 = inp$annualDTx1000,
    nonAnnualDTx1000 = inp$nonAnnualDTx1000, formulaToFit = inp$formulaToFit,
    historicalFires = inp$historicalFires, fireBufferedListDT = inp$fireBufferedListDT,
    covMinMax = inp$covMinMax, mutuallyExclusive = inp$mutuallyExclusive, Nreps = 4L, seed = 1L)
  expect_identical(names(d), c("year", "pixelID", "observed", "simulated", "p", "clim", "youngAge", "agb"))
  expect_identical(nrow(d), sum(vapply(inp$fireBufferedListDT, nrow, integer(1))))
  expect_identical(sum(d$observed), 75L)                          # three fires of 25 burned pixels
  expect_true(all(d$simulated >= 0 & d$simulated <= 1))
  expect_true(all(d$simulated * 4 == round(d$simulated * 4)))    # a count of 4 replicates
  expect_gt(sum(d$simulated), 0)
  expect_true(all(d$p >= 0.13 & d$p <= 0.26))
  ## raw values: the stored integers / 1000, before mutual exclusivity (youngAge zeroes agb)
  y1 <- d[d$year == "year2001"]
  na <- inp$nonAnnualDTx1000[[1]]
  expect_equal(y1$agb, na$agb[match(y1$pixelID, na$pixelID)] / 1000)
  expect_true(any(y1$youngAge == 1 & y1$agb > 0))
  ## p as the objective computes it: the rescaled, mutually exclusive covariates through the link
  u <- cbind(clim = y1$clim / 100, youngAge = y1$youngAge, agb = ifelse(y1$youngAge == 1, 0, y1$agb / 1e4))
  expect_equal(y1$p, as.numeric(0.13 + (0.26 - 0.13) / (1 + exp(-(u %*% c(2, -2, 3))))^1))
  info <- attr(d, "spreadFitValidation")
  expect_identical(info$climateCols, "clim")
  expect_identical(info$yearsNotSimulated, character(0))
  expect_equal(info$Nreps, 4L)
  ## the seed makes it repeatable, and the session's random numbers are left as they were
  set.seed(99); before <- .Random.seed
  d2 <- spreadFitValidationData(
    par = inp$par, landscape = inp$landscape, annualDTx1000 = inp$annualDTx1000,
    nonAnnualDTx1000 = inp$nonAnnualDTx1000, formulaToFit = inp$formulaToFit,
    historicalFires = inp$historicalFires, fireBufferedListDT = inp$fireBufferedListDT,
    covMinMax = inp$covMinMax, mutuallyExclusive = inp$mutuallyExclusive, Nreps = 4L, seed = 1L)
  expect_identical(.Random.seed, before)
  expect_equal(d2, d)
})

## a toy table of the shape spreadFitValidationData() returns, for the plots
toyValidationTable <- function() {
  set.seed(3)
  n <- 2000
  d <- data.table::data.table(year = rep(c("year2001", "year2002"), each = n / 2), pixelID = seq_len(n),
                              observed = rbinom(n, 1, 0.1), simulated = runif(n, 0, 0.3),
                              p = runif(n, 0.14, 0.25), clim = runif(n, 0, 100),
                              youngAge = rbinom(n, 1, 0.1),
                              agb = ifelse(runif(n) < 0.4, 0, runif(n, 0, 5000)))
  data.table::setattr(d, "spreadFitValidation", list(
    logisticPars = c(maxAsymptote = 0.26, hillSlope1 = 1, inflectionPoint1 = 1),
    covPars = c(clim = 2, youngAge = -2, agb = 3),
    covMinMax = list(clim = c(0, 100), youngAge = c(0, 1), agb = c(0, 1e4)), covCentre = NULL,
    lowerSpreadProb = 0.13, link = NULL, Nreps = 4L, covariates = c("clim", "youngAge", "agb"),
    climateCols = "clim", yearsNotSimulated = character(0)))
  d
}

test_that("plotSpreadFitValidation() returns a ggplot of observed and simulated shares per covariate", {
  d <- toyValidationTable()
  g <- plotSpreadFitValidation(d, nBins = 5L)
  expect_s3_class(g, "ggplot")
  expect_identical(levels(g$data$series), c("observed", "simulated by the fit"))
  bins <- g$data[g$data$series == "observed"]
  expect_identical(as.character(unique(bins$covariate)), c("clim", "youngAge", "agb"))
  expect_identical(nrow(bins[bins$covariate == "youngAge"]), 2L)   # 0/1: one bin per value
  expect_identical(nrow(bins[bins$covariate == "clim"]), 5L)
  ## fuel: only where it is present
  expect_gt(min(bins[bins$covariate == "agb"]$value), 0)
  expect_equal(sum(g$layers[[3]]$data$n[g$layers[[3]]$data$covariate == "agb"]), sum(d$agb > 0))
  expect_no_error(ggplot2::ggplot_build(g))
  expect_error(plotSpreadFitValidation(data.table::data.table(a = 1)), "spreadFitValidationData")
})

test_that("plotSpreadFitResponse() returns a ggplot of the pure-stand curves, titled as the model's response", {
  d <- toyValidationTable()
  g <- plotSpreadFitResponse(d)
  expect_s3_class(g, "ggplot")
  expect_match(g$labels$title, "not how well it fits")
  cur <- g$layers[[3]]$data
  expect_identical(as.character(unique(cur$covariate)), c("clim", "agb"))   # 0/1 covariates have no curve
  ## agb's curve: agb alone from 0, climate at its median, youngAge 0
  a <- cur[cur$covariate == "agb"]
  expect_equal(a$value[1], 0)
  eta <- 2 * stats::median(d$clim) / 100 + 3 * a$value / 1e4
  expect_equal(a$logitP, stats::qlogis(0.13 + 0.13 / (1 + exp(-eta))))
  expect_no_error(ggplot2::ggplot_build(g))
})

## With an intercept (formula "~ 1 + ..."): the intercept is a coefficient, not a covariate to plot.
test_that("the validation data omits the intercept from its covariates but keeps its coefficient", {
  skip_if_not_installed("SpaDES.tools")
  inp <- toyValidationInputs()
  inp$formulaToFit <- "~ 1 + clim + youngAge + agb"
  inp$par <- c(maxAsymptote = 0.26, `(Intercept)` = -0.4, clim = 2, youngAge = -2, agb = 3)
  set.seed(3)
  d <- spreadFitValidationData(
    par = inp$par, landscape = inp$landscape, annualDTx1000 = inp$annualDTx1000,
    nonAnnualDTx1000 = inp$nonAnnualDTx1000, formulaToFit = inp$formulaToFit,
    historicalFires = inp$historicalFires, fireBufferedListDT = inp$fireBufferedListDT,
    covMinMax = inp$covMinMax, mutuallyExclusive = inp$mutuallyExclusive,
    covCentre = list(clim = 0.5, youngAge = 0.1, agb = 0.2), Nreps = 2L)
  info <- attr(d, "spreadFitValidation")
  expect_identical(info$covariates, c("clim", "youngAge", "agb"))
  expect_true(info$intercept)
  expect_false(spreadInterceptTxt %in% names(d))
  expect_named(info$covPars, c(spreadInterceptTxt, "clim", "youngAge", "agb"))
  ## p is the hand-built logistic of the intercept and the centred, rescaled covariates
  raw <- d
  x <- list(clim = raw$clim / 100, youngAge = raw$youngAge, agb = raw$agb / 1e4)
  # (mutual exclusivity zeroes agb where youngAge is 1, as the objective does)
  x$agb[x$youngAge > 0.5] <- 0
  lp <- -0.4 + 2 * (x$clim - 0.5) - 2 * (x$youngAge - 0.1) + 3 * (x$agb - 0.2)
  expect_equal(d$p, 0.13 + (0.26 - 0.13) * plogis(lp), tolerance = 1e-6)
  ## the response curves take the intercept into account without a column for it
  expect_s3_class(plotSpreadFitResponse(d), "ggplot")
})
