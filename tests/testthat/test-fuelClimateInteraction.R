## fuelClimateInteractionScreen(): a post-fit screen of fuel x climate interactions in the fit's own fires,
## on toy fires whose size is planted to depend on the interaction, or not.

## 30 fire years (1986-2015), 3 fires a year of 10 pixels each; the vegetation tables are of data years 1985
## and 2000. `interaction` is the planted coefficient of standardised climate x fuel in log size.
toyInteractionFit <- function(interaction, seed = 11L, nfGroups = FALSE) {
  withr::local_seed(seed)
  years <- 1986:2015
  nPer <- 3L
  nPix <- length(years) * nPer * 10L
  fuelLevel <- rnorm(length(years) * nPer)                    # each fire burns its own block of pixels
  pix <- data.frame(pixelID = seq_len(nPix), dom = rep(fuelLevel, each = 10L) * 1000,
                    youngAge = rbinom(nPix, 1, 0.3) * 1000)
  if (nfGroups) {
    pix$nfLCC_1 <- rbinom(nPix, 1, 0.3) * 500
    pix$nfLCC_2 <- rbinom(nPix, 1, 0.2) * 250
  }
  nonAnnual <- list(`1985` = pix, `2000` = pix)
  climate <- rnorm(length(years))
  annual <- lapply(seq_along(years), function(i)
    data.frame(pixelID = seq_len(nPix), CMD = round(climate[i] * 1000 + rnorm(nPix, 0, 50))))
  names(annual) <- paste0("year", years)
  fireYear <- rep(seq_along(years), each = nPer)
  logSize <- 6 + 0.3 * scale(climate)[fireYear] + 0.3 * scale(fuelLevel) +
    interaction * scale(climate)[fireYear] * scale(fuelLevel) + rnorm(length(fireYear), 0, 0.2)
  obs <- round(exp(logSize))
  buffered <- lapply(seq_along(years), function(i) {
    k <- which(fireYear == i)
    data.frame(pixelID = unlist(lapply(k, function(f) ((f - 1L) * 10L + 1L):(f * 10L))),
               buffer = 1L, ids = rep(seq_along(k), each = 10L))
  })
  names(buffered) <- names(annual)
  sims <- data.table::rbindlist(lapply(seq_along(years), function(i) {
    k <- which(fireYear == i)
    data.table::data.table(member = 1L, yr = names(annual)[i], rep = rep(1:2, each = length(k)),
                           ids = rep(seq_along(k), 2L), initialLocus = rep((k - 1L) * 10L + 1L, 2L),
                           size = rep(obs[k], 2L), sim = rep(obs[k], 2L))
  }))
  list(sims = sims, annualDTx1000 = annual, nonAnnualDTx1000 = nonAnnual, fireBufferedListDT = buffered,
       formulaToFit = if (nfGroups) "~ CMD + dom + youngAge + nfLCC_1 + nfLCC_2" else "~ CMD + dom + youngAge")
}

runInteractionScreen <- function(fit, ...) do.call(fuelClimateInteractionScreen, c(fit, list(...)))

test_that("a planted fuel x climate interaction is found: positive coefficient, leave-one-year-out improves", {
  r <- runInteractionScreen(toyInteractionFit(0.6), nBoot = 199L)
  expect_identical(r$climate, "CMD")
  expect_identical(r$nYears, 30L)
  x <- r$terms[r$terms$fuel == "dom" & r$terms$at == "burn" & r$terms$response == "logSize", ]
  expect_identical(nrow(x), 1L)
  expect_gt(x$bIntMedian, 0.3)
  expect_gt(x$bIntLo, 0)
  expect_gt(x$pIntPositive, 0.95)
  expect_lt(x$mseInteraction, x$mseAdditive)
  expect_gt(x$yearsImproved, x$yearsTested / 2)
})

test_that("without a planted interaction, none is found", {
  r <- runInteractionScreen(toyInteractionFit(0), nBoot = 199L)
  x <- r$terms[r$terms$fuel == "dom" & r$terms$at == "burn" & r$terms$response == "logSize", ]
  expect_lt(x$bIntLo, 0)
  expect_gt(x$bIntHi, 0)
  expect_gt(x$mseInteraction, x$mseAdditive)
  expect_lte(x$yearsImproved, x$yearsTested / 2)
})

test_that("the model side is tested too: the simulated size, and the residual obs / sim", {
  fit <- toyInteractionFit(0.6)
  r <- runInteractionScreen(fit, nBoot = 99L)
  t <- r$terms[r$terms$fuel == "dom" & r$terms$at == "burn", ]
  expect_setequal(t$response, c("logSize", "logSim", "resid"))
  expect_gt(t$bIntMedian[t$response == "logSim"], 0.3)         # the toy simulation copies the observed size
  expect_equal(t$bIntMedian[t$response == "resid"], 0, tolerance = 1e-8)
  ## the Spearman strata are for the observed size only
  expect_false(anyNA(t$rhoDiffMedian[t$response == "logSize"]))
  expect_true(all(is.na(t$rhoDiffMedian[t$response != "logSize"])))
})

test_that("a term is reported per fuel covariate, at the burned pixels and at the start pixel, plus the non-forest share", {
  r <- runInteractionScreen(toyInteractionFit(0.6, nfGroups = TRUE), nBoot = 49L)
  expect_setequal(r$terms$fuel, c("dom", "youngAge", "nfLCC_1", "nfLCC_2", "nonForest"))
  expect_setequal(r$terms$at, c("burn", "start"))
  expect_setequal(names(r$terms), c("fuel", "at", "response", "nFires", "nYears", "bClimMedian", "bFuelMedian",
    "bIntMedian", "bIntLo", "bIntHi", "pIntPositive", "mseAdditive", "mseInteraction", "yearsImproved",
    "yearsTested", "rhoHotMedian", "rhoRestMedian", "nHot", "nRest", "rhoDiffMedian", "rhoDiffLo",
    "rhoDiffHi", "pHotAboveRest"))
})

test_that("the result is the same for the same seed, differs for another, and leaves the session's random numbers alone", {
  fit <- toyInteractionFit(0.6)
  set.seed(5); before <- runif(1); set.seed(5)
  a <- runInteractionScreen(fit, nBoot = 49L, seed = 3L)
  expect_identical(runif(1), before)
  expect_identical(a, runInteractionScreen(fit, nBoot = 49L, seed = 3L))
  expect_false(identical(a$terms$bIntLo, runInteractionScreen(fit, nBoot = 49L, seed = 4L)$terms$bIntLo))
})

test_that("a fire year reads the vegetation table of the data year strictly before it", {
  fit <- toyInteractionFit(0)
  fit$nonAnnualDTx1000[["2000"]]$dom <- fit$nonAnnualDTx1000[["2000"]]$dom + 5000   # 5 more than the 1985 table
  f <- do.call(fuelClimateFires, c(fit[c("sims", "annualDTx1000", "nonAnnualDTx1000", "fireBufferedListDT",
                                           "formulaToFit")]))
  expect_equal(f$burn_dom, f$start_dom)                                  # a fire's pixels share one fuel level
  ## fires of 1986-2000 read the 1985 table, of 2001-2015 the 2000 table, 5 higher
  level <- fit$nonAnnualDTx1000[["1985"]]$dom[(seq_len(nrow(f)) - 1L) * 10L + 1L] / 1000   # fires are in year order
  expect_equal(f$burn_dom[f$year <= 2000], level[f$year <= 2000])
  expect_equal(f$burn_dom[f$year > 2000], level[f$year > 2000] + 5)
})

test_that("the climate is the base variable: a _yr/_dev split counts as one climate, its _yr term the year's value", {
  fit <- toyInteractionFit(0.6)
  plain <- do.call(fuelClimateFires, fit[c("sims", "annualDTx1000", "nonAnnualDTx1000", "fireBufferedListDT", "formulaToFit")])
  split <- fit
  set.seed(3)
  split$annualDTx1000 <- lapply(fit$annualDTx1000, function(a) {
    yr <- stats::median(a$CMD); dev <- rnorm(nrow(a), 0, 100)
    data.frame(pixelID = a$pixelID, CMD_yr = yr, CMD_dev = round(dev))
  })
  split$formulaToFit <- "~ CMD_dev + CMD_yr + dom + youngAge"           # either term first: the same climate
  f <- do.call(fuelClimateFires, split[c("sims", "annualDTx1000", "nonAnnualDTx1000", "fireBufferedListDT", "formulaToFit")])
  expect_identical(attr(f, "climate"), "CMD")
  expect_equal(f$climYear, plain$climYear)                              # the ELF median is the same number
  want <- vapply(seq_len(nrow(f)), function(i) {                           # the burned pixels' median of _yr + _dev
    a <- split$annualDTx1000[[f$yr[i]]]; b <- split$fireBufferedListDT[[f$yr[i]]]
    px <- b$pixelID[b$ids == f$ids[i]]
    stats::median(a$CMD_yr[px] + a$CMD_dev[px]) / 1000
  }, numeric(1))
  expect_equal(f$climBurn, want)
  r <- runInteractionScreen(split, nBoot = 49L)
  expect_identical(r$climate, "CMD")
  expect_gt(r$terms$bIntMedian[r$terms$fuel == "dom" & r$terms$at == "burn" & r$terms$response == "logSize"], 0)
})

test_that("a name with a colon (an interaction term of the formula) is neither a climate nor a fuel", {
  fit <- toyInteractionFit(0.6)
  with <- modifyList(fit, list(formulaToFit = "~ CMD + dom + youngAge + CMD:dom + dom:youngAge"))
  a <- runInteractionScreen(fit, nBoot = 49L)
  b <- runInteractionScreen(with, nBoot = 49L)
  expect_identical(b, a)
  expect_identical(attr(do.call(fuelClimateFires, with[names(fit)]), "climate"), "CMD")
  expect_null(runInteractionScreen(modifyList(fit, list(formulaToFit = "~ CMD + CMD:dom")), nBoot = 9L))
})

test_that("scoreFireSizes(), scoreFireSizesByYear() and the screen share one per-fire table", {
  d <- data.table::data.table(yr = "year2001", member = 1L, rep = 1:3, ids = 7L, initialLocus = 5L, size = 10, sim = c(4, 8, 30))
  expect_equal(perFireSizes(d)$sim, 14)
  expect_equal(perFireSizes(d, stats::median)$sim, 8)
  expect_identical(perFireSizes(d)$initialLocus, 5L)
  expect_false("initialLocus" %in% names(perFireSizes(d[, !"initialLocus"])))
})

test_that("no test without a climate covariate, without a fuel covariate, with under 4 years, or with nBoot = 0", {
  fit <- toyInteractionFit(0.6)
  expect_null(runInteractionScreen(modifyList(fit, list(formulaToFit = "~ dom + youngAge")), nBoot = 9L))
  expect_null(runInteractionScreen(modifyList(fit, list(formulaToFit = "~ CMD")), nBoot = 9L))
  expect_null(runInteractionScreen(modifyList(fit, list(formulaToFit = "~ 1")), nBoot = 9L))
  few <- fit
  few$sims <- fit$sims[fit$sims$yr %in% c("year1986", "year1987", "year1988"), ]
  expect_null(runInteractionScreen(few, nBoot = 9L))
  expect_null(runInteractionScreen(fit, nBoot = 0L))
})

test_that("the ledger column name", {
  expect_identical(fireSenseUtils::spreadFitInteractionScreenTxt, "interactionScreen_spread")
})
