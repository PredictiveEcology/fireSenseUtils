#' Does fuel act differently in hot, dry years? A post-fit test of fuel x climate interactions
#'
#' A spread model is additive in the logit of the spread probability: a fuel term has the same effect in every
#' year's climate. `fuelClimateInteractionTest()` asks the fit's own fires whether that holds. For each fitted fire it
#' takes the observed size, the year's climate and the fuel the fire burned (and the fuel at its start pixel), and
#' compares `log(size) ~ climate + fuel` with `log(size) ~ climate + fuel + climate:fuel`, one model per fuel
#' covariate, with standardised predictors. It uses only what is in memory at the end of a fit and runs no simulation.
#'
#' * `fuelClimateFires()` is the per-fire table: observed size (`obs`), the median in-sample simulated size over the
#'   fit's members and runs (`simMedian`), the year's climate (`climYear`: median of the fit's climate covariate over the
#'   year's pixels, `climBurn`: median over the fire's burned pixels), and one column per fuel, `burn_<fuel>`
#'   (mean over the burned pixels) and `start_<fuel>` (at the start pixel). The fuels are the non-annual
#'   covariates of the formula and, with two or more `nfLCC_*` groups, `nonForest` (their sum, the non-forest share).
#'   A fire year reads the vegetation table of the data year strictly before it ([dataYearBefore()]), as the fit does.
#' * `fuelClimateInteractionTest()` fits the models to the observed size (`logSize`), the simulated size
#'   (`logSim`) and the residual `log(obs / simMedian)` (`resid`).
#'
#' @param sims the in-sample simulations of the fit: `attr(DE, "fitSims")` of [runDEoptim()], the result of
#'   [simulateFireSizes()] (the columns `yr`, `ids`, `initialLocus`, `size` and `sim`).
#' @param annualDTx1000,nonAnnualDTx1000,fireBufferedListDT The covariate tables the fit used, as given to
#'   [.objfunSpreadFit()] (integers x 1000; standardising removes the scale).
#' @param formulaToFit The spread formula. Its first annual (climate) covariate is the climate; its non-annual
#'   covariates are the fuels.
#' @param nBoot Number of bootstrap draws over years. `0` skips the test.
#' @param seed Seed of the draws (the session's random number stream is restored afterwards), so the result is
#'   the same for the same fit.
#'
#' @return `fuelClimateFires()`: a `data.table`, one row per fitted fire, or `NULL` when the formula has no climate
#'   covariate or no fuel covariate. `fuelClimateInteractionTest()`: `NULL` in that case, with fewer than 4 years
#'   with fires, or with `nBoot = 0`; else a list: `terms`, one row per fuel, `at` (`"burn"` or `"start"`) and
#'   `response`, with `nFires` and `nYears`; `bClimMedian`, `bFuelMedian`, and the interaction coefficient
#'   (`bIntMedian`, `bIntLo`, `bIntHi`: the median and 2.5-97.5% range over `nBoot` bootstrap draws that resample
#'   years, and `pIntPositive`, their share above 0); `mseAdditive` and `mseInteraction`, the mean squared error
#'   of leaving each year out in turn, `yearsImproved` of `yearsTested` years where the interaction model is better;
#'   and, for `logSize` only, the Spearman correlation of size with the fuel in the hottest third of years and in the
#'   rest (`rhoHotMedian`, `rhoRestMedian`, `nHot`, `nRest`), their difference `rhoDiffMedian`, `rhoDiffLo`,
#'   `rhoDiffHi` and `pHotAboveRest`, bootstrapped over years alike. The list also has `climate` (the covariate),
#'   `nFires`, `nYears`, `nBoot` and `seed`. It is what the ledger row's [spreadFitInteractionTestTxt] column holds.
#' @export
#' @rdname fuelClimateInteraction
fuelClimateInteractionTest <- function(sims, annualDTx1000, nonAnnualDTx1000, fireBufferedListDT, formulaToFit,
                                       nBoot = 999L, seed = 1L) {
  if (!isTRUE(nBoot > 0)) return(NULL)
  fires <- fuelClimateFires(sims, annualDTx1000, nonAnnualDTx1000, fireBufferedListDT, formulaToFit)
  if (is.null(fires) || data.table::uniqueN(fires$year) < 4L) return(NULL)
  withr::local_seed(seed)
  fuels <- sub("^burn_", "", grep("^burn_", names(fires), value = TRUE))
  responses <- list(logSize = log(fires$obs), logSim = log(fires$simMedian),
                    resid = log(fires$obs / fires$simMedian))
  terms <- data.table::rbindlist(lapply(c("burn", "start"), function(at) {
    data.table::rbindlist(lapply(fuels, function(fuel) {
      fv <- fires[[paste0(at, "_", fuel)]]
      keep <- which(stats::complete.cases(fires$climYear, fv, fires$obs, fires$simMedian))
      ix <- split(seq_along(keep), fires$year[keep])
      if (length(ix) < 4L || stats::sd(fv[keep]) == 0 || stats::sd(fires$climYear[keep]) == 0) return(NULL)
      draws <- replicate(nBoot, sample.int(length(ix), replace = TRUE), simplify = FALSE)
      cz <- standardise(fires$climYear[keep]); fz <- standardise(fv[keep])
      sp <- spearmanByClimate(responses$logSize[keep], fv[keep], fires$climYear[keep], fires$year[keep], nBoot)
      data.table::rbindlist(lapply(names(responses), function(resp) {
        cbind(data.table::data.table(fuel = fuel, at = at, response = resp, nFires = length(keep),
                                     nYears = length(ix)),
              interactionTerm(responses[[resp]][keep], cz, fz, ix, draws),
              if (resp == "logSize") sp)
      }), fill = TRUE)
    }))
  }))
  list(terms = terms, climate = attr(fires, "climate"), nFires = nrow(fires),
       nYears = data.table::uniqueN(fires$year), nBoot = as.integer(nBoot), seed = seed)
}

#' @export
#' @rdname fuelClimateInteraction
fuelClimateFires <- function(sims, annualDTx1000, nonAnnualDTx1000, fireBufferedListDT, formulaToFit) {
  covCols <- spreadCovCols(spreadDesignCols(formulaToFit))
  clim <- intersect(covCols, names(annualDTx1000[[1]]))[1]
  fuels <- intersect(covCols, names(nonAnnualDTx1000[[1]]))
  if (is.na(clim) || !length(fuels)) return(NULL)
  nf <- grep("^nfLCC_", fuels, value = TRUE)
  tableYears <- as.integer(names(nonAnnualDTx1000))
  d <- data.table::as.data.table(sims)[!is.na(sim)]
  fire <- d[, list(obs = size[1], simMedian = stats::median(sim), initialLocus = initialLocus[1]),
            by = c("yr", "ids")]
  ## the fuels of some pixels: each covariate, and the non-forest share
  fuelsAt <- function(tab, pix) {
    i <- match(pix, tab$pixelID)
    out <- data.table::as.data.table(lapply(stats::setNames(fuels, fuels), function(f) tab[[f]][i] / 1000))
    if (length(nf) > 1L) data.table::set(out, NULL, "nonForest", rowSums(out[, nf, with = FALSE]))
    out
  }
  out <- data.table::rbindlist(lapply(unique(fire$yr), function(y) {
    f <- fire[yr == y]
    year <- as.integer(gsub("[[:alpha:]]+", "", y))
    ann <- annualDTx1000[[y]]
    tab <- nonAnnualDTx1000[[match(dataYearBefore(tableYears, year), tableYears)]]
    fb <- fireBufferedListDT[[y]]
    k <- which(fb$buffer == 1L & fb$ids %in% f$ids)
    b <- data.table::data.table(ids = fb$ids[k], pixelID = fb$pixelID[k])
    b <- cbind(b, fuelsAt(tab, b$pixelID), climBurn = ann[[clim]][match(b$pixelID, ann$pixelID)] / 1000)
    fuelNames <- setdiff(names(b), c("ids", "pixelID", "climBurn"))
    burn <- b[, c(list(climBurn = stats::median(climBurn, na.rm = TRUE)),
                  lapply(.SD, mean, na.rm = TRUE)), by = "ids", .SDcols = fuelNames]
    data.table::setnames(burn, fuelNames, paste0("burn_", fuelNames))
    start <- fuelsAt(tab, f$initialLocus)
    data.table::setnames(start, paste0("start_", names(start)))
    cbind(f[, list(yr, ids, year = year, obs, simMedian, climYear = stats::median(ann[[clim]], na.rm = TRUE) / 1000)],
          start)[burn, on = "ids"]
  }), fill = TRUE)
  data.table::setattr(out, "climate", clim)
  out
}

standardise <- function(v) (v - mean(v)) / stats::sd(v)

## least squares coefficients; NA where the design is rank deficient (a bootstrap draw of too few years)
olsCoef <- function(X, y) {
  f <- .lm.fit(X, y)
  if (f$rank < ncol(X)) rep(NA_real_, ncol(X)) else f$coefficients
}

## One fuel, one response. `ix`: the rows of each year; `draws`: bootstrap draws, each a resample of the years
interactionTerm <- function(y, cz, fz, ix, draws) {
  X <- cbind(1, cz, fz, cz * fz)
  boot <- vapply(draws, function(d) {
    k <- unlist(ix[d], use.names = FALSE)
    olsCoef(X[k, , drop = FALSE], y[k])
  }, numeric(4))
  q <- stats::quantile(boot[4, ], c(0.025, 0.975), na.rm = TRUE, names = FALSE)
  ## leave one year out: the squared error of the additive and the interaction model on the year left out
  cv <- vapply(seq_along(ix), function(j) {
    te <- ix[[j]]; tr <- unlist(ix[-j], use.names = FALSE)
    sq <- function(cols) sum((y[te] - X[te, cols, drop = FALSE] %*% olsCoef(X[tr, cols, drop = FALSE], y[tr]))^2)
    c(n = length(te), additive = sq(1:3), interaction = sq(1:4))
  }, numeric(3))
  ok <- stats::complete.cases(t(cv))
  data.table::data.table(
    bClimMedian = stats::median(boot[2, ], na.rm = TRUE), bFuelMedian = stats::median(boot[3, ], na.rm = TRUE),
    bIntMedian = stats::median(boot[4, ], na.rm = TRUE), bIntLo = q[1], bIntHi = q[2],
    pIntPositive = mean(boot[4, ] > 0, na.rm = TRUE),
    mseAdditive = sum(cv[2, ok]) / sum(cv[1, ok]), mseInteraction = sum(cv[3, ok]) / sum(cv[1, ok]),
    yearsImproved = sum(cv[3, ok] < cv[2, ok]), yearsTested = sum(ok))
}

## Spearman correlation of size with the fuel in the hottest third of the years and in the rest: the median of
## `nBoot` draws that resample the years within each, and their difference
spearmanByClimate <- function(size, fuel, clim, year, nBoot) {
  yearClim <- tapply(clim, year, `[`, 1L)
  hot <- clim >= stats::quantile(yearClim, 2 / 3)
  rho <- function(h) {
    ix <- split(which(hot == h), year[hot == h])
    replicate(nBoot, {
      k <- unlist(ix[sample.int(length(ix), replace = TRUE)], use.names = FALSE)
      suppressWarnings(stats::cor(size[k], fuel[k], method = "spearman"))
    })
  }
  bh <- rho(TRUE); br <- rho(FALSE); dif <- bh - br
  q <- stats::quantile(dif, c(0.025, 0.975), na.rm = TRUE, names = FALSE)
  data.table::data.table(
    rhoHotMedian = stats::median(bh, na.rm = TRUE), rhoRestMedian = stats::median(br, na.rm = TRUE),
    nHot = sum(hot), nRest = sum(!hot), rhoDiffMedian = stats::median(dif, na.rm = TRUE), rhoDiffLo = q[1],
    rhoDiffHi = q[2], pHotAboveRest = mean(dif > 0, na.rm = TRUE))
}
