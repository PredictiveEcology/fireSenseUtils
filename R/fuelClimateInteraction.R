utils::globalVariables(c("climBurn", "simMedian"))

## fewer years than this have no bootstrap over years worth reading
interactionMinYears <- 4L

#' A screen of the observed fires for fuel x climate interactions
#'
#' A spread model with no interaction terms is additive in the logit of the spread probability: a fuel term has the
#' same effect in every year's climate. `fuelClimateInteractionScreen()` screens the fit's own observed fires (and the
#' fit's in-sample simulation of them) for a sign that it does not. For each fitted fire it takes the observed size,
#' the year's climate and the fuel the fire burned (and the fuel at its start pixel), and compares
#' `log(size) ~ climate + fuel` with `log(size) ~ climate + fuel + climate:fuel`, one model per fuel covariate. It uses
#' only what is in memory at the end of a fit and runs no simulation. It is a screen of fire sizes, not a spread
#' model: it neither fits nor tests a spread model with interaction terms (the `spreadInteractions` option of the
#' fit is a separate thing).
#'
#' * `fuelClimateFires()` is the per-fire table: observed size (`obs`), the median in-sample simulated size over the
#'   fit's members and runs (`simMedian`), the year's climate (`climYear`) and the fire's burned-pixel climate
#'   (`climBurn`), and one column per fuel, `burn_<fuel>` (mean over the burned pixels) and `start_<fuel>` (at the start
#'   pixel). The fuels are the non-annual covariates of the formula and, with two or more `nfLCC_*` groups,
#'   `nonForest` (their sum, the non-forest share). A fire year reads the vegetation table of the data year strictly
#'   before it ([dataYearBefore()]), as the fit does.
#' * `fuelClimateInteractionScreen()` fits the models to the observed size (`logSize`), the simulated size
#'   (`logSim`) and the residual `log(obs / simMedian)` (`resid`).
#'
#' The climate is the first annual covariate of the formula, by its base name: a `<C>_yr` / `<C>_dev` split of a
#' climate variable `C` (`fireSense_dataPrepFit`'s `spreadClimateSplit`) counts as `C`, its year's climate is the
#' `<C>_yr` value (the ELF median over the burnable pixels, the same for every pixel of the year) and `climBurn` adds
#' `<C>_dev`. Without a split it is the median of the climate column over the year's rows of `annualDTx1000`, which
#' are the pixels of that year's fire buffers (every flammable pixel only when the fit used
#' `wholeLandscapeSpreadCovariates`), not necessarily the ELF's burnable pixels. A name with a `:` (an interaction
#' term of the formula) is neither a climate nor a fuel; only columns present in the tables are used.
#'
#' The predictors are standardised once, over all fires, not within each leave-one-year-out fold. Least squares
#' predictions of the additive and the interaction model do not depend on the predictors' centring or scale (the
#' interaction model's span contains `1`, `climate`, `fuel` and their product however they are shifted and scaled), so
#' this changes only the reported coefficients (per SD of climate and fuel), not the errors compared.
#'
#' @param sims the in-sample simulations of the fit: `attr(DE, "fitSims")` of [runDEoptim()], the result of
#'   [simulateFireSizes()] (the columns `yr`, `ids`, `initialLocus`, `size` and `sim`).
#' @param annualDTx1000,nonAnnualDTx1000,fireBufferedListDT The covariate tables the fit used, as given to
#'   [.objfunSpreadFit()] (integers x 1000; standardising removes the scale).
#' @param formulaToFit The spread formula. Its first annual (climate) covariate is the climate; its non-annual
#'   covariates are the fuels (see Details).
#' @param nBoot Number of bootstrap draws over years. `0` skips the screen.
#' @param seed Seed of the draws (the session's random number stream is restored afterwards), so the result is
#'   the same for the same fit.
#'
#' @return `fuelClimateFires()`: a `data.table`, one row per fitted fire, or `NULL` when the formula has no climate
#'   covariate or no fuel covariate. `fuelClimateInteractionScreen()`: `NULL` in that case, with fewer than 4 years
#'   with fires, or with `nBoot = 0`; else a list: `terms`, one row per fuel, `at` (`"burn"` or `"start"`) and
#'   `response`, with `nFires` and `nYears`; `bClimMedian`, `bFuelMedian`, and the interaction coefficient
#'   (`bIntMedian`, `bIntLo`, `bIntHi`: the median and 2.5-97.5% range over `nBoot` bootstrap draws that resample
#'   years, and `pIntPositive`, their share above 0); `mseAdditive` and `mseInteraction`, the mean squared error
#'   of leaving each year out in turn, `yearsImproved` of `yearsTested` years where the interaction model is better;
#'   and, for `logSize` only, the Spearman correlation of size with the fuel in the hottest third of years and in the
#'   rest (`rhoHotMedian`, `rhoRestMedian`, `nHot`, `nRest`), their difference `rhoDiffMedian`, `rhoDiffLo`,
#'   `rhoDiffHi` and `pHotAboveRest`, bootstrapped over years alike. The list also has `climate` (the covariate),
#'   `nFires`, `nYears`, `nBoot` and `seed`. It is what the ledger row's [spreadFitInteractionScreenTxt] column holds.
#' @export
#' @rdname fuelClimateInteraction
fuelClimateInteractionScreen <- function(sims, annualDTx1000, nonAnnualDTx1000, fireBufferedListDT, formulaToFit,
                                       nBoot = 999L, seed = 1L) {
  if (!isTRUE(nBoot > 0) || data.table::uniqueN(sims$yr[!is.na(sims$sim)]) < interactionMinYears) return(NULL)
  fires <- fuelClimateFires(sims, annualDTx1000, nonAnnualDTx1000, fireBufferedListDT, formulaToFit)
  if (is.null(fires)) return(NULL)
  withr::local_seed(seed)
  fuels <- sub("^burn_", "", grep("^burn_", names(fires), value = TRUE))
  responses <- list(logSize = log(fires$obs), logSim = log(fires$simMedian),
                    resid = log(fires$obs / fires$simMedian))
  terms <- data.table::rbindlist(lapply(c("burn", "start"), function(at) {
    data.table::rbindlist(lapply(fuels, function(fuel) {
      fv <- fires[[paste0(at, "_", fuel)]]
      keep <- which(stats::complete.cases(fires$climYear, fv, fires$obs, fires$simMedian))
      ix <- split(seq_along(keep), fires$year[keep])
      if (length(ix) < interactionMinYears || stats::sd(fv[keep]) == 0 || stats::sd(fires$climYear[keep]) == 0) return(NULL)
      draws <- bootYears(ix, nBoot)
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
  covCols <- grep(":", spreadCovCols(spreadDesignCols(formulaToFit)), fixed = TRUE, value = TRUE, invert = TRUE)
  clim <- climateColumns(intersect(covCols, names(annualDTx1000[[1]])), names(annualDTx1000[[1]]))
  fuels <- intersect(covCols, names(nonAnnualDTx1000[[1]]))
  if (is.null(clim) || !length(fuels)) return(NULL)
  nf <- grep("^nfLCC_", fuels, value = TRUE)
  tableYears <- as.integer(names(nonAnnualDTx1000))
  d <- data.table::as.data.table(sims)[!is.na(sim)]
  fire <- perFireSizes(d, stats::median)
  data.table::setnames(fire, "sim", "simMedian")
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
    b <- cbind(b, fuelsAt(tab, b$pixelID), climBurn = climateOf(ann, clim)[match(b$pixelID, ann$pixelID)])
    fuelNames <- setdiff(names(b), c("ids", "pixelID", "climBurn"))
    burn <- b[, c(list(climBurn = stats::median(climBurn, na.rm = TRUE)),
                  lapply(.SD, mean, na.rm = TRUE)), by = "ids", .SDcols = fuelNames]
    data.table::setnames(burn, fuelNames, paste0("burn_", fuelNames))
    start <- fuelsAt(tab, f$initialLocus)
    data.table::setnames(start, paste0("start_", names(start)))
    cbind(f[, list(yr, ids, year = year, obs, simMedian, climYear = stats::median(ann[[clim$year]], na.rm = TRUE) / 1000)],
          start)[burn, on = "ids"]
  }), fill = TRUE)
  data.table::setattr(out, "climate", clim$name)
  out
}

## `nBoot` bootstrap draws over years: each the rows of years resampled with replacement. `ix`: the rows of each year.
bootYears <- function(ix, nBoot) {
  replicate(nBoot, unlist(ix[sample.int(length(ix), replace = TRUE)], use.names = FALSE), simplify = FALSE)
}

## The climate of a formula: the base name, the column of the year's climate, and the `_dev` column of a split
## (`NULL`: none). `annualCovs`: the formula's annual covariates, `tableCols`: the columns of an annual table.
climateColumns <- function(annualCovs, tableCols) {
  if (!length(annualCovs)) return(NULL)
  base <- sub("_(yr|dev)$", "", annualCovs[1])
  yr <- paste0(base, "_yr"); dev <- paste0(base, "_dev")
  if (all(c(yr, dev) %in% tableCols)) return(list(name = base, year = yr, dev = dev))
  list(name = base, year = if (base %in% tableCols) base else annualCovs[1], dev = NULL)
}

## the climate of each row of a year's annual table (x1000 removed): the year's value plus, for a split, the departure
climateOf <- function(ann, clim) {
  (ann[[clim$year]] + if (!is.null(clim$dev)) ann[[clim$dev]] else 0) / 1000
}

standardise <- function(v) (v - mean(v)) / stats::sd(v)

## least squares coefficients; NA where the design is rank deficient (a bootstrap draw of too few years)
olsCoef <- function(X, y) {
  f <- stats::.lm.fit(X, y)
  if (f$rank < ncol(X)) rep(NA_real_, ncol(X)) else f$coefficients
}

## One fuel, one response. `ix`: the rows of each year; `draws`: `bootYears()`
interactionTerm <- function(y, cz, fz, ix, draws) {
  X <- cbind(1, cz, fz, cz * fz)
  boot <- vapply(draws, function(d) {
    olsCoef(X[d, , drop = FALSE], y[d])
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
    vapply(bootYears(split(which(hot == h), year[hot == h]), nBoot), function(k)
      suppressWarnings(stats::cor(size[k], fuel[k], method = "spearman")), numeric(1))
  }
  bh <- rho(TRUE); br <- rho(FALSE); dif <- bh - br
  q <- stats::quantile(dif, c(0.025, 0.975), na.rm = TRUE, names = FALSE)
  data.table::data.table(
    rhoHotMedian = stats::median(bh, na.rm = TRUE), rhoRestMedian = stats::median(br, na.rm = TRUE),
    nHot = sum(hot), nRest = sum(!hot), rhoDiffMedian = stats::median(dif, na.rm = TRUE), rhoDiffLo = q[1],
    rhoDiffHi = q[2], pHotAboveRest = mean(dif > 0, na.rm = TRUE))
}
