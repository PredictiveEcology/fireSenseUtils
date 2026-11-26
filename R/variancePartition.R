## How much of the fire signal each part of the spread model carries (climate, fuel, young stands, non-forest, ...), for two
## signals: between years (which years burn more) and within years (where, or how large, it burns). Two ways:
## * spreadVariancePartition(): from one fitted ledger row and the fit's covariates, no new fit. The linear predictor is split
##   into each group's contribution and its variance between and within years is shared among the groups.
## * commonalityPartition(): from fits of one ELF with no covariate, climate only, fuel only and both; what each fit explains is
##   split into unique climate, unique fuel, shared and unexplained (commonalityScores() reads it from the ledger rows).

utils::globalVariables(c("part", "share", "signal"))

## The groups of spreadTermGroups(), in the order they are reported
.spreadTermGroupOrder <- c("climate", "fuel", "youngAge", "nonForest", "other")

#' Fuel covariates on the scale of the fit
#'
#' A fit made on linear fuel biomass records it in its ranges: such a fuel's entry of `covMinMax_spread` is
#' [fuelLinearRange] ([isLinearFuelRange()]). Each such column of `dt` is put on the linear scale ([fuelLogToLinear()]); the
#' other columns, and a fuel of a fit made on the log scale, are left as they are.
#'
#' @param dt A `data.table` of covariates as `fireSenseCovariatesCreate()` makes them (fuel as [logMinB()]). Changed by
#'   reference.
#' @param ranges The fit's ranges, `covMinMax_spread` ([ledgerCovMinMax()]): a `data.table` or list, one `c(min, max)` per
#'   covariate.
#' @return `dt`.
#' @export
linearFuelCovs <- function(dt, ranges) {
  linear <- names(ranges)[vapply(ranges, isLinearFuelRange, logical(1))]
  for (cn in intersect(linear, names(dt))) data.table::set(dt, NULL, cn, fuelLogToLinear(dt[[cn]]))
  dt
}

#' The groups of the spread model's terms
#'
#' Each design column of a spread fit ([spreadDesignCols()], or the coefficient names of a ledger row) in one group:
#' `"climate"` (the climate covariates, and every interaction with one of them), `"fuel"` (fuel biomass, [spreadFuelCols()]:
#' `dom_agb_*`, `sec_agb_*`, `treedWetland_agb`, or the fuel classes of a species-mode fit), `"youngAge"`, `"nonForest"`
#' (`nfLCC_*`) and `"other"` (the `treedWetland` indicator, and an interaction of two groups that are not climate). An
#' interaction of two terms of one group is in that group. The intercept is a constant and is in no group.
#'
#' @param designCols Character; the design columns.
#' @param climateCols Character; the climate covariates among them (the annual covariates that are not indicators, as
#'   [spreadClimateCols()]; [spreadDesignByYear()] records them in its attribute `"climateCols"`).
#' @return A named character vector: one group per design column (the intercept left out), named by the column.
#' @export
spreadTermGroups <- function(designCols, climateCols) {
  groupOf <- function(cn) {
    if (cn %in% climateCols) return("climate")
    if (!length(spreadIndicatorCols(cn))) return("fuel")             # not an indicator: fuel biomass (spreadFuelCols())
    if (cn == youngAgeTxt) "youngAge" else if (cn == treedWetlandTxt) "other" else "nonForest"
  }
  cols <- spreadCovCols(designCols)
  stats::setNames(vapply(cols, function(cn) {
    parts <- strsplit(cn, ":", fixed = TRUE)[[1]]
    g <- vapply(parts, groupOf, character(1))
    if ("climate" %in% g) "climate" else if (length(unique(g)) == 1L) g[[1]] else "other"
  }, character(1)), cols)
}

#' The fit's design matrix of every pixel-year, from the covariate tables of `fireSense_dataPrepFit`
#'
#' The covariates of a fitted ledger row on the scale its coefficients multiply: fuel biomass on the linear scale where the
#' row's `covMinMax_spread` says so ([isLinearFuelRange()], [fuelLogToLinear()]), rescaled with the row's ranges
#' ([ledgerCovMinMax()]), `youngAge` made mutually exclusive with fuel and the indicators ([youngAgeExclusiveCols()]),
#' centred with the row's centres ([ledgerCovCentre()]), the interactions and the intercept added. Each year goes through the
#' function the objective uses ([spreadProbFromIntegerCovs()], via `eachSpreadYearCovs()`), with the vegetation table of the
#' latest data year before it.
#'
#' The rows are those of `annualCovariates`: the pixels in the fire buffers of each year, or every burnable pixel when
#' `fireSense_dataPrepFit` ran with `wholeLandscapeSpreadCovariates`. With buffers the pixels differ between years, so the
#' between-year signal of fuel is then partly which pixels the buffers hold.
#'
#' @param row A one-row spread ledger row, as [latestSpreadFits()] or `sim$studyAreaWithSpreadParams` gives it.
#' @param annualCovariates,nonAnnualCovariates `fireSense_dataPrepFit`'s `fireSense_annualSpreadFitCovariates` (named by
#'   year) and `fireSense_nonAnnualSpreadFitCovariates` (named by data year), on their own scale (fuel as [logMinB()]). Not
#'   changed.
#' @return A `data.table`: `year` (integer), `pixelID` and the row's design columns. Attribute `"climateCols"`: its climate
#'   covariates (annual, not indicators), for [spreadTermGroups()].
#' @export
spreadDesignByYear <- function(row, annualCovariates, nonAnnualCovariates) {
  design <- ledgerDesignCols(names(as.data.frame(row$params[[1L]])))
  sources <- spreadSourceCols(design)
  covMinMax <- ledgerCovMinMax(row, 1L)
  if (is.null(covMinMax)) stop("spreadDesignByYear(): the ledger row has no covMinMax_spread", call. = FALSE)
  covMinMax <- covMinMax[, intersect(names(covMinMax), sources), with = FALSE]
  ## a copy of each table, of its columns the design uses, fuel on the scale of the fit, x1000 as the objective has them
  prepare <- function(tables) toX1000(lapply(tables, function(dt)
    linearFuelCovs(data.table::as.data.table(dt)[, intersect(c("pixelID", sources), names(dt)), with = FALSE], covMinMax)))
  annualX <- prepare(annualCovariates)
  nonAnnualX <- prepare(nonAnnualCovariates)
  climateCols <- annualCovariateCols(sources, names(annualX[[1L]]))
  mutuallyExclusive <- if (youngAgeTxt %in% sources)
    youngAgeExclusiveCols(sources, fuelCols = spreadFuelCols(setdiff(sources, climateCols)))
  byYear <- eachSpreadYearCovs(annualX, nonAnnualX, design, covMinMax, mutuallyExclusive, ledgerCovCentre(row, 1L),
                               function(covs) covs[, c("pixelID", design), with = FALSE])
  out <- data.table::rbindlist(byYear, idcol = yearTxt)
  data.table::set(out, NULL, yearTxt, yearOfName(out[[yearTxt]]))
  data.table::setattr(out, "climateCols", climateCols)
  out
}

#' How much of the between-year and the within-year variation of the spread model each group of terms carries
#'
#' For each parameter set of a fitted ledger row, the linear predictor of every pixel-year, `X = sum(coef_t * x_t)` over the
#' design columns of `design`, is split into the contribution of each group of terms ([spreadTermGroups()]), `C_g`. The
#' variance of `X` over the pixel-years is the sum of a between-year part (the variance of the year means, each year weighted
#' by its pixels) and a within-year part (around each year's mean): signal `"between"` (which years burn more) and
#' `"within"` (where it burns in a year). For each signal and group:
#'
#' * `own`: the variance of `C_g` alone, as a share of the signal's variance;
#' * `covariance`: the group's part of the covariances, `sum(cov(C_g, C_h))` over the other groups `h`, as a share (it can be
#'   negative). `own + covariance` over the groups adds to 1;
#' * `share`: the group's Shapley (LMG) share of the signal: the R-squared of `X` on the groups' contributions, averaged
#'   over the orders in which the groups can be added ([lmgShares()]). Never negative; the groups' shares add to 1. A
#'   group whose contribution does not vary has 0.
#'
#' With climate outside or in the ceiling, or the activation (`spreadClimate` other than `"predictor"`; column
#' `spreadClimate`), climate and the other terms act on different parts of the link, and `X` here is the sum of both linear
#' predictors, not one quantity the link applies to. The partition is of the linear predictor, not of the spread
#' probability or the area burned: the link saturates.
#'
#' @param row A one-row spread ledger row ([latestSpreadFits()]).
#' @param design The design of each pixel-year, [spreadDesignByYear()]: `year` and the row's design columns.
#' @param groups The group of each design column, [spreadTermGroups()].
#' @param byParamSet Logical. `FALSE` (default): the median over the row's parameter sets. `TRUE`: one set of rows per
#'   parameter set, with column `paramSet`.
#' @return A `data.frame`, one row per signal and group: `signal`, `part` (the group), `share`, `own`, `covariance`,
#'   `variance` (the signal's variance of `X`), `signalShare` (the signal's share of the variance of `X`) and
#'   `spreadClimate`.
#' @export
spreadVariancePartition <- function(row, design, groups = spreadTermGroups(names(design), attr(design, "climateCols")),
                                    byParamSet = FALSE) {
  params <- as.data.frame(row$params[[1L]])
  cols <- spreadCovCols(ledgerDesignCols(names(params)))   # the intercept is a constant, with no variance to share
  noCol <- setdiff(cols, names(design))
  if (length(noCol))
    stop("spreadVariancePartition(): `design` has no column for the coefficient(s) ", paste(noCol, collapse = ", "),
         call. = FALSE)
  noGroup <- setdiff(cols, names(groups))
  if (length(noGroup))
    stop("spreadVariancePartition(): `groups` has no group for ", paste(noGroup, collapse = ", "), call. = FALSE)
  parts <- intersect(c(.spreadTermGroupOrder, unique(groups[cols])), groups[cols])
  mat <- as.matrix(data.table::as.data.table(design)[, cols, with = FALSE])
  year <- design[[yearTxt]]
  ## one column per group: 1 where the term is in it
  member <- outer(groups[cols], parts, `==`) * 1
  each <- lapply(seq_len(NROW(params)), function(k) {
    contrib <- mat %*% (unlist(params[k, cols]) * member)          # pixel-years x groups: C_g
    out <- signalCovariances(contrib, year)
    total <- sum(vapply(out, sum, numeric(1)))
    do.call(rbind, lapply(names(out), function(signal) {
      S <- out[[signal]]
      v <- sum(S)
      ok <- v > 0
      data.frame(paramSet = k, signal = signal, part = parts, share = lmgShares(S),
                 own = if (ok) diag(S) / v else NA_real_, covariance = if (ok) (rowSums(S) - diag(S)) / v else NA_real_,
                 variance = v, signalShare = v / total)
    }))
  })
  out <- do.call(rbind, each)
  if (!isTRUE(byParamSet)) {
    num <- c("share", "own", "covariance", "variance", "signalShare")
    out <- stats::aggregate(out[num], out[c("signal", "part")], stats::median)
    out <- out[order(out$signal, match(out$part, parts)), ]
  }
  out$spreadClimate <- ledgerSpreadClimate(row, 1L)
  rownames(out) <- NULL
  out
}

## The between-year and the pooled within-year covariance matrices of the columns of `contrib` (pixel-years x parts), each
## over all pixel-years (divided by their number), so the two add to the covariance of `contrib`.
signalCovariances <- function(contrib, year) {
  n <- nrow(contrib)
  nYear <- as.vector(rowsum(rep(1, n), year))
  yearMean <- rowsum(contrib, year) / nYear
  dev <- sweep(yearMean, 2L, colMeans(contrib)) * sqrt(nYear)
  within <- contrib - yearMean[match(year, as.numeric(rownames(yearMean))), , drop = FALSE]
  list(between = crossprod(dev) / n, within = crossprod(within) / n)
}

#' Shapley (LMG) shares of a sum among its parts
#'
#' `y` is the sum of parts with covariance matrix `Sigma`. The R-squared of the regression of `y` on a subset `S` of the
#' parts is `c_S' Sigma_SS^+ c_S / var(y)`, with `c = rowSums(Sigma)` (the covariances of each part with `y`) and `^+` the
#' pseudo-inverse. A part's share is its gain in R-squared averaged over the orders in which the parts can be added (the
#' Shapley value; Lindeman, Merenda and Gold, 1980). The shares add to 1 and none is negative.
#'
#' @param Sigma A square covariance matrix of the parts, with `dimnames` naming them (optional).
#' @return A numeric vector, one share per part; `NA` for all when `y` does not vary.
#' @export
lmgShares <- function(Sigma) {
  p <- nrow(Sigma)
  vy <- sum(Sigma)
  if (!(vy > 0)) return(stats::setNames(rep(NA_real_, p), rownames(Sigma)))
  cy <- rowSums(Sigma)
  pinv <- function(A) {
    s <- svd(A)
    keep <- s$d > max(s$d) * sqrt(.Machine$double.eps)
    s$v[, keep, drop = FALSE] %*% (t(s$u[, keep, drop = FALSE]) / s$d[keep])
  }
  masks <- 0:(2^p - 1)
  inSet <- function(m) bitwAnd(m, 2^(seq_len(p) - 1)) > 0
  r2 <- vapply(masks, function(m) {
    S <- inSet(m)
    if (!any(S)) return(0)
    drop(crossprod(cy[S], pinv(Sigma[S, S, drop = FALSE]) %*% cy[S])) / vy
  }, numeric(1))
  size <- vapply(masks, function(m) sum(inSet(m)), numeric(1))
  weight <- factorial(size) * factorial(pmax(p - size - 1, 0)) / factorial(p)   # used for sets without the part
  stats::setNames(vapply(seq_len(p), function(g) {
    bit <- 2^(g - 1)
    without <- masks[bitwAnd(masks, bit) == 0]
    sum(weight[without + 1] * (r2[without + bit + 1] - r2[without + 1]))
  }, numeric(1)), rownames(Sigma))
}

#' Unique climate, unique fuel, shared and unexplained, from fits with and without each
#'
#' Commonality analysis of two groups of terms. With `R0`, `Rc`, `Rf` and `Rcf` the share of a signal explained by the fit
#' with no covariate, climate only, fuel only and both, the signal is split into `baseline` (`R0`, what the fit with no
#' covariate explains: through the ignitions, for example), unique `climate` (`Rcf - Rf`), unique `fuel` (`Rcf - Rc`),
#' `shared` (`Rc + Rf - Rcf - R0`, what either explains; negative when each makes the other more useful) and
#' `unexplained` (`1 - Rcf`). The five add to 1.
#'
#' @param scores A `data.frame`, one row per fit and signal: `climate` and `fuel` (logical: does the fit have them),
#'   `explained` (the share of the signal the fit explains) and, optionally, `signal`. [commonalityScores()] makes it.
#' @return A `data.frame`: `signal`, `part` (`"baseline"`, `"climate"`, `"fuel"`, `"shared"`, `"unexplained"`) and `share`;
#'   and `annualSource` when `scores` has it.
#' @export
commonalityPartition <- function(scores) {
  scores <- as.data.frame(scores)
  if (is.null(scores$signal)) scores$signal <- "signal"
  do.call(rbind, lapply(split(scores, factor(scores$signal, unique(scores$signal))), function(s) {
    R <- function(cl, fu) {
      v <- s$explained[s$climate == cl & s$fuel == fu]
      if (length(v) != 1L)
        stop("commonalityPartition(): signal '", s$signal[1], "' needs one fit for each of: no covariate, climate only, ",
             "fuel only, climate and fuel", call. = FALSE)
      v
    }
    R0 <- R(FALSE, FALSE); Rc <- R(TRUE, FALSE); Rf <- R(FALSE, TRUE); Rcf <- R(TRUE, TRUE)
    out <- data.frame(signal = s$signal[1], part = c("baseline", "climate", "fuel", "shared", "unexplained"),
                      share = c(R0, Rcf - Rf, Rcf - Rc, Rc + Rf - Rcf - R0, 1 - Rcf))
    if (!is.null(s$annualSource)) out$annualSource <- s$annualSource[1]
    out
  }))
}

#' What each of four fits explains of the between-year and the fire-size signal, from their ledger rows
#'
#' * signal `"between"`: [annualAreaExplained()] of the observed annual areas by the fit's simulated median ones. From the
#'   hindcasts (`annualAreas`, [hindcastAnnualAreas()]) when all four fits have one, else from each row's in-sample fire-size
#'   validation (`spreadFitInSampleSizesTxt`, `byYear`, written by `fireSense_spreadFit`) for all four, so the four are
#'   always compared on the same years; column `annualSource` says which (`"hindcast"`, `"inSample"`);
#' * signal `"within"`: the reduction of the in-sample per-fire squared error of log10 size (`sizes$fireRMSE^2`) from the fit
#'   with no covariate, `1 - fireRMSE^2 / fireRMSE_none^2`. It is 0 for that fit by definition; a fit worse than it is
#'   negative. `annualSource` is `NA`.
#'
#' @param rows A named list of four one-row ledger rows, `none`, `climate`, `fuel` and `climateFuel`: one ELF's fits with no
#'   covariate, climate only, fuel only, and both.
#' @param annualAreas `NULL`, or a list named as `rows` (a fit may be missing): each fit's annual areas from its hindcasts,
#'   [hindcastAnnualAreas()] (`obs`, `simMedian`).
#' @return A `data.frame` for [commonalityPartition()]: `signal`, `fit`, `climate`, `fuel`, `explained` and `annualSource`.
#' @export
commonalityScores <- function(rows, annualAreas = NULL) {
  fits <- c("none", "climate", "fuel", "climateFuel")
  if (!all(fits %in% names(rows)))
    stop("commonalityScores(): `rows` must be named ", paste(fits, collapse = ", "), call. = FALSE)
  inSample <- lapply(rows[fits], function(r) {
    v <- if (spreadFitInSampleSizesTxt %in% names(r)) r[[spreadFitInSampleSizesTxt]][[1L]]
    if (is.null(v$sizes) || is.null(v$byYear))
      stop("commonalityScores(): a ledger row has no in-sample sizes (", spreadFitInSampleSizesTxt, ")", call. = FALSE)
    v
  })
  fireMSE <- vapply(inSample, function(v) v$sizes$fireRMSE[[1L]]^2, numeric(1))
  hindcast <- all(fits %in% names(Filter(NROW, annualAreas)))
  years <- if (hindcast) annualAreas[fits] else lapply(inSample, `[[`, "byYear")
  between <- vapply(years, function(y) annualAreaExplained(y$obs, y$simMedian), numeric(1))
  data.frame(signal = rep(c("between", "within"), each = 4L), fit = fits, climate = fits %in% c("climate", "climateFuel"),
             fuel = fits %in% c("fuel", "climateFuel"), explained = c(between, 1 - fireMSE / fireMSE[["none"]]),
             annualSource = rep(c(if (hindcast) "hindcast" else "inSample", NA_character_), each = 4L))
}

#' Observed and simulated area burned of each hindcast year
#'
#' From the tables the hindcast writes for each year (`fireSense_summary`: `hindcast_reps.csv`, `hindcast_observed.csv`),
#' row-bound over years and parameter sets.
#'
#' @param reps The simulated runs: `year` and `areaBurnedHa`, one row per run (rep of a parameter set) and year.
#' @param observed The observed year: `year` and `areaBurnedHa` (a year read more than once counts once).
#' @return A `data.frame`, one row per year with both: `yr`, `obs` and `simMedian` (the median over the runs), as the in-sample
#'   `byYear` of [commonalityScores()].
#' @export
hindcastAnnualAreas <- function(reps, observed) {
  observed <- as.data.frame(observed)[!duplicated(observed$year), ]
  sim <- stats::aggregate(list(simMedian = reps$areaBurnedHa), list(yr = reps$year), stats::median)
  out <- merge(data.frame(yr = observed$year, obs = observed$areaBurnedHa), sim, by = "yr")
  out[order(out$yr), ]
}

#' The share of the variation of observed annual area burned that simulated annual area explains
#'
#' The squared correlation, over years, of `log(1 + observed)` and `log(1 + simulated)` (areas in any one unit; the `1 +`
#' keeps years without fire), and 0 when the correlation is not positive or either does not vary.
#'
#' @param obs,sim Numeric; the observed and the simulated (e.g. median over runs) area burned of each year.
#' @return A number from 0 to 1.
#' @export
annualAreaExplained <- function(obs, sim) {
  lo <- log1p(obs); ls <- log1p(sim)
  if (length(lo) < 2L || stats::sd(lo) == 0 || stats::sd(ls) == 0) return(0)
  r <- stats::cor(lo, ls)
  if (r > 0) r^2 else 0
}

#' Bar chart of a signal partition
#'
#' @param x [spreadVariancePartition()] or [commonalityPartition()]: `signal`, `part`, `share`.
#' @param title Optional title.
#' @return A `ggplot`: one stacked bar per signal.
#' @export
plotSignalPartition <- function(x, title = NULL) {
  x <- as.data.frame(x)
  x$part <- factor(x$part, unique(x$part))
  ggplot2::ggplot(x, ggplot2::aes(x = signal, y = share, fill = part)) +
    ggplot2::geom_col() +
    ggplot2::labs(x = NULL, y = "share of the signal", fill = NULL, title = title)
}
