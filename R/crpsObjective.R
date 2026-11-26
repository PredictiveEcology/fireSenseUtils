#' Weights of the CRPS spread objective
#'
#' The multipliers of the terms of `tests = "crps"` in [.objfunSpreadFit()]. They were set by the
#' 2026-10-28 profile of the objective on ELFs 6.6.2, 14.3 and 11.3 (the minimum of the combined score
#' against the simulated/observed area ratio, "W3" scheme: each weight is one over the median, over the
#' three ELFs, of its term's range for area ratios 0.5-2, with the total-area and top-10 terms doubled),
#' relative to the per-fire term = 1. Three ELFs and one seed: tune them with the `crpsWeights` argument.
#'
#' @return A named list: `fire` (per-fire CRPS on log size, [crpsFireSizes()]), `adTest` ([adStatistic()]),
#'   `areaDist` ([areaWeightedCvM()]), `totalArea` ([crpsTotalArea()]), `top10` ([crpsTopFires()]) and
#'   `annualArea` ([crpsAnnualArea()]), `topShare` ([crpsTopShare()]) and `climateSlope` ([crpsClimateSlope()]);
#'   the last three are 0, so their terms are left out unless a weight is given.
#' @export
crpsWeightsDefault <- function() {
  list(fire = 1, adTest = 0.34, areaDist = 2.6, totalArea = 164, top10 = 16.5, annualArea = 0, topShare = 0, climateSlope = 0)
}

#' Ensemble CRPS of a sample against one observation
#'
#' The continuous ranked probability score of the empirical distribution of `x` at `y`,
#' `E|X - y| - 0.5 E|X - X'|`: a proper score, 0 when `x` is the point mass at `y`. The same value as
#' `scoringRules::crps_sample(y, x)`; `E|X - X'|` is found from the sorted sample, not all pairs.
#'
#' @param x Numeric vector: the ensemble.
#' @param y A single number: the observation.
#' @return A single non-negative number.
#' @keywords internal
crpsEnsemble <- function(x, y) {
  n <- length(x)
  xs <- sort(x)
  mean(abs(x - y)) - sum((2 * seq_len(n) - n - 1) * xs) / n^2
}

#' CRPS of one fire's replicates on log size
#'
#' @param size Observed size of the fire (pixels).
#' @param sim Simulated sizes of its replicates (pixels).
#' @param cut Logical, as `sim`: was the replicate cut at an open buffer edge. A cut replicate would have
#'   kept burning, so its size is "at least this big": it is scored as `max(sim, size)`, never below the
#'   observed size. (The per-fire likelihood gives a cut replicate no density at the observed size; this
#'   is the same, it never counts as a hit, and one that passed the observed size counts as larger.)
#' @return A single non-negative number; a fire with no replicate is scored as 1 pixel.
#' @keywords internal
crpsFireSizes <- function(size, sim, cut = rep(FALSE, length(sim))) {
  if (!length(sim)) sim <- 1
  sim[cut] <- pmax(sim[cut], size)
  crpsEnsemble(log(sim), log(size))
}

#' CRPS of the total area burned
#'
#' @param repSizes A `data.table` with columns `rep` and `N`: every simulated fire of every fit year
#'   (pixels), one row each.
#' @param obsSizes Observed sizes of the fitted fires, all fit years (pixels).
#' @param Nreps Number of replicates; a replicate with no fire burned nothing.
#' @return A single non-negative number: [crpsEnsemble()] of log total area (at least 1 pixel).
#' @keywords internal
crpsTotalArea <- function(repSizes, obsSizes, Nreps) {
  tot <- numeric(Nreps)
  a <- repSizes[, list(A = sum(N)), by = "rep"]
  tot[a$rep] <- a$A
  crpsEnsemble(log(pmax(tot, 1)), log(sum(obsSizes)))
}

#' Each replicate's area burned in each fit year
#'
#' @param repSizes A `data.table` with columns `year`, `rep` and `N`: every simulated fire of every fit
#'   year (pixels), one row each.
#' @param years The fit years (character).
#' @param Nreps Number of replicates; a replicate with no fire in a year burned nothing.
#' @return A matrix, `Nreps` rows (replicates) by one column per year (named by year).
#' @keywords internal
annualTotals <- function(repSizes, years, Nreps) {
  out <- vapply(years, function(y) {
    tot <- numeric(Nreps)
    a <- repSizes[repSizes$year == y, list(A = sum(N)), by = "rep"]
    tot[a$rep] <- a$A
    tot
  }, numeric(Nreps))
  matrix(out, nrow = Nreps, dimnames = list(NULL, years)) # one replicate or one year: vapply gives a vector
}

#' Year weights of the annual-area term from the year's climate
#'
#' `"equal"`: 1 for every year. `"climate"`: `max(abs(z), floor)` for each year's standardised climate `z`, rescaled
#' to average 1, so years with extreme climate (either way) count more, an average year still counts (`floor`),
#' and the term's scale is unchanged on average. Climate, not area burned, because a forecast knows the climate
#' of a year and not its area.
#'
#' @param z Named numeric: the standardised climate of each fit year ([spreadYearClimate()]); can be `NULL` for `"equal"`.
#' @param type `"equal"` or `"climate"`.
#' @param floor Smallest `abs(z)` before rescaling, default 0.25.
#' @return A numeric vector, one weight per year (names of `z` kept), averaging 1. The one rule for the spread
#'   objective's annual-area term and the ignition and escape year terms.
#' @export
annualYearWeights <- function(z, type = c("equal", "climate"), floor = 0.25) {
  type <- match.arg(type)
  if (identical(type, "equal")) return(rep(1, length(z)))
  w <- pmax(abs(z), floor)
  w / mean(w)
}

#' CRPS of the area burned in each year
#'
#' For each fit year, the CRPS of the log of that year's total area (one per replicate; a replicate with no
#' fire that year burned nothing, scored as 1 pixel) against the log of the year's observed area. Unlike
#' [crpsTotalArea()], which pools the years, a fit that burns the mean every year does not score well.
#'
#' @param obsYear Named numeric: observed area burned (sum of the fitted fire sizes, pixels) per fit year,
#'   named by year.
#' @inheritParams annualTotals
#' @param yearWeights `"equal"` (default) or `"climate"` ([annualYearWeights()]); `"climate"` needs `z`.
#' @param z Named numeric: the standardised climate of each fit year, as `obsYear`; for `"climate"`.
#' @return A single non-negative number: the mean over years (not the sum, so its weight is on the scale of
#'   `totalArea`) of [crpsEnsemble()] of the log annual area (at least 1 pixel), times the year's weight.
#' @keywords internal
crpsAnnualArea <- function(repSizes, obsYear, Nreps, yearWeights = c("equal", "climate"), z = NULL) {
  yearWeights <- match.arg(yearWeights)
  yrs <- names(obsYear)
  tot <- log(pmax(annualTotals(repSizes, yrs, Nreps), 1))
  score <- vapply(seq_along(yrs), function(i) crpsEnsemble(tot[, i], log(max(obsYear[[i]], 1))), numeric(1))
  if (identical(yearWeights, "climate")) {
    if (is.null(z)) stop("yearWeights = \"climate\" needs the years' climate `z`: the year climate covariate is not in the spread data and no `yearClimateZ` was given", call. = FALSE)
    score <- score * annualYearWeights(z[yrs], yearWeights)
  }
  mean(score)
}

#' CRPS of the slope of log annual area on the year's climate
#'
#' For each replicate, the least-squares slope of the log of its annual area (at least 1 pixel) on the fit
#' years' standardised climate `z`; scored with [crpsEnsemble()] against the same slope of the observed
#' annual areas. A simulated response that is flat in climate, when the observed one rises with it, is
#' penalised; [crpsAnnualArea()] scores each year alone and does not see this.
#'
#' @inheritParams crpsAnnualArea
#' @param z Named numeric: the standardised climate of each fit year, as `obsYear`.
#' @return A single non-negative number; 0 with fewer than 2 years or no spread in `z`. An error without `z`.
#' @keywords internal
crpsClimateSlope <- function(repSizes, obsYear, Nreps, z) {
  if (is.null(z)) stop("the climate slope needs the years' climate `z`: the year climate covariate is not in the spread data and no `yearClimateZ` was given", call. = FALSE)
  yrs <- names(obsYear)
  z <- z[yrs]
  zc <- z - mean(z)
  if (length(yrs) < 2L || sum(zc^2) == 0) return(0)
  slope <- function(y) as.numeric(crossprod(zc, y - mean(y))) / sum(zc^2)
  sim <- apply(log(pmax(annualTotals(repSizes, yrs, Nreps), 1)), 1, slope)
  crpsEnsemble(sim, slope(log(pmax(as.numeric(obsYear), 1))))
}

#' CRPS of the largest fires
#'
#' The observed `k` largest fitted fires (all fit years pooled), sorted, against each replicate's `k`
#' largest, sorted: one [crpsEnsemble()] per rank on log size across replicates, summed over ranks. With
#' fewer than `k` observed fires, all of them are used; a replicate with fewer fires than ranks has 1 pixel
#' there.
#'
#' @inheritParams crpsTotalArea
#' @param k Number of fires, default 10.
#' @return A single non-negative number.
#' @keywords internal
crpsTopFires <- function(repSizes, obsSizes, Nreps, k = 10L) {
  obsTop <- utils::head(sort(obsSizes, decreasing = TRUE), k)
  k <- length(obsTop)
  top <- matrix(1, Nreps, k)
  for (r in unique(repSizes$rep)) {
    s <- utils::head(sort(repSizes$N[repSizes$rep == r], decreasing = TRUE), k)
    top[r, seq_along(s)] <- s
  }
  sum(vapply(seq_len(k), function(i) crpsEnsemble(log(top[, i]), log(obsTop[i])), numeric(1)))
}

#' CRPS of the share of the area burned by the largest fires
#'
#' For each replicate, the share of its total simulated area (all fit years) burned by its `k` largest
#' fires (pooled over the fit years), against the same share of the observed fitted fires. The CRPS is on
#' the logit of the share, so it is scale-free: multiplying every fire size by a constant leaves it
#' unchanged. It scores the weight of the tail, separately from how much burns ([crpsTotalArea()]).
#'
#' @inheritParams crpsTotalArea
#' @param k Number of fires, default 10.
#' @return A single non-negative number: [crpsEnsemble()] of the logit of the shares, each clamped to
#'   `[1e-6, 1 - 1e-6]` so the logit is finite. A replicate with no fire has share 1 (all of its, zero,
#'   area is in the largest fires), and so does an observation with fewer than `k` fires.
#' @keywords internal
crpsTopShare <- function(repSizes, obsSizes, Nreps, k = 10L) {
  share <- function(x) if (sum(x) > 0) sum(utils::head(sort(x, decreasing = TRUE), k)) / sum(x) else 1
  s <- rep(1, Nreps)
  for (r in unique(repSizes$rep)) s[r] <- share(repSizes$N[repSizes$rep == r])
  lg <- function(p) stats::qlogis(pmin(pmax(p, 1e-6), 1 - 1e-6))
  crpsEnsemble(lg(s), lg(share(obsSizes)))
}
