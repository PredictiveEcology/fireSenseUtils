#' Weights of the CRPS spread objective
#'
#' The multipliers of the terms of `tests = "crps"` in [.objfunSpreadFit()]. They were set by the
#' 2026-10-28 profile of the objective on ELFs 6.6.2, 14.3 and 11.3 (the minimum of the combined score
#' against the simulated/observed area ratio, "W3" scheme: each weight is one over the median, over the
#' three ELFs, of its term's range for area ratios 0.5-2, with the total-area and top-10 terms doubled),
#' relative to the per-fire term = 1. Three ELFs and one seed: tune them with the `crpsWeights` argument.
#'
#' @return A named list: `fire` (per-fire CRPS on log size, [crpsFireSizes()]), `adTest` ([adStatistic()]),
#'   `areaDist` ([areaWeightedCvM()]), `totalArea` ([crpsTotalArea()]) and `top10` ([crpsTopFires()]).
#' @export
crpsWeightsDefault <- function() {
  list(fire = 1, adTest = 0.34, areaDist = 2.6, totalArea = 164, top10 = 16.5)
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
