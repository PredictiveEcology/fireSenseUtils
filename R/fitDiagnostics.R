#' Diagnostics of a fitted spread model
#'
#' The checks made by hand on the FireSense phase-2 fits (September 2026), as functions, so that
#' `fireSense_SpreadFit` can make them after every fit.
#'
#' * `simulateFireSizes()` simulates the observed fires from each of several parameter sets, with
#'   the objective's size cap and "too burny" gate lifted, and returns the simulated sizes.
#' * `scoreFireSizes()` compares those with the observed sizes: per fire, per year, and by
#'   quantile. It never uses the objective's value.
#' * `linkSaturation()` reports how many pixel-years sit at the spread-probability ceiling. When
#'   most do, fire size is decided by the ceiling alone ("too many medium, too few large fires").
#' * `coefIdentifiability()` measures, for each covariate coefficient, how tightly the final
#'   population pins it down.
#' * `profileCoefficients()` sets each coefficient in turn to 0 and to values across the
#'   population's range, other coefficients held at the best member, and re-scores.
#' * `identifiedInIsolation()` combines the two: a coefficient is identified in isolation when
#'   the population pins its sign (`zPop >= zMin`) and dropping it (setting it to 0) makes the
#'   objective clearly worse.
#' * `fitConvergence()` summarises the objective across the fit's generations.
#'
#' @name fitDiagnostics
NULL

#' @param p Numeric. Spread probabilities.
#' @param ceiling The link's upper asymptote, `maxAsymptote`.
#' @param tol A probability within `tol` of `ceiling` counts as at the ceiling.
#' @param binWidth Width of the bins in which `p` is counted, so summaries of different years can
#'   be added up and still give quantiles.
#'
#' @return `spreadProbSummary()`: a list with `n`, `atCeiling` (a count), `counts` (per bin) and
#'   `binWidth`.
#' @export
#' @rdname fitDiagnostics
spreadProbSummary <- function(p, ceiling, tol = 0.002, binWidth = 0.001) {
  nBins <- ceiling(1 / binWidth) + 1L
  list(n = length(p), atCeiling = sum(p >= ceiling - tol),
       counts = tabulate(pmin(floor(p / binWidth) + 1L, nBins), nbins = nBins),
       binWidth = binWidth)
}

#' @param x A list of `spreadProbSummary()` results with the same `binWidth`.
#' @return `combineSpreadProbSummaries()`: their sum, as one summary.
#' @export
#' @rdname fitDiagnostics
combineSpreadProbSummaries <- function(x) {
  x <- Filter(Negate(is.null), x)
  if (!length(x)) return(NULL)
  list(n = sum(vapply(x, `[[`, numeric(1), "n")),
       atCeiling = sum(vapply(x, `[[`, numeric(1), "atCeiling")),
       counts = Reduce(`+`, lapply(x, `[[`, "counts")),
       binWidth = x[[1]]$binWidth)
}

#' @param pop A matrix or `data.table` of parameter sets, one per row, with column names.
#' @param fn The objective; [.objfunSpreadFit()].
#' @param fnArgs Named list of further arguments to `fn`, as for [rescorePopulation()]. The fit data
#'   (`landscape`, `annualDTx1000`, ...) are taken from the workers' global environment if absent,
#'   so a held-out prediction passes its own years here.
#' @param cl A cluster, or `NULL` to run in this session.
#' @param seed Every parameter set is evaluated with the same seeds (`seed`, or `seed + i` for
#'   replicate `i` of the profile), so differences between them are not luck.
#'
#' @return `simulateFireSizes()`: a `data.table`, one row per parameter set x fire x replicate:
#'   `member` (row of `pop`), `yr`, `rep`, `initialLocus`, `ids`, simulated `sim` and observed
#'   `size`, in pixels. A year the objective will not simulate even with its gate lifted keeps its
#'   fires with `sim` `NA`. `attr(, "spreadProb")` is a list of `spreadProbSummary()`s, one per
#'   member.
#' @export
#' @rdname fitDiagnostics
simulateFireSizes <- function(pop, fn = .objfunSpreadFit, fnArgs = list(), cl = NULL, seed = 1L) {
  pop <- as.matrix(pop)
  fnArgs <- utils::modifyList(fnArgs, list(returnSims = TRUE, capSizes = FALSE, thresh = Inf,
                                           lanscape1stQuantileThresh = Inf))
  jobList <- lapply(seq_len(NROW(pop)), function(i) list(par = pop[i, ], seed = as.integer(seed)))
  sims <- if (is.null(cl)) {
    lapply(jobList, .rescoreJob, fn = fn, fnArgs = fnArgs)
  } else {
    parallel::clusterApplyLB(cl, jobList, .rescoreJob, fn = fn, fnArgs = fnArgs)
  }
  pSum <- lapply(sims, attr, "spreadProb")
  out <- data.table::rbindlist(lapply(seq_along(sims), function(i)
    data.table::data.table(member = i, sims[[i]])), use.names = TRUE)
  data.table::setattr(out, "spreadProb", pSum)
  out
}

#' @param sims The output of `simulateFireSizes()`.
#' @param probs Quantiles of fire size to compare.
#' @param largeFire Size in pixels above which a fire counts as large (174 pixels = 1000 ha at
#'   5.76 ha per pixel).
#'
#' @return `scoreFireSizes()`: a one-row `data.table`. `fireBias` and `fireRMSE` are the mean and
#'   root mean square of log10(simulated / observed) per fire (0.3 is a factor of 2); `year*` the
#'   same for annual area burned; `yearIn90pct` the share of years whose observed area lies within
#'   the 5-95% range of simulated years; `simsOver10xObs` and `simsUnderTenthObs` the share of
#'   simulated fires more than 10 times too large or too small; `AD` the Anderson-Darling statistic
#'   between simulated and observed sizes; `qObs_*`/`qSim_*` the quantiles in `probs`.
#' @export
#' @rdname fitDiagnostics
scoreFireSizes <- function(sims, probs = c(0.01, 0.10, 0.25, 0.50, 0.75, 0.90, 0.99),
                           largeFire = 174) {
  d <- data.table::as.data.table(sims)
  nNo <- data.table::uniqueN(d[is.na(sim), list(member, yr)])
  d <- d[!is.na(sim)]
  sp <- function(x, y) stats::cor(x, y, method = "spearman")
  fire <- d[, list(obs = size[1], sim = mean(sim)), by = c("yr", "ids")]
  yrRep <- d[, list(area = sum(sim)), by = c("yr", "member", "rep")]
  yr <- merge(fire[, list(obs = sum(obs)), by = "yr"],
              yrRep[, list(sim = mean(area), lo = stats::quantile(area, 0.05),
                           hi = stats::quantile(area, 0.95)), by = "yr"], by = "yr")
  out <- data.table::data.table(
    years = nrow(yr), fires = nrow(fire), memberYearsNotSimulated = nNo,
    fireBias = mean(log10(fire$sim / fire$obs)), fireRMSE = sqrt(mean(log10(fire$sim / fire$obs)^2)),
    fireSpearman = sp(fire$obs, fire$sim),
    simsOver10xObs = mean(d$sim > 10 * d$size), simsUnderTenthObs = mean(d$sim < d$size / 10),
    yearBias = mean(log10(yr$sim / yr$obs)), yearRMSE = sqrt(mean(log10(yr$sim / yr$obs)^2)),
    yearSpearman = sp(yr$obs, yr$sim), yearIn90pct = mean(yr$obs >= yr$lo & yr$obs <= yr$hi),
    totalAreaRatio = sum(yr$sim) / sum(yr$obs), AD = adStatistic(d$sim, fire$obs),
    shareLarge_obs = mean(fire$obs >= largeFire), shareLarge_sim = mean(d$sim >= largeFire))
  qn <- sprintf("q%02d", round(100 * probs))
  out[, (paste0("qObs_", qn)) := as.list(stats::quantile(fire$obs, probs, names = FALSE))]
  out[, (paste0("qSim_", qn)) := as.list(stats::quantile(d$sim, probs, names = FALSE))]
  out[]
}

#' @param probsP Quantiles of spread probability to report.
#' @return `linkSaturation()`: a `data.table`, one row per member: `pixelYears`, `propAtCeiling`
#'   and the spread-probability quantiles `p_*` (bin midpoints).
#' @export
#' @rdname fitDiagnostics
linkSaturation <- function(sims, probsP = c(0.01, 0.10, 0.25, 0.50, 0.75, 0.90, 0.99)) {
  s <- attr(sims, "spreadProb")
  if (!is.null(s$counts)) s <- list(s) # one summary, from .objfunSpreadFit(returnSims = TRUE)
  data.table::rbindlist(lapply(seq_along(s), function(i) {
    x <- s[[i]]
    cum <- cumsum(x$counts) / x$n
    q <- (vapply(probsP, function(pr) which(cum >= pr)[1], integer(1)) - 0.5) * x$binWidth
    data.table::data.table(member = i, pixelYears = x$n, propAtCeiling = x$atCeiling / x$n,
                           t(stats::setNames(q, sprintf("p_q%02d", round(100 * probsP)))))
  }))
}

## The coefficients of the covariates: every parameter that is not a logistic parameter
covariateCoefs <- function(nms) setdiff(nms, unlist(logisticParamNames))

#' @param lower,upper The bounds the fit used, named as the columns of `pop`.
#' @param zMin The `zPop` at or above which the population is taken to pin a coefficient's sign.
#'
#' @return `coefIdentifiability()`: a `data.table`, one row per covariate coefficient: the
#'   population median `popMedian`; `spread90`, the 5-95% range across members as a share of the
#'   bound width; `atBound`, the share of members within 5% of the width from a bound;
#'   `zPop = |popMedian| / ((q95 - q5) / 3.29)`, the median in units of the population's spread
#'   (a normal's 5-95% range is 3.29 sd); `rSlope`, the correlation across members with
#'   `hillSlope1` (the link uses `hillSlope1 * X b`, so the two trade off by construction); and
#'   `signPinned = zPop >= zMin`. With three runs each of four ELFs, `zPop >= 1` found the
#'   coefficients whose sign every run agreed on with 94% precision and 44% recall.
#' @export
#' @rdname fitDiagnostics
coefIdentifiability <- function(pop, lower, upper, zMin = 1) {
  pop <- as.matrix(pop)
  h <- if ("hillSlope1" %in% colnames(pop)) pop[, "hillSlope1"]
  data.table::rbindlist(lapply(covariateCoefs(colnames(pop)), function(cn) {
    x <- pop[, cn]
    w <- upper[[cn]] - lower[[cn]]
    q <- stats::quantile(x, c(0.05, 0.95), names = FALSE)
    zPop <- abs(stats::median(x)) / (diff(q) / 3.29)
    data.table::data.table(
      coef = cn, popMedian = stats::median(x), spread90 = diff(q) / w,
      atBound = mean(x <= lower[[cn]] + 0.05 * w | x >= upper[[cn]] - 0.05 * w),
      zPop = zPop,
      rSlope = if (!is.null(h) && stats::sd(x) > 0) stats::cor(x, h) else NA_real_,
      signPinned = !is.na(zPop) && zPop >= zMin) # NaN: every member 0, so no sign at all
  }))
}

#' @param best Named numeric: the parameter set to profile around, normally the best member by
#'   replicated mean ([bestByReplicatedMean()]).
#' @param reps Evaluations of each profile point.
#' @param coefs The coefficients to profile; all covariate coefficients by default.
#' @param probsPop Population quantiles to evaluate each coefficient at, besides 0.
#'
#' @return `profileCoefficients()`: a `data.table`, one row per coefficient x profile point: `coef`,
#'   `at` (its value there; `isZero` marks the drop test), `mean` of the objective over `reps`, and
#'   `delta` and `deltaSE`: the mean paired difference from `best` (same seeds) and its standard
#'   error. A positive `delta` means the point fits worse than `best`. A reference row
#'   (`coef = ""`) holds `best` itself.
#' @export
#' @rdname fitDiagnostics
profileCoefficients <- function(best, pop, fn = .objfunSpreadFit, reps = 10L, cl = NULL, seed = 1L,
                                fnArgs = list(), coefs = covariateCoefs(names(best)),
                                probsPop = c(0.05, 0.25, 0.50, 0.75, 0.95)) {
  pop <- as.matrix(pop)
  pts <- data.table::rbindlist(c(
    list(data.table::data.table(coef = "", at = NA_real_, isZero = FALSE)),
    lapply(coefs, function(cn) data.table::data.table(
      coef = cn, at = c(0, stats::quantile(pop[, cn], probsPop, names = FALSE)),
      isZero = c(TRUE, rep(FALSE, length(probsPop)))))))
  jobs <- data.table::CJ(point = seq_len(nrow(pts)), rep = seq_len(reps))
  jobList <- lapply(seq_len(nrow(jobs)), function(i) {
    par <- best
    pt <- pts[jobs$point[i]]
    if (nzchar(pt$coef)) par[[pt$coef]] <- pt$at
    list(par = par, seed = as.integer(seed) + jobs$rep[i])
  })
  vals <- if (is.null(cl)) {
    lapply(jobList, .rescoreJob, fn = fn, fnArgs = fnArgs)
  } else {
    parallel::clusterApplyLB(cl, jobList, .rescoreJob, fn = fn, fnArgs = fnArgs)
  }
  jobs[, value := as.numeric(unlist(vals))]
  ref <- jobs[point == 1L, list(rep, refValue = value)]
  jobs <- ref[jobs, on = "rep"]
  res <- jobs[, list(mean = mean(value), delta = mean(value - refValue),
                     deltaSE = stats::sd(value - refValue) / sqrt(.N)), by = "point"]
  cbind(pts, res[, -"point"])
}

#' @param ident The output of `coefIdentifiability()`.
#' @param profile The output of `profileCoefficients()`.
#' @param seMult Dropping a coefficient must worsen the objective by more than `seMult` standard
#'   errors.
#'
#' @return `identifiedInIsolation()`: `ident` with `dropDelta`, `dropDeltaSE`, `dropMatters` and
#'   `identified` (`signPinned & dropMatters`) added.
#' @export
#' @rdname fitDiagnostics
identifiedInIsolation <- function(ident, profile, seMult = 2) {
  drop <- data.table::as.data.table(profile)[isZero == TRUE,
                                             list(coef, dropDelta = delta, dropDeltaSE = deltaSE)]
  out <- drop[data.table::as.data.table(ident), on = "coef"]
  out[, dropMatters := dropDelta > seMult * dropDeltaSE]
  out[, identified := signPinned & dropMatters %in% TRUE]
  out[]
}

#' @param DE The output of [runDEoptim()]: a list of DEoptim results, one per chunk of generations.
#' @return `fitConvergence()`: a `data.table`, one row per chunk: `generation` (its last),
#'   `bestval` (the best value so far) and the median and standard error of the median of the
#'   population's values (fail values of `1e6` and above left out).
#' @export
#' @rdname fitDiagnostics
fitConvergence <- function(DE) {
  gens <- cumsum(vapply(DE, function(d) length(d$member$bestvalit), integer(1)))
  data.table::rbindlist(lapply(seq_along(DE), function(i) {
    v <- as.numeric(DE[[i]]$member$popval)
    v <- v[is.finite(v) & v < 1e6]
    data.table::data.table(generation = gens[i], bestval = min(DE[[i]]$member$bestvalit),
                           popMedian = if (length(v)) stats::median(v) else NA_real_,
                           popMedianSE = if (length(v) > 1) 1.2533 * stats::sd(v) / sqrt(length(v))
                                         else NA_real_)
  }))
}

utils::globalVariables(c("area", "dropDelta", "dropDeltaSE", "dropMatters", "identified", "ids",
                         "isZero", "obs", "point", "refValue", "signPinned", "sim", "size", "delta",
                         "deltaSE", "yr"))
