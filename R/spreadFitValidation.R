#' Compare a spread fit with the data it was fitted to
#'
#' A response curve of `logit(p)` against a covariate cannot show misfit: its points are the model's
#' own predictions, so they sit on the fitted curve by construction. What can show misfit is the
#' share of pixels that actually burned against the share that burned in the fit's own simulations,
#' binned by each covariate.
#'
#' * `spreadFitValidationData()` runs the objective's simulations for the fitted parameters,
#'   recording which pixels burned (`.objfunSpreadFit(returnBurned = TRUE)`), and returns one row per
#'   spreadable pixel-year in the buffers of the fitted fires.
#' * `plotSpreadFitValidation()` plots, per covariate, the observed and the simulated burned share
#'   by quantile bin of that covariate.
#' * `plotSpreadFitResponse()` plots the fitted response curves (the model, not the data).
#'
#' Both plot functions take the data as their first argument and return a `ggplot`, so they can be
#' given to [SpaDES.core::Plots()] as `fn`.
#'
#' @param par The fitted parameters, as given to [.objfunSpreadFit()] (without `hillSlope1` and
#'   `inflectionPoint1`, which are fixed at 1; `yearSpreadSD` last if fitted).
#' @inheritParams .objfunSpreadFit
#' @param Nreps Simulations per fire; use the fit's `Nreps`.
#' @param seed If not `NULL`, the simulations run with this seed and the session's random number
#'   stream is restored afterwards.
#' @param ... Further arguments to [.objfunSpreadFit()], e.g. `escapeSizeHa`, `jumpTries`,
#'   `jumpMeanDist`: give it what the fit was given. The objective's gates on implausible
#'   landscapes are lifted (`lanscape1stQuantileThresh = Inf`) unless given here.
#'
#' @return `spreadFitValidationData()`: a `data.table`, one row per spreadable pixel-year inside
#'   the buffers of the fires the objective fits (`minFireSize`, `escapeSizeHa`): `year` (the name
#'   of the fire year), `pixelID`, `observed` (1 if the pixel is inside an observed fire, else 0),
#'   `simulated` (the share of the `Nreps` simulations in which it burned), `p` (its fitted spread
#'   probability, without the per-year random effect), and one column per covariate of the formula
#'   with its raw value (the stored integer / 1000, before rescaling and mutual exclusivity).
#'   `attr(, "spreadFitValidation")` holds what the plots need: the parameters, `covMinMax`,
#'   `covCentre`, `lowerSpreadProb`, `link`, `Nreps`, the names of the covariates and of the
#'   annual (climate) ones, and `yearsNotSimulated` -- years the objective declined to simulate,
#'   which have no rows.
#'
#' @export
#' @rdname spreadFitValidation
spreadFitValidationData <- function(par, landscape, annualDTx1000, nonAnnualDTx1000, formulaToFit,
                                    historicalFires, fireBufferedListDT, covMinMax = NULL,
                                    mutuallyExclusive = list("youngAge" = c("class", "nf")),
                                    covCentre = NULL, Nreps = 10, lowerSpreadProb = spreadProbFloor,
                                    maxFireSpread = spreadProbCeiling, link = NULL, fitYearSpreadSD = NULL,
                                    seed = NULL, ...) {
  if (!is.null(seed)) withr::local_seed(seed)
  objArgs <- utils::modifyList(
    list(lanscape1stQuantileThresh = Inf, doAssertions = FALSE, verbose = 0),
    list(...))
  sims <- do.call(.objfunSpreadFit, c(list(
    par = par, landscape = landscape, annualDTx1000 = annualDTx1000,
    nonAnnualDTx1000 = nonAnnualDTx1000, formulaToFit = formulaToFit,
    historicalFires = historicalFires, fireBufferedListDT = fireBufferedListDT,
    covMinMax = covMinMax, mutuallyExclusive = mutuallyExclusive, covCentre = covCentre,
    Nreps = Nreps, lowerSpreadProb = lowerSpreadProb, maxFireSpread = maxFireSpread, link = link,
    fitYearSpreadSD = fitYearSpreadSD, returnBurned = TRUE), objArgs))
  burned <- attr(sims, "burned")

  ## the parameters as the objective uses them
  if (is.null(fitYearSpreadSD)) fitYearSpreadSD <- yearSpreadSDTxt %in% names(par)
  if (isTRUE(fitYearSpreadSD)) par <- par[-length(par)]
  par <- fixLogisticPars(par)
  colsToUse <- spreadDesignCols(formulaToFit)
  covCols <- spreadCovCols(colsToUse)
  pp <- paramsSeparate(par, length(colsToUse))

  lapply(nonAnnualDTx1000, setDT)
  indexNonAnnual <- nonAnnualIndex(nonAnnualDTx1000)
  simYears <- names(Filter(Negate(is.null), burned))
  out <- rbindlist(lapply(simYears, function(y) {
    ann <- annualDTx1000[[y]]
    used <- spreadProbFromIntegerCovs(NULL, ann, nonAnnualDTx1000, indexNonAnnual, y, covMinMax,
                                      mutuallyExclusive, colsToUse, FALSE, pp$logisticPars,
                                      pp$covPars, maxFireSpread, lowerSpreadProb, covCentre = covCentre)
    p <- as.numeric(logisticAll(pp$logisticPars, as.matrix(used[, ..colsToUse]), pp$covPars,
                                lowerSpreadProb, link = link))
    raw <- spreadProbFromIntegerCovs(NULL, ann, nonAnnualDTx1000, indexNonAnnual, y, NULL, NULL,
                                     colsToUse, FALSE, pp$logisticPars, pp$covPars, maxFireSpread,
                                     lowerSpreadProb)
    ## the buffers of the fires the objective simulated this year; a pixel in two buffers once
    fb <- as.data.table(fireBufferedListDT[[y]])
    fb <- fb[fb$ids %in% sims$ids[sims$yr == y]]
    obs <- fb[, list(observed = as.integer(any(buffer == 1L))), by = "pixelID"]
    keep <- match(obs$pixelID, used$pixelID)
    obs <- obs[!is.na(keep)]
    keep <- keep[!is.na(keep)]
    nBurned <- burned[[y]][, list(n = .N), by = "pixelID"]
    d <- data.table(year = y, pixelID = obs$pixelID, observed = obs$observed,
                    simulated = 0, p = p[keep])
    m <- match(d$pixelID, nBurned$pixelID)
    set(d, which(!is.na(m)), "simulated", nBurned$n[m[!is.na(m)]] / Nreps)
    for (cn in covCols) set(d, NULL, cn, raw[[cn]][keep])
    d
  }))
  data.table::setattr(out, "spreadFitValidation", list(
    logisticPars = pp$logisticPars, covPars = pp$covPars,
    covMinMax = if (!is.null(covMinMax)) as.list(covMinMax), covCentre = covCentre,
    lowerSpreadProb = lowerSpreadProb, link = link, Nreps = Nreps, covariates = covCols,
    intercept = spreadInterceptTxt %in% colsToUse,
    climateCols = intersect(covCols, names(annualDTx1000[[1]])),
    yearsNotSimulated = setdiff(names(burned), simYears)))
  out
}

## what the plots need to know about each covariate
validationInfo <- function(d) {
  info <- attr(d, "spreadFitValidation")
  if (is.null(info)) stop("`d` must come from spreadFitValidationData()")
  info
}

## "indicator": two values or fewer; "climate": an annual covariate; "other": the rest (fuels and
## the like), where the lowest value (covMinMax's minimum, else 0) means absent
covariateKind <- function(d, cn, info) {
  if (length(unique(d[[cn]])) <= 2L) "indicator" else if (cn %in% info$climateCols) "climate" else "other"
}

absentValue <- function(cn, info) {
  if (!is.null(info$covMinMax[[cn]])) info$covMinMax[[cn]][1] else 0
}

## the pixels to show for a covariate: all of them, except absent ones for "other" covariates
presentRows <- function(d, cn, kind, info) {
  if (identical(kind, "other")) d[[cn]] > absentValue(cn, info) + 1e-3 else rep(TRUE, NROW(d))
}

validationTheme <- function() {
  ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(panel.grid.minor = ggplot2::element_blank(), plot.title.position = "plot",
                   legend.position = "bottom",
                   strip.text = ggplot2::element_text(face = "bold", hjust = 0),
                   plot.subtitle = ggplot2::element_text(size = 9),
                   plot.background = ggplot2::element_rect(fill = "white", colour = NA))
}

#' @param d The output of `spreadFitValidationData()`.
#' @param nBins Number of quantile bins per covariate.
#' @param covariates The covariates to plot; all of the formula's by default.
#' @param title Plot title; a default is used when `NULL`.
#'
#' @return `plotSpreadFitValidation()`: a `ggplot`, one panel per covariate. Pixel-years are binned
#'   by quantiles of the covariate: all of them for climate covariates, only those where it is
#'   present (above its `covMinMax` minimum, else above 0) for the others, and one bin per value
#'   for covariates with two values or fewer. Each bin shows the observed burned share and the
#'   mean simulated burned share, at the bin's mean covariate value, with its number of
#'   pixel-years along the bottom.
#' @export
#' @rdname spreadFitValidation
plotSpreadFitValidation <- function(d, nBins = 10L, covariates = NULL, title = NULL) {
  info <- validationInfo(d)
  covs <- if (is.null(covariates)) info$covariates else covariates
  binned <- rbindlist(lapply(covs, function(cn) {
    kind <- covariateKind(d, cn, info)
    keep <- presentRows(d, cn, kind, info)
    x <- d[[cn]][keep]
    if (!length(x)) return(NULL)
    ux <- sort(unique(x))
    bin <- if (length(ux) <= nBins) {
      match(x, ux)
    } else {
      br <- unique(stats::quantile(x, seq(0, 1, length.out = nBins + 1L), names = FALSE))
      findInterval(x, br, rightmost.closed = TRUE, all.inside = TRUE)
    }
    b <- data.table(bin = bin, x = x, observed = d$observed[keep], simulated = d$simulated[keep])
    b <- b[, list(value = mean(x), observed = mean(observed), simulated = mean(simulated), n = .N),
           by = "bin"]
    set(b, NULL, "covariate", cn)
    b
  }))
  set(binned, NULL, "covariate", factor(binned$covariate, levels = covs))
  series <- c("observed", "simulated by the fit")
  long <- rbind(binned[, list(covariate, value, share = observed, series = series[1])],
                binned[, list(covariate, value, share = simulated, series = series[2])])
  set(long, NULL, "series", factor(long$series, levels = series))
  set(binned, NULL, "label", formatC(binned$n, format = "d", big.mark = ","))
  if (is.null(title)) title <- "Burned share of pixel-years: observed against simulated by the fit"
  sub <- paste0(
    "Pixel-years in the buffers of the fitted fires (", formatC(NROW(d), format = "d", big.mark = ","),
    "), binned by quantiles of each covariate: only where it is present, except for climate and ",
    "0/1 covariates.\nSimulated: the mean share of ", info$Nreps, " simulations of each fire, from the ",
    "fitted parameters, in which the pixel burned. Numbers: pixel-years per bin.")
  ggplot2::ggplot(long, ggplot2::aes(value, share, colour = series)) +
    ggplot2::geom_line(linewidth = 0.8) +
    ggplot2::geom_point(size = 1.8) +
    ggplot2::geom_text(data = binned, ggplot2::aes(x = value, y = -Inf, label = label),
                       inherit.aes = FALSE, angle = 90, hjust = -0.1, size = 2.3, colour = "grey45") +
    ggplot2::facet_wrap(~covariate, scales = "free") +
    ggplot2::expand_limits(y = 0) +
    ggplot2::scale_colour_manual(values = stats::setNames(c("#1b1b1b", "#eb6834"), series), name = NULL) +
    ggplot2::labs(x = "covariate value (bin mean, raw units)", y = "share of pixel-years burned",
                  title = title, subtitle = sub) +
    validationTheme()
}

#' @param bins Number of hexagons (or squares, without \pkg{hexbin}) across each panel.
#'
#' @return `plotSpreadFitResponse()`: a `ggplot`, one panel per covariate with more than two values.
#'   The line is `logit(p)` the fitted parameters give a "pure stand": that covariate alone, every
#'   other covariate absent (0 after rescaling: other fuels, youngAge, non-forest), and climate at
#'   its median. Behind it, the pixel-years' own predicted `logit(p)` (only where the covariate is
#'   present, except for climate), up to its 99.5th percentile. Those follow the curve by
#'   construction: the figure shows the model's response, not how well it fits.
#' @export
#' @rdname spreadFitValidation
plotSpreadFitResponse <- function(d, covariates = NULL, bins = 40L, title = NULL) {
  info <- validationInfo(d)
  covs <- if (is.null(covariates)) info$covariates else covariates
  kinds <- vapply(covs, covariateKind, character(1), d = d, info = info)
  covs <- covs[kinds != "indicator"]
  kinds <- kinds[kinds != "indicator"]
  allCovs <- info$covariates
  toUsed <- function(cn, v) {
    mm <- info$covMinMax[[cn]]
    if (is.null(mm)) v else (v - mm[1]) / (mm[2] - mm[1])
  }
  base <- stats::setNames(numeric(length(allCovs)), allCovs)
  for (cn in info$climateCols) base[[cn]] <- stats::median(toUsed(cn, d[[cn]]))
  centre <- stats::setNames(numeric(length(allCovs)), allCovs)
  for (cn in intersect(names(info$covCentre), allCovs)) centre[[cn]] <- info$covCentre[[cn]]
  logitP <- function(mat) {
    ## covPars has the intercept's coefficient first when the fit has one; it multiplies a column of 1s
    if (isTRUE(info$intercept)) mat <- cbind(1, mat)
    stats::qlogis(as.numeric(logisticAll(info$logisticPars, mat, info$covPars, info$lowerSpreadProb,
                                         link = info$link)))
  }
  pts <- list()
  cur <- list()
  for (i in seq_along(covs)) {
    cn <- covs[i]
    keep <- presentRows(d, cn, kinds[i], info)
    x <- d[[cn]][keep]
    if (!length(x)) next
    xmax <- stats::quantile(x, 0.995, names = FALSE)
    lo <- if (kinds[i] == "climate") min(x) else absentValue(cn, info)
    inRange <- x <= xmax
    pts[[cn]] <- data.table(covariate = cn, value = x[inRange], logitP = stats::qlogis(d$p[keep][inRange]))
    grid <- seq(lo, xmax, length.out = 200L)
    mat <- matrix(base, nrow = length(grid), ncol = length(allCovs), byrow = TRUE,
                  dimnames = list(NULL, allCovs))
    mat[, cn] <- toUsed(cn, grid)
    mat <- sweep(mat, 2L, centre)
    cur[[cn]] <- data.table(covariate = cn, value = grid, logitP = logitP(mat))
  }
  pts <- rbindlist(pts)
  cur <- rbindlist(cur)
  lv <- intersect(covs, unique(cur$covariate))
  set(pts, NULL, "covariate", factor(pts$covariate, levels = lv))
  set(cur, NULL, "covariate", factor(cur$covariate, levels = lv))
  limits <- stats::qlogis(c(info$lowerSpreadProb, info$logisticPars[[1]]))
  geomBins <- if (requireNamespace("hexbin", quietly = TRUE)) ggplot2::geom_hex else ggplot2::geom_bin_2d
  if (is.null(title)) title <- "Model response curves: what the fitted model predicts, not how well it fits"
  sub <- paste0(
    "Line: logit(spread probability) the fitted parameters give a pure stand -- this covariate alone, ",
    "every other covariate absent (other fuels, youngAge, non-forest),\nclimate at its median. ",
    "Grey: the pixel-years' own predicted logit(p), which follows the model by construction and ",
    "cannot show misfit;\ncompare observed with simulated burning for that. ",
    "Dashed: the link's floor and ceiling.")
  ggplot2::ggplot(pts, ggplot2::aes(value, logitP)) +
    geomBins(bins = bins) +
    ggplot2::scale_fill_gradient(low = "#e6e6e3", high = "#3d3d3a", trans = "log10",
                                 name = "pixel-years (log scale)") +
    ggplot2::geom_hline(yintercept = limits, linetype = "dashed", colour = "#52514e", linewidth = 0.4) +
    ggplot2::geom_line(data = cur, colour = "#0b0b0b", linewidth = 1.1) +
    ggplot2::facet_wrap(~covariate, scales = "free_x") +
    ggplot2::labs(x = "covariate value (raw units)", y = "logit(spread probability)",
                  title = title, subtitle = sub) +
    validationTheme()
}

utils::globalVariables(c("covariate", "label", "logitP", "observed", "series", "share", "simulated"))
