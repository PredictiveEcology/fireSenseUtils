#' Multiple-parameter versions of logistic functions
#'
#' Logistic functions using 2, 3, 4, or 5 parameters; `logisticAll` dispatches
#' to one of them based on `length(logisticPars)`.
#'
#' The general form is the 5-parameter sigmoid (Richards' curve):
#' `par1 + (par2 - par1) / (1 + (exp(x) / par3)^(-par4))^par5`.
#' The 4-, 3-, and 2-parameter variants successively fix the lower asymptote
#' (`par1`), the inflection-point scale (`par3 = 1`), and the asymmetry
#' factor (`par5 = 0.5`) to reduce the parameter count.
#'
#' @param x Numeric vector. Linear predictor (e.g. `mat %*% covPars`) at which
#'   to evaluate the logistic.
#' @param par Numeric vector of logistic parameters whose length matches the
#'   variant: 2 for `logistic2p`, 3 for `logistic3p`, etc. See
#'   [logisticParamNames] for the meaning of each position.
#' @param logisticPars Numeric vector of logistic parameters. Unless `link` names the form, its
#'   length selects the variant `logisticAll` dispatches to -- except that a vector with an element
#'   named `upperTail1` is always `logistic3pUpper`, whatever its length.
#' @param mat Numeric matrix of covariate values, one row per observation.
#' @param covPars Numeric vector of covariate coefficients (same length as
#'   `ncol(mat)`); `mat %*% covPars` forms the linear predictor.
#' @param lowerSpreadProb Numeric scalar in `[0, 1]`. The lower asymptote
#'   (`par1`) used for the 2- and 3-parameter forms.
#' @param link `NULL` (the default) to choose the form from `logisticPars` as above,
#'   `"logistic3pUpper"` to use the 3-parameter form with Stukel's upper tail, or `"powerLaw"`,
#'   `"powerLawSmooth"` or `"richards"` for [powerLawLink()], [powerLawSmoothLink()] or [richardsLink()]
#'   (`logisticPars[1]` is then `p_c`). Pass it explicitly
#'   wherever `logisticPars` may arrive without names, as inside the objective function.
#' @param gamma The `"powerLaw"` and `"powerLawSmooth"` links' scale `gamma` ([powerLawLink()]); not used by the other links
#'   (`"richards"` has its own fixed scale, [richardsLink()]).
#' @param shape The `"powerLaw"` and `"powerLawSmooth"` links' curvature exponent `shape` ([powerLawLink()]); default 1; not used by the other links.
#' @param lowerAsymptote The lower asymptote `L` of the `"powerLawSmooth"` and `"richards"` links ([powerLawSmoothLink()],
#'   [richardsLink()]), [spreadLowerAsymptoteMin] `<= L <= spreadProbFloor` ([checkLinkParams()]); default [spreadProbFloor]; not used by the other links.
#' @param richardsNu The `"richards"` link's `nu` ([richardsLink()]); default [spreadRichardsNu]; not used by the other links.
#'
#' @return Numeric vector of logistic values; same length as `x` (or
#'   `nrow(mat)` for `logisticAll`).
#'
#' @export
#' @rdname logistic
logistic4p <- function(x, par) {
  par[1L] + (par[2L] - par[1L]) / (1 + exp(x)^(-par[3L]))^par[4L]
}

#' @export
#' @rdname logistic
logistic5p <- function(x, par) {
  par[1L] + (par[2L] - par[1L]) / (1 + (exp(x) / par[3L])^(-par[4L]))^par[5L]
}

#' @param par1 Numeric scalar. Lower asymptote (replaces `par[1]` of the
#'   4-parameter form).
#' @export
#' @rdname logistic
logistic3p <- function(x, par, par1 = 0.1) {
  par1 + (par[1L] - par1) / (1 + exp(x)^(-par[2L]))^par[3L]
}

#' @details `logistic3pUpper` is `logistic3p` with Stukel's (1988) generalized-logistic upper tail
#'   applied to the linear predictor: `par1 + (par[1] - par1) / (1 + exp(-H(par[2] * x)))^par[3]`,
#'   where `H` leaves negative values alone and, for `u = par[2] * x >= 0` and `a = par[4]`, is
#'   `(exp(a * u) - 1) / a` for `a > 0`, `u` for `a = 0` and `-log(1 - a * u) / a` for `a < 0`
#'   (Stukel, T. A. (1988), Generalized logistic models, JASA 83:426-431; as `sirt::pgenlogis()`).
#'   `par[4]` changes only how the curve approaches its upper asymptote: negative values slow the
#'   approach, so pixels with a high linear predictor are no longer all pressed against the ceiling;
#'   `par[4] = 0` is `logistic3p` exactly. In `logistic3p` the approach to the ceiling is set by the
#'   slope `par[2]` alone -- `par[3]` shapes only the lower end -- so this is the one parameter that
#'   changes the upper end independently.
#' @export
#' @rdname logistic
logistic3pUpper <- function(x, par, par1 = 0.1) {
  par1 + (par[1L] - par1) / (1 + exp(-upperTail(par[2L] * x, par[4L])))^par[3L]
}

#' Stukel's upper-tail transformation
#'
#' The upper half of the generalized logistic link of Stukel (1988): values below zero are returned
#' unchanged; values `u >= 0` become `(exp(a * u) - 1) / a` (`a > 0`), `u` (`a = 0`) or
#' `-log(1 - a * u) / a` (`a < 0`). Matches the upper branch of `sirt::pgenlogis()`.
#'
#' @param u Numeric vector.
#' @param a Numeric scalar, the upper-tail parameter.
#' @return Numeric vector the length of `u`.
#' @keywords internal
upperTail <- function(u, a) {
  if (a == 0) return(u)
  pos <- u >= 0
  u[pos] <- if (a > 0) (exp(a * u[pos]) - 1) / a else -log(1 - a * u[pos]) / a
  u
}

#' @param par4 Numeric scalar. Asymmetry factor (replaces `par[4]` of the
#'   4-parameter form).
#' @export
#' @rdname logistic
logistic2p <- function(x, par, par1 = 0.1, par4 = 0.5) {
  par1 + (par[1L] - par1) / (1 + exp(x)^(-par[2L]))^par4
}

#' @export
#' @rdname logistic
logisticAll <- function(logisticPars, mat, covPars, lowerSpreadProb, link = NULL, gamma = spreadPowerLawGamma,
                        shape = spreadPowerLawShape, lowerAsymptote = spreadProbFloor, richardsNu = spreadRichardsNu) {
  if (is.null(link) && "upperTail1" %in% names(logisticPars)) link <- "logistic3pUpper"
  if (!is.null(link)) {
    link <- match.arg(link, c("logistic3pUpper", spreadLinearLinks))
    if (link %in% spreadLinearLinks)
      return(linearLinkProb(mat %*% covPars, link, linearLinkPars(link, logisticPars[[1L]], gamma, lowerSpreadProb,
                                                                  shape, lowerAsymptote, richardsNu)))
    return(logistic3pUpper(mat %*% covPars, logisticPars, par1 = lowerSpreadProb))
  }
  if (length(logisticPars) == 4) {
    stop("logistic with 4 parameters not tested yet")
    logistic4p(mat %*% covPars, logisticPars)
  } else if (length(logisticPars) == 3) {
    logistic3p(mat %*% covPars, logisticPars, par1 = lowerSpreadProb)
  } else if (length(logisticPars) == 2) {
    logistic2p(mat %*% covPars, logisticPars, par1 = lowerSpreadProb)
  }
}


#' The power-law spread link
#'
#' `p = pc - (pc - floor) * exp(-sign(x) * abs(x / gamma)^shape)`, from a percolation calibration (log mean fire size is linear in
#' `log(pc - p)` with slope `-gamma`, `pc` about 0.25 and `gamma` about 1.9). `log(pc - p)` falls linearly in
#' the linear predictor `x`, so a year effect added to `x` multiplies a year's mean fire size by about
#' `exp(eps)`; and `p` cannot pass `pc`, so no year percolates. `floor` is `p` at `x = 0`, not a lower asymptote:
#' `x < 0` gives less than `floor`, and a value below 0 is set to 0, so a cell with a low enough linear predictor
#' cannot burn. [powerLawSmoothLink()] has the same curve for `x >= 0` and a lower asymptote below 0. `pc` takes the slot
#' of `maxAsymptote` and is fitted; `gamma` is
#' fixed, not fitted: it is exactly confounded with the scale of the coefficients. Used by [logisticAll()] with
#' `link = "powerLaw"` (so by [spreadProbFromCovs()] and the objective) and by the objective's year effect.
#' `shape` (default 1, the formula `p = pc - (pc - floor) * exp(-x / gamma)`) is fixed, not fitted. With `shape > 1`
#' a moderate `x` gives a lower `p` and a large `x` gets closer to `pc` (a heavier upper tail); with `shape < 1` the
#' reverse. The exact `log(pc - p)` linearity in `x` holds only at `shape = 1`.
#'
#' @param x Linear predictor (a vector or a one-column matrix).
#' @param pc The percolation threshold, `p` as `x` goes to infinity.
#' @param gamma The scale; default [spreadPowerLawGamma].
#' @param shape The curvature exponent; default [spreadPowerLawShape] (1).
#' @param floor `lowerSpreadProb`: `p` at `x = 0` (not a lower asymptote: see above).
#' @return A numeric vector, one value per `x`.
#' @export
powerLawLink <- function(x, pc, gamma = spreadPowerLawGamma, floor = spreadProbFloor,
                         shape = spreadPowerLawShape) {
  x <- as.vector(x)
  pmax(pc - (pc - floor) * exp(-sign(x) * abs(x / gamma)^shape), 0)
}

#' The power-law spread link with a lower asymptote
#'
#' For `x >= 0` exactly [powerLawLink()]: `p = pc - (pc - floor) * exp(-(x / gamma)^shape)`, so `p = floor` at `x = 0`.
#' For `x < 0`, `p = lower + (floor - lower) * exp(-((pc - floor) / (floor - lower)) * (abs(x) / gamma)^shape)`: it
#' falls from `floor` towards `lower` and never below it, so no cell becomes unburnable (the power law's `x < 0` branch
#' falls to 0). The curve is continuous at 0 and has no clamp. With `lower = floor` the `x < 0` side is the constant
#' `floor`. `gamma` is part of the power law's definition and is calibrated as for [powerLawLink()]
#' ([calibratePowerLaw()]).
#'
#' @inheritParams powerLawLink
#' @param lower The lower asymptote `L`, `0 <= L <= floor`; default `floor`. `pc`, `gamma`, `floor`, `shape` and `lower`
#'   may each be one value or one per `x`.
#' @return A numeric vector, one value per `x`.
#' @export
powerLawSmoothLink <- function(x, pc, gamma = spreadPowerLawGamma, floor = spreadProbFloor,
                               shape = spreadPowerLawShape, lower = floor) {
  x <- as.vector(x)
  p <- powerLawLink(pmax(x, 0), pc, gamma, floor, shape)
  neg <- which(x < 0)
  if (length(neg)) {
    ## the parameters may be one per x (a cell's own, as fireSense_burn gives them)
    par <- function(v) rep_len(v, length(x))[neg]
    pcN <- par(pc); fN <- par(floor); lN <- par(lower)
    p[neg] <- ifelse(lN == fN, fN, lN + (fN - lN) * exp(-((pcN - fN) / (fN - lN)) * abs(x[neg] / par(gamma))^par(shape)))
  }
  p
}

#' The Richards (generalized logistic) spread link
#'
#' `p = lower + (pc - lower) / (1 + nu * exp(-B * x))^(1 / nu)`: a curve between the lower asymptote `lower` and `pc`.
#' `nu = 0` is the Gompertz curve, `p = lower + (pc - lower) * exp(-exp(-B * x))`; `nu = 1` is the logistic.
#' `nu` is fixed, not fitted. `B` (default [spreadRichardsB], `1 / spreadPowerLawGamma`) and the location (0) are fixed conventions, because
#' they are confounded with the scale of the coefficients and the intercept. It does not use `gamma` or the calibration
#' of [calibratePowerLaw()].
#'
#' @inheritParams powerLawLink
#' @param lower The lower asymptote `L`, `0 <= L`; default [spreadProbFloor].
#' @param nu The shape `nu >= 0`; default [spreadRichardsNu].
#' @param B The slope; default [spreadRichardsB].
#' @return A numeric vector, one value per `x`.
#' @export
richardsLink <- function(x, pc, lower = spreadProbFloor, nu = spreadRichardsNu, B = spreadRichardsB) {
  e <- exp(-B * as.vector(x))
  nu <- rep_len(nu, length(e))
  ## (1 + nu * e)^(1 / nu) is exp(log1p(nu * e) / nu); at nu = 0 that is exp(e), the Gompertz curve
  lower + (pc - lower) * exp(-ifelse(nu == 0, e, log1p(nu * e) / nu))
}

#' Check the parameters of the links with a lower asymptote
#'
#' The one check of `lowerAsymptote` and `richardsNu`, used by `fireSense_spreadFit` for its parameters and by
#' `fireSense_spreadPredict` for a ledger row. `lower` must be one number from [spreadLowerAsymptoteMin] to [spreadProbFloor] with
#' `"powerLawSmooth"` and `"richards"`; `nu` one number, 0 or more, with `"richards"`. Other links have nothing to check.
#'
#' @param link The spread link.
#' @param lower,nu The `lowerAsymptote` and `richardsNu` to check; `NULL` is not a number.
#' @return `TRUE`, invisibly; an error otherwise.
#' @export
checkLinkParams <- function(link, lower = NULL, nu = NULL) {
  one <- function(v) is.numeric(v) && length(v) == 1L && is.finite(v)
  if (any(link %in% spreadAsymptoteLinks) &&
      !(one(lower) && lower >= spreadLowerAsymptoteMin && lower <= spreadProbFloor))
    stop("lowerAsymptote must be one number from ", spreadLowerAsymptoteMin, " to ", spreadProbFloor, "; got ", format(lower),
         ". It is at least ", spreadLowerAsymptoteMin, ", never 0: a fire under its minSize keeps trying to spread until it ",
         "has reached it, and with a spread probability of about 0 that is an unbounded number of generations.", call. = FALSE)
  if (any(link %in% "richards") && !(one(nu) && nu >= 0))
    stop("richardsNu must be one number, 0 or more; got ", format(nu), call. = FALSE)
  invisible(TRUE)
}

#' The parameters of a linear-predictor link, in the order of `spreadCpp(linkPars =)`
#'
#' The parameter vector of one of [spreadLinearLinks], named as [spreadLinkParNames] says. It is what
#' `SpaDES.tools::spreadCpp(link = link, linkPars = )` takes, and [linearLinkProb()] turns the same vector into `p`, so the
#' R and the C++ link are given the same numbers. `"richards"` uses `lowerAsymptote`, `richardsNu` and its fixed
#' `B = ` [spreadRichardsB], and not `gamma`, `floor` or `shape`.
#'
#' @param link One of [spreadLinearLinks].
#' @param pc The percolation threshold (`maxAsymptote`).
#' @param gamma,floor,shape As in [powerLawLink()].
#' @param lowerAsymptote,richardsNu As in [logisticAll()].
#' @return A named numeric vector.
#' @export
linearLinkPars <- function(link, pc, gamma = spreadPowerLawGamma, floor = spreadProbFloor,
                           shape = spreadPowerLawShape, lowerAsymptote = spreadProbFloor, richardsNu = spreadRichardsNu) {
  all <- c(pc = pc, gamma = gamma, floor = floor, shape = shape, lower = lowerAsymptote, nu = richardsNu,
           B = spreadRichardsB)
  all[spreadLinkParNames[[link]]]
}

#' Spread probability from the linear predictor, for a link of [spreadLinearLinks]
#'
#' The one place the link is chosen: [powerLawLink()], [powerLawSmoothLink()] or [richardsLink()].
#'
#' @param x The linear predictor.
#' @param link One of [spreadLinearLinks].
#' @param pars The link's parameters: a named numeric vector ([linearLinkPars()]), or a matrix with those names as columns
#'   and one row per element of `x`.
#' @return A numeric vector, one value per `x`.
#' @export
linearLinkProb <- function(x, link, pars) {
  par <- function(nm) if (is.matrix(pars)) pars[, nm] else pars[[nm]]
  switch(link,
         powerLaw = powerLawLink(x, par("pc"), par("gamma"), par("floor"), par("shape")),
         powerLawSmooth = powerLawSmoothLink(x, par("pc"), par("gamma"), par("floor"), par("shape"), par("lower")),
         richards = richardsLink(x, par("pc"), par("lower"), par("nu"), par("B")),
         stop("`link` must be one of ", paste(dQuote(spreadLinearLinks, FALSE), collapse = ", "), "; got ", dQuote(link, FALSE)))
}

#' The lower bound the spread-probability gates use
#'
#' With `"powerLawSmooth"` or `"richards"` it is the link's lower asymptote, so a year in which most cells sit near it
#' is not refused for a median below `lowerSpreadProb`; with any other link it is `lowerSpreadProb`, as always.
#'
#' @param link The spread link.
#' @param lowerSpreadProb,lowerAsymptote As in [.objfunSpreadFit()].
#' @return A number.
#' @keywords internal
spreadLinkLower <- function(link, lowerSpreadProb, lowerAsymptote)
  if (any(link %in% spreadAsymptoteLinks)) lowerAsymptote else lowerSpreadProb

#' Which arguments a `link`'s cache key leaves out
#'
#' `powerLawGamma` and `powerLawShape` matter only with `link = "powerLaw"` or `"powerLawSmooth"`; `lowerAsymptote` only
#' with `"powerLawSmooth"` or `"richards"`; `richardsNu` only with `"richards"`. An argument that does not matter for
#' the link is left out of the key, so a fit made without it is found as before.
#'
#' With `link = "powerLaw"` or `"powerLawSmooth"`, `powerLawShape` is also left out when it is the default, 1 (the original
#' link), so a fit made before `powerLawShape` existed is found.
#'
#' @param link The spread link: `NULL`, `"logistic3pUpper"`, `"powerLaw"`, `"powerLawSmooth"` or `"richards"`.
#' @param shape The `powerLawShape` of the call; default [spreadPowerLawShape].
#' @return A character vector of argument names, for `omitArgs` of `Cache()`.
#' @export
powerLawOmitArgs <- function(link, shape = spreadPowerLawShape) {
  c(if (!any(link %in% spreadPowerLawLinks)) c("powerLawGamma", "powerLawShape")
    else if (isTRUE(shape == spreadPowerLawShape)) "powerLawShape",
    if (!any(link %in% spreadAsymptoteLinks)) "lowerAsymptote",
    if (!any(link %in% "richards")) "richardsNu") |> as.character()
}

#' Logistic parameter names
#'
#' @export
#' @return A named list of length 4, with "2p", "3p", "4p", "5p" as names representing the
#'   2-parameter etc. logistic curve.
logisticParamNames <- list("2p" = c("maxAsymptote", "hillSlope1"),
                           "3p" = c("maxAsymptote", "hillSlope1", "inflectionPoint1"),
                           "3pUpper" = c("maxAsymptote", "hillSlope1", "inflectionPoint1", "upperTail1"),
                           "4p" = c("minAsymptote", "maxAsymptote", "inflectionPoint1", "hillSlope1"),
                           "5p" = c("minAsymptote", "maxAsymptote", "inflectionPoint1", "hillSlope1", "asymmetryFactor"))
# logisticAll <- function(logisticPars, covsDT, mat, covPars, lowerSpreadProb) {
#   if (length(logisticPars) == 4) {
#     stop("logistic with 4 parameters not tested yet")
#     set(covsDT, NULL, "spreadProb", logistic4p(mat %*% covPars, logisticPars))
#   } else if (length(logisticPars) == 3) {
#     set(covsDT, NULL, "spreadProb", logistic3p(mat %*% covPars, logisticPars, par1 = lowerSpreadProb))
#   } else if (length(logisticPars) == 2) {
#     set(covsDT, NULL, "spreadProb", logistic2p(mat %*% covPars, logisticPars, par1 = lowerSpreadProb, par4 = 0.5))
#   }
# }

#' Replace `NA`s in a `data.table` with zeros
#'
#' Modifies `DT` in place: for each column in `colsToUse`, every `NA` is
#' overwritten with `0` via [data.table::set()]. Useful when a covariate
#' table is later passed to a model that does not tolerate `NA`s.
#'
#' @param DT A `data.table` to be modified by reference.
#' @param colsToUse Character vector of column names in `DT` to process.
#'   `NULL` (default) means every column.
#'
#' @return `DT`, invisibly modified by reference; columns are unchanged
#'   where no `NA`s were present.
#'
#' @export
#' @importFrom data.table set
dtReplaceNAwith0 <- function(DT, colsToUse = NULL) {
  if (is.null(colsToUse)) {
    colsToUse <- names(DT)
  }
  for (i in colsToUse) {
    nas <- which(is.na(DT[[i]]))
    if (length(nas)) {
      set(DT, nas, i, 0)
    }
  }
  DT
}

#' Convert list of annual `SpatRaster` to `data.table`
#'
#' @param x `RasterStack` or list of rasters to convert to `data.table`
#'   and multiply by 1000 to save space.
#' @param whNotNA Pixel indexes that should go through this process (i.e. not NA)
#' @param ... Not currently used
#'
#' @return `data.table` of the `SpatRaster` or the list
#' @export
#' @rdname annualStackToDTx1000
#'
#' @examples
#' withr::local_package("raster")
#'
#' r1 <- raster(extent(0, 10, 0, 10), vals = 1:100)
#' r2 <- raster(extent(0, 10, 0, 10), vals = 100:1)
#' r3 <- raster(extent(0, 10, 0, 10), vals = 200:101)
#' r4 <- raster(extent(0, 10, 0, 10), vals = 300:201)
#'
#' # list of Rasters
#' lRast <- list(r1, r2, r3)
#' lRast[[1]][5] <- NA
#' whNotNA <- setdiff(1:ncell(r1), 5)
#'
#' # unnamed -- should error
#' try(out1 <- annualStackToDTx1000(lRast, whNotNA))
#'
#' # named
#' names(lRast) <- c("OneToHun", "HunToOne", "TwoHunToOneHun")
#' out1 <- annualStackToDTx1000(lRast, whNotNA)
#'
#' # RasterStack
#' out2 <- annualStackToDTx1000(raster::stack(lRast), whNotNA)
#'
#' # List of RasterStacks
#' s1 <- raster::stack(r1, r2)
#' names(s1) <- names(lRast)[1:2]
#' s2 <- raster::stack(r4, r3)
#' names(s2) <- c(names(lRast)[3], "ThreeHunToTwoHun")
#' out3 <- annualStackToDTx1000(list(s1 = s1, s2 = s2), whNotNA) ## named list required
#'
#' # With duplicated names -- to remove duplicates;
#' #  actually, this doesn't make sense: RasterStack can't have duplicated names
#' names(lRast) <- c("OneToHun", "OneToHun", "TwoHunToOneHun")
#' out4 <- annualStackToDTx1000(raster::stack(lRast), whNotNA)
#'
#' ## cleanup
#' withr::deferred_run()
annualStackToDTx1000 <- function(x, whNotNA, ...) {
  UseMethod("annualStackToDTx1000")
}

#' @export
#' @importFrom data.table as.data.table set
#' @importFrom LandR asInteger
#' @importFrom terra values
#' @rdname annualStackToDTx1000
annualStackToDTx1000.SpatRaster <- function(x, whNotNA, ...) {
  layDT <- as.data.table(values(x))[whNotNA]
  layDT <- dtReplaceNAwith0(layDT)
  set(layDT, NULL, 1L, asInteger(layDT[[1L]] * 1000))
  names(layDT) <- names(x)
  message("Layer ", names(layDT), " converted to data.table")
  layDT
}

#' @export
#' @rdname annualStackToDTx1000
annualStackToDTx1000.Raster <- function(x, whNotNA, ...) {
  annualStackToDTx1000(rast(x), whNotNA, ...)
}

#' @export
#' @importFrom data.table as.data.table
#' @rdname annualStackToDTx1000
annualStackToDTx1000.list <- function(x, whNotNA, ...) {
  # check for names
  # check for rasters
  if (is.null(names(x))) {
    stop("x must be a named list (or stack)")
  }
  out <- lapply(x, whNotNA = whNotNA, annualStackToDTx1000, ...)
  rastersDT <- as.data.table(out)
  # rastersDT <- lapply(names(x), whNotNA = whNotNA, function(x, whNotNA) {
  #   if (any(is(x, "list"), is(x, "RasterStack"))) {
  #     lay <- x[[x]]
  #     } else {
  #       stop("x must be either a list or a RasterStack")
  #     }
  #   layDT <- as.data.table(lay[])[whNotNA]
  #   layDT <- dtReplaceNAwith0(layDT)
  #   names(layDT) <- names(lay)
  #   message("Layer ", names(layDT), " converted to data.table")
  #   return(layDT)
  # })
  # if (is(x, "RasterStack")) { # Should be the prediction, raster stack
  #   rastersDT <- cbindFromList(rastersDT)
  #   rastersDT[ , (names(rastersDT)) := lapply(X = .SD, FUN = function(column){
  #     column <- asInteger(column*1000)
  #     return(column)
  #   }), .SDcols = names(rastersDT)]
  # } else {
  #   lapply(rastersDT, function(x) { # Should be the fitting, list of years. Doesn't work for stack
  #     for (col in colnames(x)) {
  #       set(x, NULL, col, asInteger(x[[col]]*1000))
  #       message("Layer ", col, " converted to integer")
  #     }
  #   })
  # }
  return(rastersDT)
}

#' Generate random beta variates between 2 values and a mean
#'
#' @inheritParams stats::Beta
#' @param shape1 non-negative parameter of the Beta distribution.
#' @param shape2 If provided, passed to [stats::rbeta()]. If not, `m` must be.
#' @param l scalar numeric for the lower bound.
#' @param u scalar numeric for the upper bound.
#' @param m scalar numeric for the mean.
#'
#' @export
#' @importFrom stats rbeta
#' @seealso [stats::rbeta]
rbetaBetween <- function(n, l, u, m, shape1, shape2 = NULL) {
  if (is.null(shape2)) {
    m1 <- ((1) / (u - l) * (m - l))
    shape2 <- (shape1 - shape1 * m1) / m1
  }
  out <- rbeta(n, shape1, shape2)
  out * (u - l) + (l)
}

#' Split the character vector of parameters into `covPars` and `logisticPars`
#'
#' [DEoptim::DEoptim] does not differentiate between the logistic parameters and the covariates.
#' This splits the vector into the correct components.
#' The split is based on the number of covariates.
#' Therefore the number of logistic parameters is deduced from `length(pars) - parsModel`.
#'
#' @param par Numeric vector of all parameters. The covariate parameters must be the
#'   second group.
#' @param parsModel Integer. The number of covariates.
#' @return list of 2 numeric vectors `covPars` and `logisticPars`, representing the
#'   parameters for the covariates and the logistic equation, respectively.
#' @export
paramsSeparate <- function(par, parsModel) {
  covPars <- tail(x = par, n = parsModel)
  logisticPars <- head(x = par, n = length(par) - parsModel)
  list(covPars = covPars, logisticPars = logisticPars)
}

#' The columns of the spread model's design, in the order of its coefficients
#'
#' The terms of the formula, preceded by `spreadInterceptTxt` when the formula has an intercept
#' (`~ 1 + a + b`; the module formulas are `~ 0 + ...` and have none). The intercept is a column of 1s that
#' `spreadProbFromIntegerCovs()` adds after rescaling and centring, so it is multiplied by its
#' coefficient like every covariate, and `paramsSeparate(par, length(spreadDesignCols(f)))` finds it as
#' the first covariate coefficient. [spreadCovCols()] drops it.
#'
#' @param formulaToFit The spread formula, a character string or a formula; `NULL` gives no columns.
#' @return A character vector.
#' @export
spreadDesignCols <- function(formulaToFit) {
  if (is.null(formulaToFit)) return(character(0))
  tt <- terms(as.formula(formulaToFit, env = .GlobalEnv))
  c(if (attr(tt, "intercept") == 1L) spreadInterceptTxt, attr(tt, "term.labels"))
}

#' @param colsToUse A character vector from `spreadDesignCols()`.
#' @return `spreadCovCols()`: `colsToUse` without the intercept, i.e. the columns that are covariates
#'   (rescaled, centred, and in `covMinMax`).
#' @export
#' @rdname spreadDesignCols
spreadCovCols <- function(colsToUse) setdiff(colsToUse, spreadInterceptTxt)

#' Names of the arguments that are `NULL`, to give to `Cache(omitArgs =)`
#'
#' An argument added to a cached call changes its key even when it is `NULL`. Omitting it while it is
#' `NULL` keeps the key of a call made before the argument existed, so the cached result is still found.
#'
#' @param ... Named arguments.
#' @return A character vector: the names of those that are `NULL`.
#' @export
omitNullArgs <- function(...) names(Filter(is.null, list(...)))

## hillSlope1 and inflectionPoint1 (the logistic link's slope and Richards exponent) are fixed at 1,
## not fitted. With the linear predictor x = covariates %*% beta, hillSlope1 enters
## logistic3p()/logistic3pUpper() only as hillSlope1 * x, so scaling every covariate coefficient by k
## and dividing hillSlope1 by k leaves every prediction unchanged -- it is not identifiable, and letting
## DEoptim fit it let every coefficient drift along that ridge. inflectionPoint1 is not a location but
## the exponent in (1 + exp(-x))^-inflectionPoint1: in the 2026-10-04 fits its best values were bimodal
## (at the lower bound or at 2-4), it spanned the bounds within a final population at no cost in the
## objective, and the coefficients explained a median 50% of its variance. With both at 1 the link is
## the plain logistic par1 + (maxAsymptote - par1) * plogis(x).
## `par` (as DEoptim, or a caller mimicking it, supplies it) is always the logistic parameters
## (maxAsymptote first, upperTail1 next if used) followed by the covariates, so the fixed values are
## inserted right after maxAsymptote, in the order of logisticParamNames[["3pUpper"]], before
## `.objfunSpreadFit()` splits and evaluates it. A named `par` that already has one of them (a stored
## full parameter set from before they were fixed) is not given a second.
#' The logistic parameters fixed at 1, not fitted
#'
#' `hillSlope1` and `inflectionPoint1`, in the position they take in the link, right after
#' `maxAsymptote`. `DEoptim`'s `par`, `lower` and `upper` do not contain them; the objective, the
#' validation and the gates insert them (`fixLogisticPars()`), and `fireSense_spreadFit` restores them
#' in the ledger's parameter sets from this vector.
#'
#' @format A named numeric vector.
#' @export
fixedLogisticPars <- c(hillSlope1 = 1, inflectionPoint1 = 1)

fixLogisticPars <- function(par) {
  missingFixed <- if (is.null(names(par))) fixedLogisticPars
                  else fixedLogisticPars[!names(fixedLogisticPars) %in% names(par)]
  append(par, missingFixed, after = 1L)
}

#' The spread-probability ceiling and floor
#'
#' One value each, used wherever the spread model needs them. `spreadProbCeiling` is the upper bound
#' `fireSense_spreadFit` gives `maxAsymptote`, and the default `maxFireSpread` of the objective, the
#' gates and the validation: a fitted `maxAsymptote` above it draws a warning. How a year's random
#' effect can take a pixel past it is described at `maxFireSpread` in [.objfunSpreadFit()].
#' `spreadProbFloor` is the default `lowerSpreadProb` in the objective, the
#' gates, the validation and `fireSense_spreadPredict`. It is the lower asymptote of the logistic links. With
#' `"powerLaw"` it is `p` at `x = 0` and there is no lower asymptote: `p` falls below it for `x < 0` and is set to 0
#' ([powerLawLink()]). With `"powerLawSmooth"` it is `p` at `x = 0` and the default `lowerAsymptote`, the lower
#' asymptote, which may be anywhere from `spreadLowerAsymptoteMin` to it ([powerLawSmoothLink()]). `spreadPowerLawPcBounds` are the bounds
#' `fireSense_spreadFit` gives `p_c`, the `maxAsymptote` slot, with the `"powerLaw"` link ([powerLawLink()]):
#' `p_c` must stay above the floor (0.01 above it), and its upper bound is 0.35. These are the outer limits:
#' `fireSense_spreadFit` bounds `p_c` above by the ELF's own calibrated threshold ([powerLawPcBounds()]).
#' Runaway-curve simulations on real ELF landscapes put the threshold at 0.258-0.288 (0.248 on a uniform
#' landscape), higher where less of the landscape can burn, so `spreadProbCeiling` (0.276) would bind.
#'
#' @format Numeric scalars.
#' @name spreadProbBounds
#' @aliases spreadProbCeiling spreadProbFloor spreadPowerLawPcBounds
#' @export
spreadProbCeiling <- 0.276

#' @rdname spreadProbBounds
#' @export
spreadProbFloor <- 0.13

#' @rdname spreadProbBounds
#' @export
spreadPowerLawPcBounds <- c(spreadProbFloor + 0.01, 0.35)

## An integer seed from a digest of `x`: the same inputs always give the same seed, so a random step
## that is seeded with it (withr::with_seed()) is reproducible without fixing a seed by hand.
## `x`'s digest, or NULL when caching is off: the digests only feed Cache(.cacheExtra =), which is unused then
.digestIfCaching <- function(x, useCache) if (isFALSE(useCache)) NULL else reproducible::.robustDigest(x)

.digestSeed <- function(x) strtoi(substr(reproducible::.robustDigest(x), 1, 7), 16L)

#' Log with a minimum
#'
#' Used for transforming Biomass to the log scale
#'
#' @param x Any value to be adjusted with log and a minimum B
#'
#' @return The original vector, logged with a minimum.
#'
#' @export
logMinB <- function(x) {
  minimumB <- exp(.logMinBFloor)
  x[x < minimumB] <- minimumB
  x <- log(x)
}

## the floor logMinB() applies, on the log scale (3.6), and how close to it a value is "on" it; the one
## place they are defined, for logMinB(), fuelLogToLinear() and emptySpreadCovariates()
.logMinBFloor <- log(100) - 1
.logMinBTol <- 1e-3

#' Range that puts linear fuel biomass on the scale the spread fit uses
#'
#' Fuel covariates of the spread model are biomass on the LINEAR scale, divided by a fixed
#' `1e4`. The division is not done to the data: it is the `covMinMax` given to
#' [rescaleKnown2()], which is affine and does not clamp, so `c(0, 1e4)` is exactly
#' `biomass / 1e4`. The constant is fixed, not `max(biomass)`, so that
#' `fireSense_spreadPredict` reproduces the fit's scaling in every year from the stored
#' `covMinMax_spread` alone. It is also how a fit is recognised as linear: see
#' [isLinearFuelRange()].
#'
#' @export
fuelLinearRange <- c(0, 1e4)

#' Fuel biomass from the log scale back to the linear scale
#'
#' `fireSenseCovariatesCreate()` returns fuel-class biomass as [logMinB()]: logged, with
#' everything below `exp(log(100) - 1)` (36.8) raised to that floor. On that scale the spread
#' model cannot use the fuel gradient: on ELF 5.3.2, 45% of pixels sit on the floor and the
#' treed ones fall in 16% of the covariate range (rescaled sd 0.05), so the coefficient
#' estimates little more than treed against treeless. On the linear scale the treed pixels
#' span twice the range and an absent fuel is exactly 0, so it contributes exactly nothing.
#'
#' This undoes the log where the spread model needs it, instead of removing it at the source,
#' because `fireSenseCovariatesCreate()` also builds the ignition covariates and its output is
#' cached for every fitted polygon. Both `fireSense_spreadFit` and `fireSense_spreadPredict`
#' call this function, so the fit and the prediction cannot drift apart.
#'
#' @param x Numeric vector of fuel biomass as returned by [logMinB()].
#'
#' @return Biomass on the linear scale; values on the [logMinB()] floor become 0.
#' @export
fuelLogToLinear <- function(x) {
  out <- exp(x)
  out[x <= .logMinBFloor + .logMinBTol] <- 0
  out
}

#' Was this covariate rescaled as linear fuel biomass?
#'
#' @param range Length-2 numeric: a covariate's entry in `covMinMax`.
#'
#' @return `TRUE` when `range` is [fuelLinearRange], i.e. the fit used linear fuel biomass, so
#'   a prediction must apply [fuelLogToLinear()] to that covariate before rescaling it. `FALSE`
#'   for a fit made on the log scale, whose covariates must be left as they are.
#' @export
isLinearFuelRange <- function(range) {
  length(range) == 2L && isTRUE(all.equal(as.numeric(range), fuelLinearRange))
}
