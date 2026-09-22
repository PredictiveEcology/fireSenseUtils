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
#' @param link `NULL` (the default) to choose the form from `logisticPars` as above, or
#'   `"logistic3pUpper"` to use the 3-parameter form with Stukel's upper tail. Pass it explicitly
#'   wherever `logisticPars` may arrive without names, as inside the objective function.
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
logisticAll <- function(logisticPars, mat, covPars, lowerSpreadProb, link = NULL) {
  if (is.null(link) && "upperTail1" %in% names(logisticPars)) link <- "logistic3pUpper"
  if (!is.null(link)) {
    link <- match.arg(link, c("logistic3pUpper"))
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
  minimumB <- exp(log(100) - 1)
  x[x < minimumB] <- minimumB
  x <- log(x)
}

#' Range that puts linear fuel biomass on the scale the spread fit uses
#'
#' Fuel covariates of the spread model are biomass on the LINEAR scale, divided by a fixed
#' `1e4`. The division is not done to the data: it is the `covMinMax` given to
#' [rescaleKnown2()], which is affine and does not clamp, so `c(0, 1e4)` is exactly
#' `biomass / 1e4`. The constant is fixed, not `max(biomass)`, so that
#' `fireSense_SpreadPredict` reproduces the fit's scaling in every year from the stored
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
#' cached for every fitted polygon. Both `fireSense_SpreadFit` and `fireSense_SpreadPredict`
#' call this function, so the fit and the prediction cannot drift apart.
#'
#' @param x Numeric vector of fuel biomass as returned by [logMinB()].
#'
#' @return Biomass on the linear scale; values on the [logMinB()] floor become 0.
#' @export
fuelLogToLinear <- function(x) {
  floorLog <- log(100) - 1 # the floor logMinB() applies
  out <- exp(x)
  out[x <= floorLog + 1e-3] <- 0
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
