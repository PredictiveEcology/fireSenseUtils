#' Re-simulate the fires cut at an open buffer edge on the whole landscape
#'
#' A replicate fire that burns an open edge pixel of its buffer is cut: it would have kept burning, and the
#' objective scores it at the size it reached in the buffer. With whole-landscape covariate tables
#' (`annualDTx1000Whole`, `nonAnnualDTx1000Whole` of [.objfunSpreadFit()]) each such fire is run again on the
#' whole landscape, from its ignition and with its replicate's year effect, so that every term scores the size
#' the fire reaches unbounded.
#'
#' This is exact by coupling. Inside the buffer, the buffered fire and the same fire on the whole landscape follow
#' one process until the fire burns an open edge pixel; so "the buffered fire is cut" is "the whole-landscape fire
#' burns an open edge pixel". A cut fire is therefore a whole-landscape fire given that it burns one, which is
#' sampled by rejection: re-run until it does, at most `maxTries` times. A fire still short of an open edge pixel
#' after `maxTries` keeps its buffered size and stays cut.
#'
#' Departures from running every fire on the whole landscape: the cells burned by a replicate's uncut fires are
#' unburnable for its re-simulated fires, so a cut fire no longer competes with the fires simulated beside it; and
#' re-simulated fires of one replicate are accepted one try at a time, as they burn an open edge pixel.
#'
#' @param spreadState The replicates' burned pixels: `rep`, `initialLocus`, `indices` (landscape cells).
#' @param simSizes One row per replicate fire: `rep`, `initialLocus`, `N` (size) and `cut` (burned an open edge pixel).
#' @param openPixels Cells of the year's open buffer edge pixels (`openEdge$pixelID`).
#' @param loci Ignition cells of the year's fires.
#' @param eps The year effect drawn in each replicate (`0` without one), in replicate order.
#' @param wholeLp A function returning the whole landscape's `pixelID` and linear predictor `lp` for the year, without the year effect; it is called only when a fire is cut.
#' @param toP A function from a linear predictor (with the year effect added) to a spread probability.
#' @param r The landscape, only its geometry is used.
#' @param escapeMinPx,jumpTries,jumpMeanDist As for [spreadFire()].
#' @param maxTries The most re-simulations of a fire to find one that burns an open edge pixel.
#' @return A list: `simSizes`, with the re-simulated `N` and `cut = FALSE` for every accepted fire and
#'   `Nbuffered`, the buffered size; `resimulated`, the number of fires accepted; `capHit`, the number that
#'   used up `maxTries` and stay cut.
#' @keywords internal
resimulateCutFires <- function(spreadState, simSizes, openPixels, loci, eps, wholeLp, toP, r,
                               escapeMinPx, jumpTries, jumpMeanDist, maxTries) {
  simSizes <- data.table::copy(simSizes)
  set(simSizes, NULL, "Nbuffered", simSizes$N)
  if (!any(simSizes$cut)) return(list(simSizes = simSizes, resimulated = 0L, capHit = 0L))
  w <- wholeLp()
  crop <- cropToCells(r, c(w$pixelID, loci))
  pixCrop <- crop$toCrop(w$pixelID)
  lociCrop <- crop$toCrop(loci)
  resimulated <- 0L
  for (rp in unique(simSizes$rep[simSizes$cut])) {
    cutLoci <- simSizes$initialLocus[simSizes$rep == rp & simSizes$cut]
    sp <- numeric(crop$ncell)
    sp[pixCrop] <- toP(w$lp + eps[rp])
    sp[lociCrop] <- 1
    ## the replicate's uncut fires keep their cells
    uncut <- spreadState$indices[spreadState$rep == rp & !spreadState$initialLocus %in% cutLoci]
    sp[crop$toCrop(uncut)] <- 0
    todo <- cutLoci
    nTry <- 0L
    while (length(todo) && nTry < maxTries) {
      nTry <- nTry + 1L
      o <- spreadFire(crop$r, crop$toCrop(todo), sp, escapeMinPx, jumpTries, jumpMeanDist)
      o <- data.table::data.table(initialLocus = crop$toFull(o$initialLocus), pixelID = crop$toFull(o$indices))
      g <- o[, list(N = .N, hit = any(pixelID %in% openPixels)), by = "initialLocus"]
      ok <- g[which(hit)]
      if (!NROW(ok)) next
      row <- which(simSizes$rep == rp & simSizes$initialLocus %in% ok$initialLocus)
      set(simSizes, row, "N", ok$N[match(simSizes$initialLocus[row], ok$initialLocus)])
      set(simSizes, row, "cut", FALSE)
      sp[crop$toCrop(o$pixelID[o$initialLocus %in% ok$initialLocus])] <- 0 # an accepted fire keeps its cells
      todo <- setdiff(todo, ok$initialLocus)
      resimulated <- resimulated + length(row)
    }
  }
  list(simSizes = simSizes, resimulated = resimulated, capHit = sum(simSizes$cut))
}

#' One spread of fires from `loci`, as the objective runs them
#'
#' The `spreadCpp()` call of the spread objective ([.objfunSpreadFit()]), shared with
#' [resimulateCutFires()]. With `escapeMinPx` an escaped fire burns its first `escapeMinPx` cells whatever
#' `spreadProb`.
#'
#' @param landscape A (cropped) `SpatRaster`.
#' @param loci Ignition cells of `landscape`.
#' @param spreadProb Spread probability of each cell of `landscape`.
#' @param escapeMinPx,jumpTries,jumpMeanDist Minimum fire size (pixels), or `NULL` for none, and the jump
#'   settings that go with it.
#' @param stopArgs `NULL`, or a list with `stopCells`, `stopEvent` and `stopAt` for `spreadCpp()`.
#' @keywords internal
spreadFire <- function(landscape, loci, spreadProb, escapeMinPx = NULL, jumpTries = 0, jumpMeanDist = 0,
                       stopArgs = NULL) {
  if (is.null(escapeMinPx)) {
    SpaDES.tools::spreadCpp(
      landscape = landscape,
      loci = loci,
      spreadProb = spreadProb,
      stopCells = stopArgs$stopCells, stopEvent = stopArgs$stopEvent, stopAt = stopArgs$stopAt
    )
  } else { ## an escaped fire burns its first escapeMinPx cells whatever spreadProb
    SpaDES.tools::spreadCpp(
      landscape = landscape,
      loci = loci,
      spreadProb = spreadProb,
      minSize = escapeMinPx,
      jumpTries = jumpTries,
      jumpMeanDist = jumpMeanDist,
      stopCells = stopArgs$stopCells, stopEvent = stopArgs$stopEvent, stopAt = stopArgs$stopAt
    )
  }
}

#' The year's spread probability, from its covariates, as the objective computes it
#'
#' `spreadProbCovsOfYear()`, the design matrix and `spreadProbFromCovs()` of [.objfunSpreadFit()], shared with the
#' re-simulation of cut fires on the whole landscape.
#' @return A list: `dt`, the year's covariate table with a `spreadProb` column, and `mat`, its design matrix.
#' @keywords internal
yearSpreadProb <- function(yr, annDTx1000, nonAnnualDTx1000, indexNonAnnual, covMinMax, mutuallyExclusive,
                           colsToUse, doAssertions, logisticPars, covPars, maxFireSpread, lowerSpreadProb,
                           covCentre, link, climateRef, powerLawGamma, powerLawShape) {
  dt <- spreadProbCovsOfYear(
    yr, annDTx1000, nonAnnualDTx1000, indexNonAnnual, covMinMax, mutuallyExclusive, colsToUse,
    doAssertions, logisticPars, covPars, maxFireSpread, lowerSpreadProb, covCentre = covCentre
  )
  mat <- as.matrix(dt[, ..colsToUse])
  set(dt, NULL, "spreadProb", spreadProbFromCovs(mat, logisticPars, covPars, lowerSpreadProb, link,
                                                 climateRef = climateRef, powerLawGamma = powerLawGamma,
                                                 powerLawShape = powerLawShape))
  list(dt = dt, mat = mat)
}

## The year effect is added to the linear predictor with the powerLaw link, and to the logit of p otherwise:
## the linear predictor, and the function from it (with the year effect) to p
yearEffectLink <- function(link, mat, covPars, spreadProb, logisticPars, powerLawGamma, powerLawShape,
                           lowerSpreadProb) {
  powerLaw <- identical(link, "powerLaw")
  list(lp = if (powerLaw) as.vector(mat %*% covPars) else stats::qlogis(spreadProb),
       toP = if (powerLaw) function(z) powerLawLink(z, logisticPars[1L], powerLawGamma, lowerSpreadProb, powerLawShape) else stats::plogis)
}

## Both whole-landscape tables or neither, and not with the penalty that stops cut fires at the buffer on purpose
checkCutFireTables <- function(annualDTx1000Whole, nonAnnualDTx1000Whole, penaliseRunaways, annualDTx1000 = NULL,
                               nonAnnualDTx1000 = NULL) {
  if (is.null(annualDTx1000Whole) && is.null(nonAnnualDTx1000Whole)) return(invisible(FALSE))
  if (is.null(annualDTx1000Whole) || is.null(nonAnnualDTx1000Whole))
    stop("`annualDTx1000Whole` and `nonAnnualDTx1000Whole` go together: give both or neither", call. = FALSE)
  if (isTRUE(penaliseRunaways))
    stop("The whole-landscape tables (`annualDTx1000Whole`, `nonAnnualDTx1000Whole`) re-simulate cut fires, ",
         "while `penaliseRunaways = TRUE` stops them at the buffer on purpose: use `penaliseRunaways = FALSE`",
         call. = FALSE)
  if (!is.null(annualDTx1000) && !all(names(annualDTx1000) %in% names(annualDTx1000Whole)))
    stop("`annualDTx1000Whole` lacks the fire year(s) ",
         paste(setdiff(names(annualDTx1000), names(annualDTx1000Whole)), collapse = ", "), call. = FALSE)
  if (!is.null(nonAnnualDTx1000) && !identical(names(nonAnnualDTx1000Whole), names(nonAnnualDTx1000)))
    stop("`nonAnnualDTx1000Whole` must have the tables of `nonAnnualDTx1000`, by the same names", call. = FALSE)
  invisible(TRUE)
}
