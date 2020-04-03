utils::globalVariables(c("..colsToUse", ".N", "buffer", "N", "pixelID", "spreadProb"))

#' Objective function for \code{fireSense_spreadFit} module
#'
#' @param par parameters
#' @param landscape DESCRIPTION NEEDED
#' @param annualDTx1000 DESCRIPTION NEEDED
#' @param nonAnnualDTx1000 DESCRIPTION NEEDED
#' @param pixelIndices DESCRIPTION NEEDED
#' @param formula DESCRIPTION NEEDED
#' @param historicalFires DESCRIPTION NEEDED
#' @param fireBufferedListDT DESCRIPTION NEEDED
#' @param wADtest DESCRIPTION NEEDED
#' @param verbose DESCRIPTION NEEDED
#' @param maxFireSpread DESCRIPTION NEEDED
#'
#' @return DESCRIPTION NEEDED
#'
#' @export
#' @importFrom data.table := rbindlist set setDT setDTthreads setnames
#' @importFrom kSamples ad.test
#' @importFrom purrr transpose
#' @importFrom raster ncell
#' @importFrom SpaDES.tools spread
#' @importFrom stats dbinom median terms
#' @importFrom utils tail
.objfun <- function(par,
                    landscape,
                    annualDTx1000,
                    nonAnnualDTx1000,
                    pixelIndices,
                    formula,
                    historicalFires,
                    fireBufferedListDT,
                    wADtest = 1, # not used yet
                    verbose = TRUE,
                    maxFireSpread = 0.255) {
  # Optimization's objective function
  data.table::setDTthreads(1)
  if (missing(landscape)) {
    landscape <- get("landscape", envir = .GlobalEnv)
  }
  if (missing(annualDTx1000)) {
    annualDTx1000 <- get("annualDTx1000", envir = .GlobalEnv)
  }
  if (missing(nonAnnualDTx1000)) {
    nonAnnualDTx1000 <- get("nonAnnualDTx1000", envir = .GlobalEnv)
  }
  if (missing(historicalFires)) {
    historicalFires <- get("historicalFires", envir = .GlobalEnv)
  }
  if (missing(fireBufferedListDT)) {
    fireBufferedListDT <- get("fireBufferedListDT", envir = .GlobalEnv)
  }
  lapply(nonAnnualDTx1000, setDT)
  colsToUse <- attributes(terms(formula))[["term.labels"]]

  # How many of the parameters belong to the model?
  parsModel <- length(colsToUse)
  ncells <- ncell(landscape)
  r <- raster(landscape)
  years <- as.character(names(annualDTx1000))
  names(years) <- years
  cells <- integer(ncells)
  Nreps <- 10
  yearSplit <- strsplit(names(nonAnnualDTx1000), "_")
  names(yearSplit) <- as.character(seq_along(nonAnnualDTx1000))
  indexNonAnnual <- rbindlist(
    Map(
      ind = seq_along(nonAnnualDTx1000), date = yearSplit,
      function(ind, date) data.table(ind = ind, date = date)
    )
  )
  results <- Map(
    annDTx1000 = annualDTx1000,
    yr = years,
    annualFires = historicalFires,
    annualFireBufferedDT = fireBufferedListDT,
    MoreArgs = list(
      par = par, parsModel = parsModel,
      # pixelIndices = pixelIndices,
      verbose = verbose,
      nonAnnualDTx1000 = nonAnnualDTx1000,
      indexNonAnnual = indexNonAnnual,
      colsToUse = colsToUse
    ),
    function(yr, annDTx1000, par, parsModel,
             annualFires, nonAnnualDTx1000, annualFireBufferedDT, # pixelIndices,
             indexNonAnnual, colsToUse,
             verbose = TRUE) {
      # Rescale to numerics and /1000
      setDT(annDTx1000)
      # needed because data.table objects were recovered from disk
      shortAnnDTx1000 <- nonAnnualDTx1000[[indexNonAnnual[date == yr]$ind]][annDTx1000, on = "pixelID"]
      mat <- as.matrix(shortAnnDTx1000[, ..colsToUse]) / 1000
      # matrix multiplication
      covPars <- tail(x = par, n = parsModel)
      logisticPars <- par[1:4]
      set(shortAnnDTx1000, NULL, "spreadProb", logistic4p(mat %*% covPars, logisticPars))

      # logistic multiplication
      set(annDTx1000, NULL, "spreadProb", logistic4p(annDTx1000$pred, par[1:4]))

      medSP <- median(shortAnnDTx1000[, mean(spreadProb, na.rm = TRUE)], na.rm = TRUE)
      if (medSP <= maxFireSpread & medSP >= 0.16) {
        if (verbose) {
          print(paste0(
            Sys.getpid(), "-- year: ", yr, ", spreadProb raster: median in buffered pixels = ",
            round(medSP, 3)
          ))
        }
        cells[as.integer(shortAnnDTx1000$pixelID)] <- shortAnnDTx1000$spreadProb
        maxSizes <- rep(annualFires$size, times = Nreps) * 2
        lociAll <- rep(annualFires$cells, times = Nreps)

        spreadState <- SpaDES.tools::spread(
          landscape = r,
          maxSize = maxSizes,
          loci = rep(annualFires$cells, times = Nreps),
          spreadProb = cells,
          returnIndices = TRUE,
          allowOverlap = TRUE,
          quick = TRUE
        )
        fireSizes <- tabulate(spreadState[["id"]]) # Here tabulate() is equivalent to table() but faster
        if (length(fireSizes) == 0) {
          message("We have a fire size == 0. Entering debug mode.")
          browser()
        }
        burnedProb <- spreadState[, .N, by = "indices"]
        setnames(burnedProb, "indices", "pixelID")
        setDT(annualFireBufferedDT)
        out <- burnedProb[annualFireBufferedDT, on = "pixelID"]
        out[is.na(N), N := 0]

        # THis is a work around for cases where "initial pixels" are not actually burned in
        #   the polygon database
        out[pixelID %in% annualFires$cells, buffer := 1]
        # Add a very small number so that no pixel has exactly zero probability -- creating Inf
        SNLL <- -sum(dbinom(
          prob = pmin(out$N / Nreps + 0.001, 0.99), size = 1,
          x = out$buffer,
          log = TRUE
        ), na.rm = TRUE) # Sum of the negative log likelihood
      } else {
        SNLL <- 1e7
        fireSizes <- sample(1:3, 1)
      }

      list(fireSizes = fireSizes, SNLL = SNLL)
    }
  )
  results <- purrr::transpose(results)
  historicalFiresTr <- purrr::transpose(historicalFires)

  adTest <- try(ad.test(unlist(results$fireSizes), unlist(historicalFiresTr$size))[["ad"]][1L, 1L])
  SNLLTest <- sum(unlist(results$SNLL))
  objFunRes <- adTest + SNLLTest / 1e3 # wADtest is the weight for the AD test
  # gc()
  # Figure out what we want from these. This is potentially correct (i.e. we want the smallest ad.test and the smallest SNLL)
  return(objFunRes)
}
