utils::globalVariables(c("..colsToUse", ".N", "buffer", "burned", "burnedClass", "N",
                         "pixelID", "prob", "spreadProb"))

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
#' @param covMinMax DESCRIPTION NEEDED
#' @param maxFireSpread DESCRIPTION NEEDED
#' @param tests One or more of "mad", "adTest", or "SNLL" as a character vector. Default is "mad"
#' @param verbose DESCRIPTION NEEDED
#'
#' @return DESCRIPTION NEEDED
#'
#' @export
#' @importFrom data.table := rbindlist set setDT setDTthreads setnames
#' @importFrom kSamples ad.test
#' @importFrom purrr transpose pmap
#' @importFrom raster ncell
#' @importFrom SpaDES.tools spread
#' @importFrom stats dbinom median terms
#' @importFrom utils tail
.objfun <- function(par,
                    landscape,
                    annualDTx1000,
                    nonAnnualDTx1000,
                    pixelIndices,
                    formula, #loci, sizes,
                    historicalFires,
                    fireBufferedListDT,
                    covMinMax = NULL,
                    maxFireSpread = 0.25, # 0.257 makes gigantic fires
                    minFireSize = 2,
                    tests = "mad",
                    #bufferedRealHistoricalFiresList,
                    verbose = TRUE){ #fireSense_SpreadFitRaster
  # Optimization's objective function
  # lapply(historicalFires, setDT)
  data.table::setDTthreads(1)
  doMADTest <- grepl("mad", tolower(tests))
  doSNLLTest <- grepl("snll", tolower(tests))
  doADTest <- grepl("adtest", tolower(tests))
  if (missing(landscape))
    landscape <- get("landscape", envir = .GlobalEnv)
  if (missing(annualDTx1000))
    annualDTx1000 <- get("annualDTx1000", envir = .GlobalEnv)
  if (missing(nonAnnualDTx1000))
    nonAnnualDTx1000 <- get("nonAnnualDTx1000", envir = .GlobalEnv)
  if (missing(historicalFires))
    historicalFires <- get("historicalFires", envir = .GlobalEnv)
  if (missing(fireBufferedListDT))
    fireBufferedListDT <- get("fireBufferedListDT", envir = .GlobalEnv)
  #lapply(annualDTx1000, setDT)
  lapply(nonAnnualDTx1000, setDT)
  #lapply(fireBufferedListDT, setDT)
  # dtThreadsOrig <- data.table::setDTthreads(1)
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
    Map(ind = seq_along(nonAnnualDTx1000), date = yearSplit,
        function(ind, date) data.table(ind = ind, date = date))
  )
  historicalFiresAboveMin <- lapply(historicalFires, function(x) {
    x <- x[x$size >= minFireSize,]
    x <- x[!duplicated(x$cells),]
    x
    }
  )
  results <- purrr::pmap(
    list(annDTx1000 = annualDTx1000,
         yr = years,
         annualFires = historicalFiresAboveMin,
         annualFireBufferedDT = fireBufferedListDT),
    par = par, parsModel = parsModel,
    #pixelIndices = pixelIndices,
    verbose = verbose,
    nonAnnualDTx1000 = nonAnnualDTx1000,
    indexNonAnnual = indexNonAnnual,
    colsToUse = colsToUse,
    covMinMax = covMinMax,
    .f = function(yr, annDTx1000, par, parsModel,
                  annualFires, nonAnnualDTx1000, annualFireBufferedDT, #pixelIndices,
                  indexNonAnnual, colsToUse, covMinMax,
                  verbose = TRUE) {
      # needed because data.table objects were recovered from disk

      # Rescale to numerics and /1000
      #setDT(nonAnnDTx1000)
      setDT(annDTx1000)
      shortAnnDTx1000 <- nonAnnualDTx1000[[indexNonAnnual[date == yr]$ind]][annDTx1000, on = "pixelID"]
      # browser(expr = yr == 1991)
      if (!is.null(covMinMax)) {
        for (cn in colnames(covMinMax)) {
          set(shortAnnDTx1000, NULL, cn,
              rescaleKnown(shortAnnDTx1000[[cn]], 0, 1000, covMinMax[[cn]][1], covMinMax[[cn]][2]))
        }
      }
      mat <- as.matrix(shortAnnDTx1000[, ..colsToUse])/1000
      # matrix multiplication
      covPars <- tail(x = par, n = parsModel)
      logisticPars <- par[1:4]
      set(shortAnnDTx1000, NULL, "spreadProb", logistic4p(mat %*% covPars, logisticPars))
      # logistic multiplication
      # set(annDTx1000, NULL, "spreadProb", logistic4p(annDTx1000$pred, par[1:4])) ## 5-parameters logistic
      #set(annDTx1000, NULL, "spreadProb", logistic5p(annDTx1000$pred, par[1:5])) ## 5-parameters logistic
      #actualBurnSP <- annDTx1000[annualFireBufferedDT, on = "pixelID"]
      medSP <- median(shortAnnDTx1000[, mean(spreadProb, na.rm = TRUE)], na.rm = TRUE)
      # browser(expr = yr == 2014)
      if (medSP <= maxFireSpread & medSP >= 0.16) {
        if (verbose) {
          print(paste0(Sys.getpid(), "-- year: ",yr, ", spreadProb raster: median in buffered pixels = ",
                       round(medSP, 3)))
        }
        cells[as.integer(shortAnnDTx1000$pixelID)] <- shortAnnDTx1000$spreadProb
        #maxSizes <- rep(annualFires$size, times = Nreps)

        # this will make maxSizes be a little bit larger for large fires, but a lot bigger for small fires
        #maxSizes <- maxSizes * 1.5#(1.1+pmax(0,5-log10(maxSizes)))

        maxSizes <- multiplier(annualFires$size)
        #maxSizes <- annualFires$size * 2
        loci <- annualFires$cells
        if (any(cells[loci] == 0)) {
          cells[loci] <- 1
        }
        dups <- duplicated(annualFires$cells)
        if (any(dups)) {
          maxSizes <- maxSizes[!dups]
          loci <- annualFires$cells[!dups]
        }

        spreadState <- lapply(seq_len(Nreps), function(i) {
            SpaDES.tools::spread(
              landscape = r,
              maxSize = maxSizes,
              loci = loci,
              spreadProb = cells,
              returnIndices = TRUE,
              allowOverlap = FALSE,
              quick = TRUE)
          })
        spreadState <- rbindlist(spreadState, idcol = "rep")
        fireSizes <- round(tabulate(spreadState[["id"]])/Nreps,0) # Here tabulate() is equivalent to table() but faster
        # if (length(fireSizes) == 0) browser()
        burnedProb <- spreadState[, .N, by = "indices"]
        setnames(burnedProb, "indices", "pixelID")
        setDT(annualFireBufferedDT)
        out <- burnedProb[annualFireBufferedDT, on = "pixelID"]

        # fix the out
        # 1 -- set pixels that had not simulated fires to N = 0
        out[is.na(N), N := 0]
        # 2 -- rescale probability surface between 0.001 and 0.99
        #      so probabilities can be calculated
        out[, prob := pmin(out$N/Nreps + 0.001, 0.99)]
        # 3 -- convert buffer (which has 1 in buffer) to burned = 1 - buffer
        out[, burned := buffer]
        # 4 -- Set initial pixels to burned = 2 -- is a work around for cases where "initial pixels" are not actually burned in
        #   the polygon database
        out[, burnedClass := burned]
        out[pixelID %in% annualFires$cells, burnedClass := 2]

        if (FALSE) { # THIS IS PLOTTING STUFF
          r <- raster(landscape)
          r[out$pixelID] <- out$prob
          clearPlot();Plot(r)
          # ex <- clickExtent()
          ex <- new("Extent", xmin = -1130927.72835113, xmax = -1029209.34163701,
                    ymin = 8098144.00948992, ymax = 8224186.35824437)
          # ex <- new("Extent", xmin = -1098283.46889952, xmax = -1037633.32535885,
          #            ymin = 7969991.96172249, ymax = 8030642.10526316)
          predictedFireProb <- crop(r, ex)
          # clearPlot();Plot(r1)
          actualFire <- raster(r)
          actualFire[out$pixelID] <- out$burnedClass
          actualFire[]
          actualFire <- crop(actualFire, ex)
          levels(actualFire) <- data.frame(ID = 0:2, class = c("unburned", "burned", "ignited"))
          predictedLiklihood <- dbinom(prob = out$prob,
                                       size = 1,
                                       x = out$burned,
                                       log = TRUE
          )
          spreadProbMap <- raster(r)
          spreadProbMap[out$pixelID] <- cells[out$pixelID]
          spreadProbMap <- crop(spreadProbMap, ex)
          predLiklihood <- raster(r)
          predLiklihood[out$pixelID] <- predictedLiklihood
          predLiklihood <- crop(predLiklihood, ex)
          clearPlot(); Plot(actualFire, predictedFireProb, predLiklihood, spreadProbMap)
          Plot(predLiklihood, cols = "RdYlGn", new = TRUE, legendRange = range(round(predLiklihood[], 0), na.rm = TRUE))
        }
        # Add a very small number so that no pixel has exactly zero probability -- creating Inf
        # SNLL <- -sum(dbinom(prob = out$prob,
        #                     size = 1,
        #                     x = out$burned,
        #                     log = TRUE
        # ), na.rm = TRUE) # Sum of the negative log likelihood
      } else {
        #SNLL <- 1e7
        fireSizes <- sample(1:3, 1)
      }

      list(fireSizes = fireSizes)#, SNLL = SNLL)
    })
  results <- purrr::transpose(results)

  mess <- character()
  objFunRes <- 0
  if (isTRUE(doMADTest)) {
    a <- purrr::map2(historicalFiresAboveMin, results$fireSizes, function(x, y) {
      data.table(x, simFireSize = y)
    })

    a <- rbindlist(a)
    a[, dev := abs(size - simFireSize)]
    mad <- round(mean(a$dev), 1)
    objFunRes <- objFunRes + mad #+ SNLLTest
    mess <- paste(" mad:", mad, "; ")
  }
  if (isTRUE(doADTest)) {
    historicalFiresTr <- unlist(purrr::transpose(historicalFiresAboveMin)$size)
    simulatedFires <- unlist(results$fireSizes)
    adTest <- try(ad.test(simulatedFires, historicalFiresTr)[["ad"]][1L, 1L])
    objFunRes <- objFunRes + adTest #+ SNLLTest
    mess <- paste(mess, " adTest:", adTest, "; ")
  }
  if (isTRUE(doSNLLTest)) {
    rescalor <- 6e4
    SNLLTest <- sum(unlist(results$SNLL)) / rescalor
    mess <- paste(mess, " SNLLTest:", SNLLTest, "; ")
    objFunRes <- objFunRes + SNLLTest #+ SNLLTest

  }
  print(paste0("  ", Sys.getpid(), mess))
  # gc()
  # Figure out what we want from these. This is potentially correct (i.e. we want the smallest ad.test and the smallest SNLL)
  return(objFunRes)
}

rescaleKnown <- function(x, minNew, maxNew, minOrig, maxOrig) {
  a1 <- x - minOrig # brings min to zero
  a2 <- a1 * maxNew/max(a1)
  a2
}
