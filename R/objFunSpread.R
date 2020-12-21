utils::globalVariables(c(
  "..colsToUse", ".N", "buffer", "burned", "burnedClass", "id", "ids", "N",
  "pixelID", "prob", "simFireSize", "size", "spreadProb"))

#' Objective function for \code{fireSense_spreadFit} module
#'
#' @param par parameters
#' @param landscape A RasterLayer with extent, res, proj used for SpaDES.tools::spread
#' @param annualDTx1000 A list of data.table class objects. Each list element is
#'   data from a single calendar year, and whose name is "yearxxxx" where xxxx is the 4 number
#'   year. The columns in the data.table must integers, that are 1000x their actual values as
#'   this function will divide by 1000.
#' @param nonAnnualDTx1000 Like \code{annualDTx1000}, but with where each list element will be
#'   used for >1 year. The names of the list elements must be "yearxxxx_yearyyyy_yearzzzz" where the
#'   xxxx, yyyy, or zzzz represent the calendar years for which that list element should be used.
#'   The columns are variables that are used for more than 1 year.
#' @param FS_formula Formula, put provided as a character string, not class \code{formula}.
#'   (if it is provided as a class \code{formula}, then it invariably will have an
#'   enormous amount of data hidden in the formula environment; this is bad for DEoptim)
#' @param historicalFires DESCRIPTION NEEDED
#' @param fireBufferedListDT DESCRIPTION NEEDED
#' @param covMinMax This is a 2 row by multiple column data.frame indicating
#'   the minimum and maximum values of the original covariate data values. These will
#'   be used to rescale the covariates internally so that they are all between 0 and 1. It is important
#'   to not simply rescale internally here because only 1 year is run at a time; all years
#'   must be rescaled for a given covariate by the same amount.
#' @param maxFireSpread A value for spreadProb that is considered impossible to go above.
#'   Default 0.28, which is overly generous unless there are many non-flammable pixels (e.g., lakes).
#' @param minFireSize DESCRIPTION NEEDED
#' @param tests One or more of \code{"mad"}, \code{"adTest"}, \code{"SNLL"}, or \code{"SNLL_FS"}.
#'              Default: \code{"mad"}.
#' @param Nreps Integer. The number of replicates, per ignition, to run.
#' @param plot.it DESCRIPTION NEEDED
#' @param objFunCoresInternal Internally, this function can use \code{mcmapply} to run multiple
#'   parallel \code{spread} function calls. This should only be >1L if there are spare threads.
#'   It is highly likely that there won't be. However, sometimes the \code{DEoptim} is
#'   particularly inefficient, it starts X cores, and immediately several of them are
#'   stopped inside this function because the parameters are so bad, only 2 year are attempted.
#'   Then the core will stay idle until all other cores for the \code{DEoptim} iteration are complete.
#'   Similarly, if only physical cores are used for \code{DEoptim}, the additional use of
#'   hyperthreaded cores here, internally will speed things up (i.e., this maybe could be 2L or 3L).
#' @param verbose DESCRIPTION NEEDED
#'
#' @return DESCRIPTION NEEDED
#'
#' @export
#' @importFrom data.table := rbindlist set setDT setDTthreads setnames setorderv
#' @importFrom EnvStats demp
#' @importFrom graphics abline axis hist mtext
#' @importFrom kSamples ad.test
#' @importFrom purrr map2 pmap transpose
#' @importFrom quickPlot clearPlot dev gpar Plot
#' @importFrom raster buffer crop extent ncell raster trim xyFromCell
#' @importFrom sp SpatialPoints
#' @importFrom SpaDES.tools spread
#' @importFrom stats as.formula dbinom median terms quantile
#' @importFrom utils tail
.objfunSpreadFit <- function(par,
                             landscape,
                             annualDTx1000,
                             nonAnnualDTx1000,
                             FS_formula, #loci, sizes,
                             historicalFires,
                             fireBufferedListDT,
                             covMinMax = NULL,
                             maxFireSpread = 0.28, # 0.257 makes gigantic fires
                             minFireSize = 2,
                             tests = "mad",
                             Nreps = 10,
                             plot.it = FALSE,
                             objFunCoresInternal = 1,
                             #bufferedRealHistoricalFiresList,
                             verbose = TRUE) { #fireSense_SpreadFitRaster
  # Optimization's objective function
  # lapply(historicalFires, setDT)

  data.table::setDTthreads(1)

  doMADTest <- any(grepl("mad", tolower(tests)))
  doSNLLTest <- any(grepl("snll$", tolower(tests)))
  doSNLL_FSTest <- any(grepl("snll_fs", tolower(tests)))
  doADTest <- any(grepl("adtest", tolower(tests)))
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
  if (is(FS_formula, "formula"))
    stop("FS_formula must be provided as a charater string because it takes too much RAM otherwise.")
  FS_formula <- as.formula(FS_formula)
  colsToUse <- attributes(terms(FS_formula))[["term.labels"]]
  # How many of the parameters belong to the model?
  parsModel <- length(colsToUse)
  ncells <- ncell(landscape)

  r <- raster(landscape)
  years <- as.character(names(annualDTx1000))
  names(years) <- years
  cells <- integer(ncells)
  # Nreps <- 10
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
  lowerSpreadProb <- 0.13
  fireSizesByYear <- unlist(lapply(historicalFiresAboveMin, function(x) sum(x$size)))
  largest <- head(sort(fireSizesByYear, decreasing = TRUE), 2)#max(2, objFunCoresInternal))
  smallest <- setdiff(names(fireSizesByYear), names(largest))
  lrgSmallFireYears <- list(large = names(largest),
                            small = smallest)
  objFunResList <- list() # will hold objective function values --> which is now >1 for large, then small fires
  for (ii in seq(lrgSmallFireYears)) {
    yrs <- lrgSmallFireYears[[ii]]
    if (length(yrs)) {
      results <- parallel::mcmapply(
        mc.cores = min(length(years[yrs]), objFunCoresInternal),
        mc.preschedule = FALSE,
        SIMPLIFY = FALSE,
        # results <- purrr::pmap(
        annDTx1000 = annualDTx1000[yrs],
        yr = years[yrs],
        annualFires = historicalFiresAboveMin[yrs],
        annualFireBufferedDT = fireBufferedListDT[yrs],
        MoreArgs = list(par = par, parsModel = parsModel,
                        verbose = verbose,
                        nonAnnualDTx1000 = nonAnnualDTx1000,
                        indexNonAnnual = indexNonAnnual,
                        colsToUse = colsToUse,
                        covMinMax = covMinMax),
        # .f = function(yr, annDTx1000, par, parsModel,
        function(yr, annDTx1000, par, parsModel,
                 annualFires, nonAnnualDTx1000, annualFireBufferedDT,
                 indexNonAnnual, colsToUse, covMinMax,
                 verbose = TRUE) {

          # needed because data.table objects were recovered from disk
          # Rescale to numerics and /1000
          #setDT(nonAnnDTx1000)
          setDT(annDTx1000)
          shortAnnDTx1000 <- nonAnnualDTx1000[[indexNonAnnual[date == yr]$ind]][annDTx1000, on = "pixelID"]
          if (!is.null(covMinMax)) {
            for (cn in colnames(covMinMax)) {
              set(shortAnnDTx1000, NULL, cn,
                  rescaleKnown(shortAnnDTx1000[[cn]], 0, 1000, covMinMax[[cn]][1], covMinMax[[cn]][2]))
            }
          }
          mat <- as.matrix(shortAnnDTx1000[, ..colsToUse])/1000
          # matrix multiplication
          covPars <- tail(x = par, n = parsModel)
          logisticPars <- head(x = par, n = length(par) - parsModel)
          if (logisticPars[1] > maxFireSpread) warning("The first parameter of the logistic is > ", maxFireSpread, ". The parameter should be lowered.")
          if (length(logisticPars) == 4) {
            set(shortAnnDTx1000, NULL, "spreadProb", logistic4p(mat %*% covPars, logisticPars))
          } else if (length(logisticPars) == 3) {
            set(shortAnnDTx1000, NULL, "spreadProb", logistic3p(mat %*% covPars, logisticPars, par1 = lowerSpreadProb))
          } else if (length(logisticPars) == 2) {
            set(shortAnnDTx1000, NULL, "spreadProb", logistic2p(mat %*% covPars, logisticPars, par1 = lowerSpreadProb, par4 = 0.5))
          }
          # logistic multiplication
          # set(annDTx1000, NULL, "spreadProb", logistic4p(annDTx1000$pred, par[1:4])) ## 5-parameters logistic
          #set(annDTx1000, NULL, "spreadProb", logistic5p(annDTx1000$pred, par[1:5])) ## 5-parameters logistic
          #actualBurnSP <- annDTx1000[annualFireBufferedDT, on = "pixelID"]
          medSP <- median(shortAnnDTx1000[, mean(spreadProb, na.rm = TRUE)], na.rm = TRUE)
          out <- tryCatch(if (medSP <= maxFireSpread & medSP >= lowerSpreadProb) {}, error = function(x) "error")
          if (identical(out, "error")) {
            # browser()
          }

          if (medSP <= maxFireSpread & medSP >= lowerSpreadProb) {

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
              annualFires <- annualFires[which(!dups),]#
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
            ret <- list()
            if (isTRUE(doSNLL_FSTest)) {
              emp <- spreadState[, list(N = .N, id = id[1]), by = c("rep", "initialLocus")]
              emp <- emp[annualFires, on = c("initialLocus" = "cells")]
              if (plot.it) {
                maxX <- log(max(emp$N))
                par(mfrow = c(7,7), omi = c(0.5, 0, 0, 0),
                    mai = c(0.2,0.3,0.4,0.1));
                emp <- setorderv(emp, c("size"), order = -1L)
                emp[, {
                  dat <- round(log(N))
                  h <- hist(dat, breaks = 30,
                            main = paste(as.character(.BY)),
                            # main = "",
                            axes = FALSE, xlim = c(0, maxX))
                  seqq <- seq(0, ceiling(maxX), by = 1)
                  axis(1, at = seqq, labels = round(exp(seqq),0))
                  abline(v = log(size[1]), col = "red")
                }, by = "initialLocus"]
                mtext(outer = TRUE,
                      paste("Fire year:",yr,"; Simulated fire sizes (# pixels); Actual Fire (red); Sorted by actual fire size"),
                      line = 2, side = 1)
              }

              # only use fires that escaped --> i.e., greater than 1 pixel
              emp <- emp[N > 1, list(lik = EnvStats::demp(x = size[1], obs = N)), by = "initialLocus"]
              minLik <- 1e-7#min(emp$lik[emp$lik > 0])
              set(emp, NULL, "lik", log(pmax(minLik, emp$lik)))
              SNLL_FS <- -sum(emp$lik)
              ret <- append(ret, list(SNLL_FS = SNLL_FS))
            }

            if (isTRUE(doMADTest)) {
              fireSizes <- round(tabulate(spreadState[["id"]])/Nreps,0) # Here tabulate() is equivalent to table() but faster
              ret <- append(ret, list(fireSizes = fireSizes))
            }
            if (isTRUE(doSNLLTest)) {
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
              if (isTRUE(plot.it)) { # THIS IS PLOTTING STUFF
                bigFire1 <- raster(r)
                bigFire1[out$pixelID] <- out$ids
                keepFire <- tail(sort(table(out$ids)),4)
                theseFires <- annualFires[ids %in% names(keepFire)]
                clearPlot();
                firesToDo <- theseFires$ids
                names(firesToDo) <- firesToDo
                out2 <- lapply(firesToDo, function(id) {
                  keepFire <- as.numeric(id)
                  # keepFire <- 65
                  bigFire <- bigFire1
                  bigFire[bigFire != keepFire] <- NA
                  bf <- trim(bigFire)
                  ex <- extent(bf)

                  thisFire <- annualFires[ids == keepFire]

                  r <- raster(landscape)
                  r[out$pixelID] <- out$prob

                  predictedFireProb <- crop(r, ex)
                  # clearPlot();Plot(r)
                  actualFire <- raster(r)
                  actualFire[out$pixelID] <- out$burnedClass
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
                  spreadProbMap[spreadProbMap >= par[1]] <- par[1]
                  ccc <- cells[out$pixelID];
                  ccc <- ccc[ccc > 0];
                  lowerLim <- quantile(ccc, 0.05);
                  ccc <- ccc[ccc > lowerLim];
                  spreadProbMap[spreadProbMap <= lowerLim] <- lowerLim
                  predLiklihood <- raster(r)
                  predLiklihood[out$pixelID] <- predictedLiklihood
                  predLiklihood <- crop(predLiklihood, ex)
                  spIgnits <- SpatialPoints(coords = xyFromCell(r, thisFire$cells))
                  spIgnits <- buffer(spIgnits, width = 5000)
                  spIgnits <- crop(spIgnits, ex)
                  list(spIgnits = spIgnits, predictedFireProb = predictedFireProb,
                       predLiklihood = predLiklihood,
                       spreadProbMap = spreadProbMap)
                })
                out3 <- purrr::transpose(out2)
                notSp <- grep("spIgnits", names(out3), value = TRUE, invert = TRUE)

                out4 <- unlist(out3[notSp], recursive = F)
                a <- Plot(out4)
                nn <- lapply(names(out3$spIgnits), function(id) {
                  spDat <- out3$spIgnits[[id]]
                  Plot(spDat, addTo = grep(id, names(out4), value = TRUE)[2],
                       gp = gpar(fill = rep("transparent", 10), col = "black"), title = "")
                })
                # grid::grid.newpage()

                # Plot(predictedFireProb, predLiklihood, spreadProbMap, title = "")
                # Plot(predictedFireProb, title = paste0("fire prob, date: ",yr, ", id: ", thisFire$cells), new = TRUE)
                # # Plot(predLiklihood, title = paste0("likelihood, date: ",yr, ", id: ", thisFire$cells), new = TRUE)
                # Plot(spreadProbMap, title = paste0("spreadProb, date: ",yr, ", id: ", thisFire$cells), new = TRUE)
                # # clearPlot(); Plot(actualFire, predictedFireProb, predLiklihood, spreadProbMap)
                # Plot(spIgnits, addTo = "spreadProbMap", gp = gpar(fill = rep("transparent", 10), col = "black"), title = "")
                # # Plot(spIgnits, addTo = "actualFire", gp = gpar(fill = rep("transparent", 10), col = "black"), title = "")
                # Plot(spIgnits, addTo = "predictedFireProb", gp = gpar(fill = rep("transparent", 10), col = "black"), title = "")
                # Plot(spIgnits, addTo = "predLiklihood", gp = gpar(fill = rep("transparent", 10), col = "black"), title = "")
                # Plot(predLiklihood, cols = "RdYlGn", new = TRUE,
                #      title = paste0("likelihood, date: ",yr, ", id: ", thisFire$cells),
                #      legendRange = range(round(predLiklihood[], 0), na.rm = TRUE))

              }
              # Add a very small number so that no pixel has exactly zero probability -- creating Inf
              # SNLL <- -sum(dbinom(prob = out$prob,
              #                     size = 1,
              #                     x = out$burned,
              #                     log = TRUE
              # ), na.rm = TRUE) # Sum of the negative log likelihood
            }
          } else {
            stop("encountered error with spreadProb - contact module developers")
            #Ian added this stop. Unclear what is supposed to happen. Object ret doesn't exist
            SNLL <- 1e7
            fireSizes <- sample(1:3, 1)
          }

          return(ret)
        }
      )
      results <- purrr::transpose(results)

      mess <- character()
      objFunRes <- 0

      if (isTRUE(doMADTest)) {
        a <- purrr::map2(historicalFiresAboveMin[yrs], results$fireSizes, function(x, y) {
          data.table(x, simFireSize = y)
        })

        a <- rbindlist(a)
        a[, dev := abs(size - simFireSize)]
        #a[, devLog := sqrt(abs(size - simFireSize))]
        mad <- round(mean(a$dev), 1)
        #mad <- round(mean(a$devLog), 1)
        objFunRes <- objFunRes + mad #+ SNLLTest
        mess <- paste(" mad:", mad, "; ")
        objFunResList[ii] <- list(list(objFunRes = objFunRes, nFires = NROW(a)))
        # if (mad > 2700) {
        #   print(paste0("  ", Sys.getpid(), mess))
        #   break
        # }

      }
      if (isTRUE(doADTest)) {
        historicalFiresTr <- unlist(purrr::transpose(historicalFiresAboveMin)$size)
        simulatedFires <- unlist(results$fireSizes)
        adTest <- try(ad.test(simulatedFires, historicalFiresTr)[["ad"]][1L, 1L])
        objFunRes <- objFunRes + adTest #+ SNLLTest
        mess <- paste(mess, " adTest:", adTest, "; ")
      }
      if (isTRUE(doSNLL_FSTest)) {
        SNLL_FSTest <- round(sum(unlist(results$SNLL)),1)
        failVal <- 1e5L
        threshold <- 550 * length(results$SNLL_FS) #parameterize?
        if (SNLL_FSTest > threshold && ii == 1) { #fine tune this threshold
          SNLL_FSTestOrig <- SNLL_FSTest
          SNLL_FSTest <- failVal
          mess <- paste('SNLL threshold:', threshold,'; ', "SNLL_FSTestOrig:", SNLL_FSTestOrig, "; ")
        }
        mess <- paste(mess, " SNLL_FSTest:", SNLL_FSTest, "; ")
        objFunRes <- SNLL_FSTest #+ SNLL_FSTest
        objFunResList[ii] <- list(list(objFunRes = objFunRes))#, nFires = NROW(a)))
        print(Sys.time())
        print(paste0("  ", Sys.getpid(), mess))
        if (SNLL_FSTest == failVal && ii == 1) {
          break
        }
      }
      if (isTRUE(doSNLLTest)) {
        rescalor <- 6e4
        SNLLTest <- sum(unlist(results$SNLL)) / rescalor
        mess <- paste(mess, " SNLLTest:", SNLLTest, "; ")
        objFunRes <- objFunRes + SNLLTest #+ SNLLTest

      }

    }
  }# run through 2nd batch of smaller fires

  bb <- purrr::transpose(objFunResList)
  bb <- purrr::map(bb, unlist)

  if (isTRUE(doMADTest) && !isTRUE(doSNLL_FSTest)) {
    totalNFires <- sum(bb$nFires)
    objFunRes <- do.call(sum, purrr::map2(bb$objFunRes, bb$nFires,
                                          function(.x, .y) .x * .y))/totalNFires
  }
  if (isTRUE(doSNLL_FSTest)) {
    objFunRes <- sum(unlist(bb$objFunRes))
    mess <- paste(" SNLL_FSTest total:", objFunRes, "; ", mess)
  }
  # Figure out what we want from these. This is potentially correct (i.e. we want the smallest ad.test and the smallest SNLL)
  return(objFunRes)
}

rescaleKnown <- function(x, minNew, maxNew, minOrig, maxOrig) {
  a1 <- x - minOrig # brings min to zero
  if (any(a1 != 0))
    a2 <- a1 * maxNew/max(a1)
  else
    a2 <- a1
  a2
}
