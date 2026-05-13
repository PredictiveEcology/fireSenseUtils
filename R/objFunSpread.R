utils::globalVariables(c(
  "..colsToKeep", "..colsToUse", ".N", "buffer", "burned", "burnedClass",
  "id", "ids", "initialLocus", "N", "numAvailPixels", "pixelID", "prob",
  "simFireSize", "size", "spreadProb", "value"
))

#' Objective function for `fireSense_spreadFit` module
#'
#' @param par parameters
#'
#' @param landscape A `SpatRaster` with extent, resolution, and projection (crs) used for
#'   [SpaDES.tools::spread2].
#'
#' @param annualDTx1000 A list of data.table class objects. Each list element is
#'   data from a single calendar year, and whose name takes the form `"yearxxxx"`,
#'   where `xxxx` is the four-digit year. The columns in the data.table must integers,
#'   that are `1000x` their actual values as this function will divide by 1000.
#'
#' @param nonAnnualDTx1000 Like `annualDTx1000`, but with where each list element will be
#'   used for >1 year. The names of the list elements must be of the form
#'   `"yearxxxx_yearyyyy_yearzzzz"` where `xxxx`, `yyyy`, and `zzzz` represent the four-digit
#'   calendar years for which that list element should be used.
#'   The columns are variables that are used for more than 1 year.
#'
#' @param formulaToFit Formula, put provided as a character string, not class `formula`.
#'   (if it is provided as a class `formula`, then it invariably will have an
#'   enormous amount of data hidden in the formula environment; this is bad for [DEoptim::DEoptim])
#'
#' @param historicalFires TODO: DESCRIPTION NEEDED
#'
#' @param fireBufferedListDT TODO: DESCRIPTION NEEDED
#'
#' @param covMinMax This is a 2 row by multiple column data.frame indicating
#'   the minimum and maximum values of the original covariate data values.
#'   These will be used to rescale the covariates internally so that they are all between 0 and 1.
#'   It is important to not simply rescale internally here because only 1 year is run at a time;
#'   all years must be rescaled for a given covariate by the same amount.
#'
#' @param maxFireSpread A value for `spreadProb` that is considered impossible to go above.
#'   Default 0.28, which is overly generous unless there are many non-flammable pixels (e.g., lakes).
#'
#' @param minFireSize TODO: DESCRIPTION NEEDED
#'
#' @template mutuallyExclusive
#'
#' @param doAssertions Logical. If `TRUE`, the default, the function will test a few minor things
#'   for consistency. This should be set to `FALSE` for operational situations, as the assertions
#'   take some small amount of time.
#'
#' @param tests One or more of `"mad"`, `"adTest"`, `"SNLL"`, or `"snll_fs"`. Default: `"snll_fs"`.
#'
#' @param Nreps Integer. The number of replicates, per ignition, to run.
#'
#' @param plot.it Passed to [SpaDES.core::Plots()], so will show (TRUE or "screen") or save files
#'   of several plots, using \pkg{ggplot2}.
#'
#' @param objFunCoresInternal Internally, this function can use `mcmapply` to run multiple
#'   parallel `spread` function calls. This should only be > `1L` if there are spare threads.
#'   It is highly likely that there won't be. However, sometimes the [DEoptim::DEoptim] is
#'   particularly inefficient, it starts X cores, and immediately several of them are
#'   stopped inside this function because the parameters are so bad, only two years are attempted.
#'   Then the core will stay idle until all other cores for the [DEoptim::DEoptim] iteration are complete.
#'   Similarly, if only physical cores are used for [DEoptim::DEoptim], the additional use of
#'   hyperthreaded cores here, internally will speed things up (i.e., this maybe could be `2L` or `3L`).
#'
#' @param thresh Threshold multiplier used in SNLL fire size (`"snll_fs"`) test. Default 550.
#'   Lowering the threshold value will be more restrictive, but being too restrictive will result
#'   in [DEoptim::DEoptim] rejecting more tests and using the "fail value" of 10000.
#'   Too high a threshold, and more years will be run and it will take longer to find values.
#'
#' @param lanscape1stQuantileThresh A `spreadProb` value that represents a threshold for the
#'   1st quantile of the `spreadProbs` on the landscape; if that quantile is above this
#'   number, then the `.objFunSpredFit` will bail because it is "too burny" a landscape.
#'   Default 0.265, meaning if only 25%% of the pixels on the landscape are below
#'   this `spreadProb`, then it will bail.
#'
#' @param weighted Logical. Should empirical likelihood be weighted by log of the actual fire size?
#'    This will give large fires more influence on the SNLL.
#'
#' @param verbose If >= 2, then this will show more information about `spreadProb` fitting.
#'
#' @param lowerSpreadProb Numeric. Lower bound for `spreadProb`; if a candidate
#'   `spreadProb` falls at or below this, the call exits early to avoid wasted
#'   simulation work. Default `0.13`.
#'
#' @param mutuallyExclusive Named list of vectors describing groups of model
#'   terms that must not be active together (mutually exclusive in the
#'   formula). Default `list("youngAge" = c("class", "nf"))`.
#'
#' @param ... This is not used here, but allows for extraneous arguments to not break this function.
#'
#' @return
#' Attempting a weighted likelihood,
#' <https://stats.stackexchange.com/questions/267464/algorithms-for-weighted-maximum-likelihood-parameter-estimation>.
#' With `log(fireSize) * likelihood` for each fire.
#'
#' @export
#' @importFrom data.table := rbindlist set setDT setDTthreads setnames setorderv
#' @importFrom EnvStats demp
#' @importFrom graphics abline axis hist mtext
#' @importFrom kSamples ad.test
#' @importFrom purrr map2 pmap
#' @importFrom quickPlot clearPlot dev gpar Plot
#' @importFrom terra buffer crop ext ncell rast trim xyFromCell
#' @importFrom sp SpatialPoints
#' @importFrom SpaDES.tools spread2
#' @importFrom stats as.formula dbinom median terms quantile
#' @importFrom utils tail
.objfunSpreadFit <- function(par,
                             landscape,
                             annualDTx1000,
                             nonAnnualDTx1000,
                             formulaToFit, # loci, sizes,
                             historicalFires,
                             fireBufferedListDT,
                             covMinMax = NULL,
                             maxFireSpread = 0.28, # 0.257 makes gigantic fires
                             lowerSpreadProb = 0.13,
                             minFireSize = 2,
                             tests = "snll_fs",
                             Nreps = 10,
                             mutuallyExclusive = list("youngAge" = c("class", "nf")),
                             doAssertions = TRUE,
                             plot.it = FALSE, # TODO parameterize this line for plotting
                             objFunCoresInternal = 1,
                             lanscape1stQuantileThresh = 0.265,
                             thresh = 550,
                             weighted = TRUE,
                             # bufferedRealHistoricalFiresList,
                             verbose = 2,
                             ...) { # fireSense_SpreadFitRaster
  # Optimization's objective function
  # lapply(historicalFires, setDT)

  data.table::setDTthreads(1)

  doMADTest <- any(grepl("mad", tolower(tests)))
  doSNLLTest <- any(grepl("snll$", tolower(tests)))
  doSNLL_FSTest <- any(grepl("snll_fs", tolower(tests)))
  doADTest <- any(grepl("adtest", tolower(tests)))
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
  # lapply(annualDTx1000, setDT)
  lapply(nonAnnualDTx1000, setDT)
  # lapply(fireBufferedListDT, setDT)
  # dtThreadsOrig <- data.table::setDTthreads(1)
  if (is(formulaToFit, "formula")) {
    stop("formulaToFit must be provided as a charater string because it takes too much RAM otherwise.")
  }
  formulaToFit <- as.formula(formulaToFit, env = .GlobalEnv)
  colsToUse <- attributes(terms(formulaToFit))[["term.labels"]]
  # How many of the parameters belong to the model?
  parsModel <- length(colsToUse)
  ncells <- ncell(landscape)

  r <- rast(landscape)
  years <- as.character(names(annualDTx1000))
  names(years) <- years
  cells <- integer(ncells)
  # Nreps <- 10
  yearSplit <- strsplit(names(nonAnnualDTx1000), "_")
  names(yearSplit) <- as.character(seq_along(nonAnnualDTx1000))
  indexNonAnnual <- rbindlist(
    Map(
      ind = seq_along(nonAnnualDTx1000), date = yearSplit,
      function(ind, date) data.table(ind = ind, date = date)
    )
  )
  historicalFiresAboveMin <- lapply(historicalFires, function(x) {
    x <- x[x$size >= minFireSize, ]
    x <- x[!duplicated(x$cells), ]
    x
  })

  ## can't fit fires for years with no data; drop these years
  omitYears <- names(historicalFiresAboveMin[which(lapply(historicalFiresAboveMin, nrow) == 0)])
  if (length(omitYears > 0)) {
    historicalFiresAboveMin[omitYears] <- NULL
  }

  fireSizesByYear <- unlist(lapply(historicalFiresAboveMin, function(x) sum(x$size)))
  largest <- head(sort(fireSizesByYear, decreasing = TRUE), 2) # max(2, objFunCoresInternal))
  smallest <- setdiff(names(fireSizesByYear), names(largest))
  lrgSmallFireYears <- list(large = names(largest), small = smallest)
  objFunResList <- list() # will hold objective function values --> which is now >1 for large, then small fires
  for (ii in seq(lrgSmallFireYears)) {
    yrs <- lrgSmallFireYears[[ii]]
    if (length(yrs)) {
      # results <- parallel::mcmapply(                             # normal
      # mc.cores = min(length(years[yrs]), objFunCoresInternal), # normal
      # mc.preschedule = FALSE,                                  # normal
      # SIMPLIFY = FALSE,                                        # normal
      results <- purrr::pmap( # interactive debugging
        .l = list( # interactive debugging
          annDTx1000 = annualDTx1000[yrs],
          yr = years[yrs],
          annualFires = historicalFiresAboveMin[yrs],
          annualFireBufferedDT = fireBufferedListDT[yrs] # interactive debugging
          # annualFireBufferedDT = fireBufferedListDT[yrs],            # normal
        ), # interactive debugging
        # MoreArgs = list(                                         # normal
        par = par, parsModel = parsModel,
        verbose = verbose,
        nonAnnualDTx1000 = nonAnnualDTx1000,
        indexNonAnnual = indexNonAnnual,
        colsToUse = colsToUse,
        mutuallyExclusive = mutuallyExclusive,
        doAssertions = doAssertions,
        maxFireSpread = maxFireSpread,
        lowerSpreadProb = lowerSpreadProb,
        lanscape1stQuantileThresh = lanscape1stQuantileThresh,
        Nreps = Nreps,
        plot.it = plot.it,
        r = r, weighted = weighted,
        doSNLL_FSTest = doSNLL_FSTest,
        doMADTest = doMADTest, doADTest = doADTest,
        cells = cells,
        covMinMax = covMinMax, # interactive debugging
        # covMinMax = covMinMax                              # normal
        # ),                                                   # normal
        .f = objFunInner # interactive debugging
        # objFunInner#(yr, annDTx1000, par, parsModel,             # normal
        #  annualFires, nonAnnualDTx1000, annualFireBufferedDT,
        #  indexNonAnnual, colsToUse, covMinMax,
        #  verbose = TRUE)
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
        # a[, devLog := sqrt(abs(size - simFireSize))]
        mad <- round(mean(a$dev), 1)
        # mad <- round(mean(a$devLog), 1)
        objFunRes <- objFunRes + mad #+ SNLLTest
        mess <- paste(" mad:", mad, "; ")
        objFunResList[ii] <- list(list(objFunRes = objFunRes, nFires = NROW(a)))
        # if (mad > 2700) {
        #   print(paste0("  ", Sys.getpid(), mess))
        #   break
        # }
      }

      if (isTRUE(doADTest)) {
        if (ii == 2) {
          historicalFiresTr <- unlist(purrr::transpose(historicalFiresAboveMin)$size)
          simulatedFires <- unlist(results$fireSizes)
          adTest <- try(ad.test(simulatedFires, historicalFiresTr)[["ad"]][1L, 1L])
          if (is(adTest, "try-error")) {
            adTest <- 1e6L
          }
          adTest <- adTest * 50
          objFunRes <- objFunRes + adTest #+ SNLLTest
          mess <- paste(mess, " adTest:", adTest, "; ")
        }
      }
      if (isTRUE(doSNLL_FSTest)) {
        thresh <- round(thresh, 0)
        SNLL_FSTest <- round(sum(unlist(results$SNLL)), 0)
        failVal <- 1e6L
        numYrsDone <- length(results$SNLL_FS)
        threshold <- thresh * numYrsDone ## lower is _more_ restrictive; too high takes too long
        mess <- character()
        annualSNLL <- round(SNLL_FSTest / numYrsDone, 0)
        if (any(SNLL_FSTest > threshold) && ii == 1) {
          SNLL_FSTestOrig <- SNLL_FSTest
          SNLL_FSTest <- failVal
          mess <- paste0(
            " Fail! Bailing after ", numYrsDone, " yrs; SNLL threshold: ", thresh, "; ",
            "Avg annual SNLL: ", annualSNLL, "; "
          )
        } else {
          if (ii == 1) {
            mess <- paste0(
              " Decent in 1st ", numYrsDone, " years -- continuing. ", mess, " SNLL threshold: ", thresh, ", Avg annual: ",
              annualSNLL, "; "
            )
          }
        }
        objFunRes <- objFunRes + SNLL_FSTest #+ SNLL_FSTest
        objFunResList[ii] <- list(list(objFunRes = objFunRes)) # , nFires = NROW(a)))
        if (length(mess) > 0) {
          print(paste0("  ", Sys.getpid(), mess))
        }
        if (SNLL_FSTest == failVal && ii == 1) {
          break
        }
      }
      if (isTRUE(doSNLLTest)) {
        rescalor <- 6e4
        SNLLTest <- sum(unlist(results$SNLL)) / rescalor
        mess <- paste(mess, " SNLLTest:", SNLLTest, "; ")
        objFunRes <- objFunRes + SNLLTest
      }
    }
  } # run through 2nd batch of smaller fires

  bb <- purrr::transpose(objFunResList)
  bb <- purrr::map(bb, unlist)

  if (isTRUE(doMADTest) && !isTRUE(doSNLL_FSTest)) {
    totalNFires <- sum(bb$nFires)
    objFunRes <- do.call(sum, purrr::map2(
      bb$objFunRes, bb$nFires,
      function(.x, .y) .x * .y
    )) / totalNFires
  }
  if (isTRUE(doSNLL_FSTest)) {
    objFunRes <- sum(unlist(bb$objFunRes))
  }
  if (length(objFunResList) > 1) {
    print(paste0(Sys.getpid(), "; FINISHED! ", Sys.time(), "; Objective Final: ", round(objFunRes, 0)))
  }
  ## Figure out what we want from these.
  ## This is potentially correct (i.e. we want the smallest ad.test and the smallest SNLL)
  return(objFunRes)
}

rescaleKnown <- function(x, minNew, maxNew, minOrig, maxOrig) {
  a1 <- x - minOrig # brings min to zero
  if (any(a1 != 0)) {
    a2 <- a1 * maxNew / max(a1)
  } else {
    a2 <- a1
  }
  a2
}

#' rescale function no.2
#'
#' @param x a vector to be rescaled
#' @param minNew the minimum of the new range
#' @param maxNew the max of the new range
#' @param minOrig the minimum of the original data
#' @param maxOrig the maximum of the original data
#' @return the rescaled vector
#'
#' @export
rescaleKnown2 <- function(x, minNew, maxNew, minOrig, maxOrig) {
  A <- maxOrig - minOrig # range of original
  b <- maxNew - minNew # range of new
  z <- b / A # ratio of range size
  C <- x - minOrig # make it have a new minimum above the minOrig
  D <- C * z
  return(D)
}

#' @keywords internal
#'
#' @importFrom data.table set setDT
#' @importFrom ggplot2 aes facet_wrap geom_histogram ggplot
#' @importFrom SpaDES.core anyPlotting Plots
#' @importFrom SpaDES.tools spread
#' @importFrom tidyr gather
objFunInner <- function(yr, annDTx1000, par, parsModel, # normal
                        annualFires, nonAnnualDTx1000, shortAnnDTx1000 = NULL,
                        annualFireBufferedDT,
                        indexNonAnnual, colsToUse, covMinMax, mutuallyExclusive,
                        doAssertions, maxFireSpread, lowerSpreadProb, cells, lanscape1stQuantileThresh,
                        weighted,
                        r, Nreps, doSNLL_FSTest, doMADTest, doADTest,
                        plot.it, verbose = 2) {
  if (isTRUE(plot.it)) plot.it <- "screen"

  # needed because data.table objects were recovered from disk
  # Rescale to numerics and /1000
  # setDT(nonAnnDTx1000)
  # matrix multiplication
  parsList <- paramsSeparate(par, parsModel)
  logisticPars <- parsList[["logisticPars"]]
  covPars <- parsList[["covPars"]]

  shortAnnDT <- spreadProbFromIntegerCovs(
    shortAnnDTx1000 = NULL, annDTx1000, nonAnnualDTx1000,
    indexNonAnnual, yr, covMinMax, mutuallyExclusive, colsToUse,
    doAssertions, logisticPars, covPars, maxFireSpread, lowerSpreadProb
  )

  set(shortAnnDT, NULL, "spreadProb",
      logisticAll(logisticPars, mat = as.matrix(shortAnnDT[, ..colsToUse]), covPars, lowerSpreadProb))
  doFitting <- any(c(doSNLL_FSTest, doMADTest, doADTest))
  if (isTRUE(doFitting)) {
    
    # shortAnnDTx1000 <- rescaleAllCovsFromX1000(annDTx1000 = annDTx1000,
    #                                   nonAnnualDTx1000 = nonAnnualDTx1000,
    #                                   indexNonAnnual = indexNonAnnual,
    #                                   yr = yr, covMinMax = covMinMax,
    #                                   mutuallyExclusive = mutuallyExclusive,
    #                                   colsToUse = colsToUse)
    # setDT(annDTx1000)
    # if (is.null(shortAnnDTx1000))
    #   shortAnnDTx1000 <- nonAnnualDTx1000[[indexNonAnnual[date == yr]$ind]][annDTx1000, on = "pixelID"]
    # if (!is.null(covMinMax)) {
    #   for (cn in colnames(covMinMax)) {
    #     set(
    #       shortAnnDTx1000, NULL, cn,
    #       rescaleKnown2(
    #         shortAnnDTx1000[[cn]], 0, 1000,
    #         covMinMax[[cn]][1] * 1000,
    #         covMinMax[[cn]][2] * 1000
    #       )
    #     )
    #   }
    # }
    # if (!is.null(mutuallyExclusive)) {
    #   shortAnnDTx1000 <- makeMutuallyExclusive(
    #     dt = shortAnnDTx1000,
    #     mutuallyExclusiveCols = mutuallyExclusive
    #   )
    # }
    # mat <- as.matrix(shortAnnDTx1000[, ..colsToUse]) / 1000
    # mat <- as.matrix(shortAnnDTx1000[, ..colsToUse])
    # if (doAssertions) {
    #   test1 <- sum(apply(round(mat[, colsToUse], 3), 2, min) < 0) == 0
    #   test2 <- sum(apply(round(mat[, colsToUse], 3), 2, max) > 1) == 0
    #   if (!all(test1, test2)) {
    #     stop("Covariates are not all between 0 and 1, which they should be")
    #   }
    # }
    # # matrix multiplication
    # parsList <- paramsSeparate(par, parsModel)
    # logisticPars <- parsList[["logisticPars"]]
    # covPars <- parsList[["covPars"]]
    #
    # # covPars <- tail(x = par, n = parsModel)
    # # logisticPars <- head(x = par, n = length(par) - parsModel)
    # if (logisticPars[1] > maxFireSpread) {
    #   warning(
    #     "The first parameter of the logistic is > ", maxFireSpread, ".",
    #     "The parameter should be lowered."
    #   )
    # }
    #
    # set(shortAnnDTx1000, NULL, "spreadProb", logisticAll(logisticPars, mat, covPars, lowerSpreadProb))
    
    if (SpaDES.core::anyPlotting(plot.it)) {
      par(
        mfrow = c(5, 6), omi = c(0.5, 0, 0, 0),
        mai = c(0.2, 0.3, 0.4, 0.1)
      )
      
      # suppressMessages(Require::Require("tidyr"))
      # a <- gather(as.data.frame(mat)) |> ggplot(aes(value)) +
      #   geom_histogram(bins = 10) +
      #   facet_wrap(~key, scales = 'free_y')
      # Plots(a)
      # Map(dat = as.data.frame(mat), nam = colnames(mat), function(dat, nam) {
      #   hist(dat, main = nam)})
    }
    
    # shortAnnDTx1000 <- logisticAll(logisticPars, shortAnnDTx1000, mat, covPars, lowerSpreadProb)
    # if (length(logisticPars) == 4) {
    #   stop("logistic with 4 parameters not tested yet")
    #   set(shortAnnDTx1000, NULL, "spreadProb", logistic4p(mat %*% covPars, logisticPars))
    # } else if (length(logisticPars) == 3) {
    #   set(shortAnnDTx1000, NULL, "spreadProb", logistic3p(mat %*% covPars, logisticPars, par1 = lowerSpreadProb))
    # } else if (length(logisticPars) == 2) {
    #   set(shortAnnDTx1000, NULL, "spreadProb", logistic2p(mat %*% covPars, logisticPars, par1 = lowerSpreadProb, par4 = 0.5))
    # }
    # logistic multiplication
    # set(annDTx1000, NULL, "spreadProb", logistic4p(annDTx1000$pred, par[1:4])) ## 5-parameters logistic
    # set(annDTx1000, NULL, "spreadProb", logistic5p(annDTx1000$pred, par[1:5])) ## 5-parameters logistic
    # actualBurnSP <- annDTx1000[annualFireBufferedDT, on = "pixelID"]
    medSP <- median(shortAnnDT$spreadProb, na.rm = TRUE)
    cells[as.integer(shortAnnDT$pixelID)] <- shortAnnDT$spreadProb
    
    nonEdgeValues <- cells[cells > (lowerSpreadProb * 1.025) | cells > (logisticPars[1] * 0.99)]
    sdSP <- diff(quantile(nonEdgeValues, c(0.1, 0.9)))
    if (is.na(sdSP)) sdSP <- 0
    
    medSPRight <- medSP <= maxFireSpread & medSP >= lowerSpreadProb
    spreadOutEnough <- sdSP / medSP > 0.025
    ret <- list()
    minLik <- 1e-29 # min(emp$lik[emp$lik > 0])
    loci <- annualFires$cells
    summ <- summary(nonEdgeValues)
    lowSPLowEnough <- summ[2] < lanscape1stQuantileThresh
    
    if (verbose > 1) {
      if (isTRUE(!spreadOutEnough)) {
        print(paste0(
          "  ",
          Sys.getpid(), " FAIL! ", yr, "; Not spread out enough; bailing: ",
          paste(names(summ), round(summ, 3), collapse = ", ")
        ))
      }
      if (isTRUE(!lowSPLowEnough)) {
        print(paste0(
          "  ",
          Sys.getpid(), " FAIL! ", yr, "; Too burny a landscape; bailing: ",
          paste(names(summ), round(summ, 3), collapse = ", ")
        ))
      }
    }
    
    if (medSPRight && spreadOutEnough && lowSPLowEnough) {
      if (verbose > 1) {
        ww <- if (isTRUE(weighted)) "weighted" else "unweighted"
        print(paste0(
          " ",
          Sys.getpid(), ": ", yr, ", ", ww, ", spreadProbs: ",
          paste(names(summ), round(summ, 3), collapse = ", ")
        ))
      }
      # maxSizes <- rep(annualFires$size, times = Nreps)
      
      # this will make maxSizes be a little bit larger for large fires, but a lot bigger for small fires
      # maxSizes <- maxSizes * 1.5#(1.1+pmax(0,5-log10(maxSizes)))
      setDT(annualFireBufferedDT)
      minSize <- 100
      if (doAssertions || SpaDES.core::anyPlotting(plot.it)) {
        tableOfBufferedMaps <- annualFireBufferedDT[, list(numAvailPixels = .N), by = "ids"]
        tableOfBufferedMaps <- tableOfBufferedMaps[annualFires, on = "ids"]
        setnames(tableOfBufferedMaps, old = "cells", new = "initialLocus")
        minSizes <- tableOfBufferedMaps$numAvailPixels
        minSize <- quantile(minSizes, 0.3)
        if (minSize < 2000) {
          warning(
            "The fireSizeBufferDT has too many fires < 2000 burned + unburned pixels;",
            " needs larger buffers."
          )
        }
      }
      maxSizes <- fireSenseUtils::multiplier(annualFires$size, minSize = minSize)
      # maxSizes <- annualFires$size * 2
      # if (any(cells[loci] == 0)) {
      cells[loci] <- 1
      # }
      dups <- duplicated(annualFires$cells)
      if (any(dups)) {
        annualFires <- annualFires[which(!dups), ] #
        maxSizes <- maxSizes[!dups]
        loci <- annualFires$cells[!dups]
      }
      st <- system.time(spreadState <- lapply(seq_len(Nreps), function(i) {
        SpaDES.tools::spread(
          # SpaDES.tools::spread2(
          landscape = r,
          maxSize = maxSizes,
          # start = loci,
          loci = loci,
          spreadProb = cells,
          # asRaster = FALSE,
          returnIndices = TRUE,
          allowOverlap = FALSE,
          skipChecks = TRUE
        )
      }))
      if (SpaDES.core::anyPlotting(plot.it)) {
        # par(
        #   mfrow = c(7, 7), omi = c(0.5, 0, 0, 0),
        #   mai = c(0.2, 0.3, 0.4, 0.1)
        # )
        
        lapply(colsToUse, function(cn) hist(shortAnnDT[[cn]], main = cn))
        mtext(side = 3, "Histograms of distribution of rescaled variables", outer = TRUE, line = -1)
        hist(nonEdgeValues, main = "spreadProb")
        sam <- sample(NROW(shortAnnDT), NROW(shortAnnDT) / 100)
        val <- as.matrix(shortAnnDT[sam, ..colsToUse]) %*% covPars
        # val <- mat[sam, ] %*% covPars
        plot(val, shortAnnDT$spreadProb[sam],
             pch = ".",
             main = paste0("logits: ", paste(round(logisticPars, 2), collapse = ", "))
        )
      }
      
      spreadState <- rbindlist(spreadState, idcol = "rep")
      if (isTRUE(doSNLL_FSTest)) {
        emp <- spreadState[, list(N = .N), by = c("rep", "initialLocus")] # N is "size of simulated fire"
        emp <- emp[annualFires, on = c("initialLocus" = "cells")]
        if (SpaDES.core::anyPlotting(plot.it)) {
          colsToKeep <- c(setdiff(colnames(tableOfBufferedMaps), colnames(emp)), "initialLocus")
          emp <- tableOfBufferedMaps[, ..colsToKeep][emp, on = c("initialLocus"), nomatch = NULL]
          maxX <- log(max(c(annualFires$size, emp$N, emp$numAvailPixels)))
          emp <- setorderv(emp, c("size"), order = -1L)
          numLargest <- 4
          numHists <- 49 - numLargest - length(par) - 12 - 1 # 12 for rasters
          uniqueEmpIds <- unique(emp$initialLocus)
          sam <- # if (length(uniqueEmpIds) >= (numHists)) {
            # try(c(
            unique(emp$ids)[1:numLargest]#,
          #sample(unique(emp$initialLocus)[-(1:numLargest)],
          #  size = min(length(unique(emp$initialLocus)) - numLargest, numHists)
          #)
          #   ))
          # } else {
          #   uniqueEmpIds
          # }
          emp[ids %in% sam,
              {
                dat <- round(log(N))
                h <- hist(dat,
                          breaks = 30,
                          main = paste(as.character(.BY)),
                          # main = "",
                          axes = FALSE, xlim = c(0, maxX)
                )
                seqq <- seq(0, ceiling(maxX), by = 1)
                axis(1, at = seqq, labels = round(exp(seqq), 0))
                abline(v = log(size[1]), col = "red")
                abline(v = log(unique(numAvailPixels)), col = "green")
              },
              by = "initialLocus"
          ]
          mtext(
            outer = TRUE,
            paste(
              yr, "; sample of fires (incl. 4 largest);",
              "Simulated fire sizes (# pixels);",
              "Actual Fire (red);",
              "Available pixels to burn (green - should be well right of hist bars);",
              "Sorted by actual fire size."
            ),
            line = 2, side = 1
          )
        }
        
        # only use fires that escaped --> i.e., greater than 1 pixel
        # print(quantile(emp$size))
        ## TODO: occasional errors during fitting because `obs` < 2; suggests N <2 (#8)
        # Eliot: If this errors, it means that only one simulated fire had 2 pixels (i.e., N>1); if only one
        #   fire has 2+ pixels, need more fires, or call it a fail
        emp <- emp[N > 1, list(size = size[1],
                               lik = if(.N > 1) {
                                 EnvStats::demp(x = size[1], obs = sqrt(N))
                               } else {
                                 0
                               }), by = "ids"]
        # emp <- emp[N > 1, list(size = size[1],
        #                            lik = ifelse(.N > 1,
        #                                         EnvStats::demp(x = size[1], obs = sqrt(N)),
        #                                         0)), by = "ids"]
        # emp <- emp[N > 1, list(size = size[1], lik = EnvStats::demp(x = size[1], obs = N)), by = "ids"]
        if (isTRUE(weighted)) {
          set(emp, NULL, "lik", log(pmax(minLik, emp$lik * log(emp$size))))
        } else {
          set(emp, NULL, "lik", log(pmax(minLik, emp$lik)))
        }
        
        # set(emp3, NULL, "lik", log(pmax(minLik, emp3$lik)))
        SNLL_FS <- -sum(emp$lik)
        # SNLL_FS3 <- -sum(emp3$lik)
        # print(paste0("Sqrt: ", round(SNLL_FS, 0), ", Normal: ", round(SNLL_FS3, 0)))
        ret <- append(ret, list(SNLL_FS = SNLL_FS))
      }
      
      if (isTRUE(doMADTest) || isTRUE(doADTest)) {
        fireSizes <- round(spreadState[, .N, .(initialLocus)][["N"]] / Nreps, 0) # Here tabulate() is equivalent to table() but faster
        ret <- append(ret, list(fireSizes = fireSizes))
      }
      if (SpaDES.core::anyPlotting(plot.it)) { # THIS IS PLOTTING STUFF
        # if (isTRUE(doSNLLTest)) {
        
        burnedProb <- spreadState[, .N, by = "indices"]
        setnames(burnedProb, "indices", "pixelID")
        out <- burnedProb[annualFireBufferedDT, on = "pixelID"]
        
        # fix the out
        # 1 -- set pixels that had not simulated fires to N = 0
        out[is.na(N), N := 0]
        # 2 -- rescale probability surface between 0.001 and 0.99
        #      so probabilities can be calculated
        out[, prob := pmin(out$N / Nreps + 0.001, 0.99)]
        # 3 -- convert buffer (which has 1 in buffer) to burned = 1 - buffer
        out[, burned := buffer]
        # 4 -- Set initial pixels to burned = 2 -- is a work around for cases where "initial pixels" are not actually burned in
        #   the polygon database
        out[, burnedClass := burned]
        out[pixelID %in% annualFires$cell, burnedClass := 2]
        bigFire1 <- rast(r)
        bigFire1[out$pixelID] <- out$ids
        keepFire <- tail(sort(table(out$ids)), 4)
        setDT(annualFires)
        theseFires <- annualFires[ids %in% names(keepFire)]
        # clearPlot()
        firesToDo <- theseFires$ids
        names(firesToDo) <- firesToDo
        out2 <- lapply(firesToDo, function(id) {
          keepFire <- as.numeric(id)
          # keepFire <- 65
          bigFire <- bigFire1
          bigFire[bigFire != keepFire] <- NA
          bf <- trim(bigFire)
          ex <- ext(bf)
          
          thisFire <- annualFires[ids == keepFire]
          
          r <- rast(r)
          r[out$pixelID] <- out$prob
          
          predictedFireProb <- crop(r, ex)
          # clearPlot();Plot(r)
          actualFire <- rast(r)
          actualFire[out$pixelID] <- out$burnedClass
          actualFire <- crop(actualFire, ex)
          levels(actualFire) <- data.frame(ID = 0:2, class = c("unburned", "burned", "ignited"))
          
          predictedLiklihood <- dbinom(
            prob = out$prob,
            size = 1,
            x = out$burned,
            log = TRUE
          )
          spreadProbMap <- rast(r)
          spreadProbMap[out$pixelID] <- cells[out$pixelID]
          spreadProbMap <- crop(spreadProbMap, ex)
          spreadProbMap[spreadProbMap >= par[1]] <- par[1]
          ccc <- cells[out$pixelID]
          ccc <- ccc[ccc > 0]
          lowerLim <- quantile(ccc, 0.05)
          ccc <- ccc[ccc > lowerLim]
          spreadProbMap[spreadProbMap <= lowerLim] <- lowerLim
          predLiklihood <- rast(r)
          predLiklihood[out$pixelID] <- predictedLiklihood
          predLiklihood <- crop(predLiklihood, ex)
          spIgnits <- terra::vect(xyFromCell(r, thisFire$cells))
          spIgnits <- buffer(spIgnits, width = 5000)
          spIgnits <- crop(spIgnits, ex)
          list(
            spIgnits = spIgnits, predictedFireProb = predictedFireProb,
            predLiklihood = predLiklihood,
            spreadProbMap = spreadProbMap
          )
        })
        out3 <- purrr::transpose(out2)
        notSp <- grep("spIgnits", names(out3), value = TRUE, invert = TRUE)
        
        out4 <- unlist(out3[notSp], recursive = FALSE)
        lapply(out4, function(x) terra::plot(x, col = RColorBrewer::brewer.pal(9, "Paired")))
        # clearPlot()
        # clearPlot()
        # a <- Plot(out4, cols = "Paired", new = TRUE, visualSqueeze = 0.85)
        # nn <- lapply(names(out3$spIgnits), function(id) {
        #   spDat <- out3$spIgnits[[id]]
        #   Plot(spDat,
        #        addTo = grep(id, names(out4), value = TRUE)[2],
        #        gp = gpar(fill = rep("transparent", 10), col = "black"), title = ""
        #   )
        # })
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
        # }
        # Add a very small number so that no pixel has exactly zero probability -- creating Inf
        # SNLL <- -sum(dbinom(prob = out$prob,
        #                     size = 1,
        #                     x = out$burned,
        #                     log = TRUE
        # ), na.rm = TRUE) # Sum of the negative log likelihood
      }
    } else {
      llik <- rep(log(minLik), length(loci))
      SNLL_FS <- -sum(llik)
      ret <- append(ret, list(SNLL_FS = SNLL_FS))
      # stop("encountered error with spreadProb - contact module developers")
      # Ian added this stop. Unclear what is supposed to happen. Object ret doesn't exist
      # SNLL <- 1e7
      # fireSizes <- sample(1:3, 1)
    }
  }

  return(ret)
}

#' Convert covariates from their `x1000` integer to usable by spread
#'
#' @inheritParams .objfunSpreadFit
#'
#' @param shortAnnDTx1000 Optional if annDTx1000 and nonAnnualDTx1000 are supplied. Otherwise
#'   it must be a data.frame/data.table that has all covariates (>= 1 forests, >=0 non-forests, 
#'   climate, youngAge), whose values are 1000x their original values so they can be stored as
#'   integers.
#'
#' @param annDTx1000 TODO: use description of `annualDTx1000` parameter in `.objfunSpreadFit`.
#'   If `shortAnnDTx1000` is not supplied, this must be supplied.
#'   It must be a data.frame/data.table that has all the annual covariates (e.g., `>=1` forest biomass, 
#'   `>=0` non forest binary), whose values are 1000x their original values so they can be stored
#'   as integers.
#' @param nonAnnualDTx1000 If `shortAnnDTx1000` is not supplied, this must be supplied.
#'   A named list of data.frames/data.tables that has all the non-annual covariates (e.g., `>=1` forest biomass, 
#'   `>=0` non forest binary), whose values are 1000x their original values so they can be stored
#'   as integers. These non annual covariates must have a `date` column that can be assessed
#'   against `yr` (via `nonAnnualDTx1000[[whKeep]][annDTx1000, on = "pixelID"]`)
#'
#' @param indexNonAnnual `data.table` with 2 columns: `ind` and `date` which links the
#'   list elements in `nonAnnualDTx1000` with their date (in case `nonAnnualDTx1000` is 
#'   not a named list. This is not necessary if the `nonAnnualDTx1000` is named with 
#'   equivalent convention (e.g., character or numeric with year) as `yr`.
#'
#' @param yr A single character or numeric/integer representing the full year (e.g., 2020) to
#'   but used. This year will be extracted from both the `annDTx1000` and `nonAnnualDTx1000` if
#'   they are supplied
#' @param covMinMax A data.table, one column for each column in `shortAnnDTx1000` (or the
#'   ann and nonAnnual alternatives), where the two rows represent the minimum and maximum
#'   values in the original fitting dataset. This MUST be supplied if this is prediction
#'   scenario so that the new covariates values are not scaled to themselves. I.e., they
#'   must be rescaled compared to the fitted data or else their rescaled values will be 
#'   incorrect.
#'
#' @template mutuallyExclusive
#'
#' @param colsToUse Optional. If this is supplied, it must be a character vector indicating
#'   the column names to use in `shortAnnDTx1000`, i.e., it must include everything, 
#'   annual or nonAnnual, that will be used.
#'
#' @inheritParams logisticAll
#'
#' @return This returns the full data.table to be used for fireSense, on the original numeric 
#'   i.e., real scale, with the `x1000` reversed.
#'
#' @export
#' @importFrom data.table set setDT
spreadProbFromIntegerCovs <- function(shortAnnDTx1000 = NULL, annDTx1000, nonAnnualDTx1000,
                                      indexNonAnnual, yr, covMinMax = NULL, mutuallyExclusive, colsToUse,
                                      doAssertions, logisticPars, covPars, maxFireSpread,
                                      lowerSpreadProb) {
  # rescaleA <- function(annDTx1000, shortAnnDTx1000, nonAnnualDTx1000, indexNonAnnual,
  #                      yr, covMinMax, mutuallyExclusive, colsToUse, doAssertions, logisticPars,
  #                      maxFireSpread, covPars, lowerSpreadProb) {

  if (!missing(annDTx1000)) {
    setDT(annDTx1000)
  }
  if (is.character(yr)) {
    yr <- gsub("[[:alpha:]]+", "", yr) |> as.integer()
  }
  if (is.null(shortAnnDTx1000)) {
    if (!is.null(names(nonAnnualDTx1000))) {
      whElement <- which(as.integer(names(nonAnnualDTx1000)) <= yr)
    } 
    
    if (is.numeric(whElement) && all(!is.na(whElement))) {
      whKeep <- tail(whElement, 1)
    } else {
      whKeep <- tail(indexNonAnnual$ind[indexNonAnnual$date <= yr], 1)
    }
    shortAnnDTx1000 <- nonAnnualDTx1000[[whKeep]][annDTx1000, on = "pixelID"]
  }
  if (!is.null(covMinMax)) {
    for (cn in colnames(covMinMax)) {
      set(
        shortAnnDTx1000, NULL, cn,
        rescaleKnown2(
          shortAnnDTx1000[[cn]], 0, 1000,
          covMinMax[[cn]][1] * 1000,
          covMinMax[[cn]][2] * 1000
        )
      )
    }
  }
  if (!is.null(mutuallyExclusive)) {
    shortAnnDTx1000 <- makeMutuallyExclusive(
      dt = shortAnnDTx1000,
      mutuallyExclusiveCols = mutuallyExclusive
    )
  }
  # mat <- as.matrix(shortAnnDTx1000[, ..colsToUse]) / 1000
  for (cn2 in colsToUse) {
    set(shortAnnDTx1000, NULL, cn2, shortAnnDTx1000[[cn2]] / 1000)
  }
  
  # rename because it is no longer x1000
  shortAnnDT <- shortAnnDTx1000

  if (doAssertions) {
    test1 <- sum(apply(round(shortAnnDT[, ..colsToUse], 3), 2, min) < 0) == 0
    test2 <- sum(apply(round(shortAnnDT[, ..colsToUse], 3), 2, max) > 1) == 0

    # test1 <- sum(apply(round(mat[, colsToUse], 3), 2, min) < 0) == 0
    # test2 <- sum(apply(round(mat[, colsToUse], 3), 2, max) > 1) == 0
    if (!all(test1, test2)) {
      stop("Covariates are not all between 0 and 1, which they should be")
    }
    if (logisticPars[1] > maxFireSpread) {
      warning(
        "The first parameter of the logistic is > ", maxFireSpread, ".",
        "The parameter should be lowered."
      )
    }
  }

  # covPars <- tail(x = par, n = parsModel)
  # logisticPars <- head(x = par, n = length(par) - parsModel)
  shortAnnDT
}
