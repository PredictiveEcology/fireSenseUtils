fireSpreadSDTxt <- "fireSpreadSD"

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
#' @param historicalFires A named list of `data.frame`s (one per year, names
#'   matching those of `annualDTx1000`), each with columns `cells` (pixel
#'   indices of ignitions) and `size` (fire size in pixels). The observed
#'   fire record against which the simulation is compared.
#'
#' @param fireBufferedListDT A named list of `data.table`s (one per year,
#'   names matching `annualDTx1000`), each with columns `pixelID`, `buffer`
#'   (`1` if the pixel is inside a buffered fire footprint, `0` otherwise),
#'   and `ids` (a unique identifier per fire). Restricts the set of pixels
#'   used when fitting and scoring spread.
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
#' @param minFireSize Integer. Minimum fire size (in pixels) to include when
#'   scoring simulated against historical fires; fires smaller than this are
#'   filtered out of the comparison. Default `1`.
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
#' @param covCentre Optional named list or vector of values subtracted from the rescaled covariates
#'   of the same name. `NULL` (default) does no centring. Must be the same at fit and predict.
#' @param pruneAbove Numeric. An upper bound on the accumulated SNLL after the *first* batch of
#'   fire years, past which the evaluation stops early and returns the fail value. Default `Inf`,
#'   which leaves `thresh * <years done>` as the only bound -- i.e. no change in behaviour.
#'
#'   A caller that runs one DEoptim generation per call knows the current population, and can pass
#'   the worst objective value that population would accept. Doing so is *exact*, not a heuristic:
#'   the second batch contributes a non-negative SNLL, so the first batch's accumulated value is a
#'   lower bound on the evaluation's final value. A trial exceeding `pruneAbove` would therefore
#'   have finished at or above a value every parent already beats, and selection would have
#'   discarded it. Pruning cannot change which trials DEoptim keeps; it only avoids finishing an
#'   evaluation whose outcome is already decided. This matters because a generation is
#'   synchronous -- its wall time is the slowest of its `NP` evaluations, not the median one.
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
#' @param weighted Should each fire's log-likelihood be multiplied by a weight that grows with its
#'   observed size, to give large fires more influence on the SNLL? `FALSE`: no. `TRUE` or `"log"`:
#'   `log(size)`. `"sqrt"`: `sqrt(size)`. Weights are divided by their mean over the fitted fires, so
#'   the SNLL keeps its scale against the other tests.
#'
#' @param sizeLik How a fire's likelihood is taken from its `Nreps` simulated sizes, on the
#'   square-root scale. `"kde"` (default): their kernel density at the observed size
#'   ([EnvStats::demp()]); it is zero away from the simulated sizes, so a fire they never reach
#'   scores the `minLik` floor however far off it is. `"t"`: a Student-t centred on their mean with
#'   their standard deviation as scale; no floor, and a miss costs more the further it is.
#'
#' @param sizeLikDf Degrees of freedom of the Student-t when `sizeLik = "t"`. Default 5. Lower
#'   values forgive a missed fire more.
#'
#' @param adWeight What the `"adTest"` statistic is multiplied by before it is added to the fire-size
#'   SNLL. `"auto"` (default) uses [adWeightAuto()], which gives the two terms about equal influence
#'   for any `sizeLik` and `weighted`. A number is used as given; the former fixed default, `50`,
#'   balanced the two only for `sizeLik = "kde"` with `weighted = FALSE`.
#'
#' @param link Passed to [logisticAll()]: `NULL` for the form implied by the logistic parameters, or
#'   `"logistic3pUpper"` for the 3-parameter logistic with Stukel's upper tail, whose fourth logistic
#'   parameter is `upperTail1`. Pass it explicitly when fitting: `par` may arrive unnamed.
#'
#' @param returnSims If `TRUE`, return the simulated fires instead of the objective value: a
#'   `data.table` with one row per fire year x fire x replicate -- `yr`, `rep`, `initialLocus`, `ids`,
#'   the simulated size `sim` and the observed `size`, in pixels. A year the objective declines to
#'   simulate (its spread probabilities fail the sanity checks) keeps its fires with `rep` and `sim`
#'   `NA`, so the caller can count them. Every year is simulated: the early stop (`thresh`) is not
#'   applied. For diagnostics and validation; the value is not an objective.
#'   The result carries `attr(, "spreadProb")`: the spread probabilities of every fitted pixel-year,
#'   summarised by [spreadProbSummary()], for [linkSaturation()].
#' @param fitFireSpreadSD If `TRUE`, the last element of `par` is the sd of a per-fire random effect on
#'   the logit of spread probability: in each replicate every fire draws `eps ~ N(0, sd^2)` and every pixel
#'   of its buffer spreads with `plogis(qlogis(p) + eps)`. Weather during a fire is not in the model, so
#'   without it every fire of a year sees identical probabilities and the simulated fire sizes are too
#'   alike: medians too large and the largest fires too small at once (2026-09-23). `sd = 0` is the model
#'   without it, exactly (no random numbers are drawn). `NULL` (default): `TRUE` when `par` is named and
#'   includes `fireSpreadSD`. DEoptim passes `par` unnamed, so `runDEoptim()` sets it explicitly.
#' @param capSizes If `FALSE`, simulated fires are not capped at [multiplier()] of their observed
#'   size. The fit needs the cap; validation should not have it, or a model that predicts fires far
#'   too large cannot show it.
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
                             pruneAbove = Inf,
                             covCentre = NULL,
                             weighted = TRUE,
                             sizeLik = "kde",
                             sizeLikDf = 5,
                             adWeight = "auto",
                             link = NULL,
                             returnSims = FALSE,
                             capSizes = TRUE,
                             fitFireSpreadSD = NULL,
                             # bufferedRealHistoricalFiresList,
                             verbose = 2,
                             ...) { # fireSense_SpreadFitRaster
  # Optimization's objective function
  # lapply(historicalFires, setDT)

  data.table::setDTthreads(1)
  sizeLik <- match.arg(sizeLik, c("kde", "t"))
  ## the per-fire random effect's sd, fitted as the LAST element of `par`
  fireSpreadSD <- 0
  if (is.null(fitFireSpreadSD)) fitFireSpreadSD <- fireSpreadSDTxt %in% names(par)
  if (isTRUE(fitFireSpreadSD)) {
    if (!is.null(names(par)) && !identical(names(par)[length(par)], fireSpreadSDTxt))
      stop("`", fireSpreadSDTxt, "` must be the last parameter")
    fireSpreadSD <- unname(par[length(par)])
    par <- par[-length(par)]
  }

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
  ## numeric, not integer: it receives spreadProb (double), and assigning doubles into an
  ## integer vector coerces the whole landscape-length vector on every fire year.
  cells <- numeric(ncells)
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

  sizeWeightMean <- mean(sizeWeight(unlist(lapply(historicalFiresAboveMin, function(x) x$size)), weighted))

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
  fireSizesList <- list() # simulated fire sizes, pooled across batches for the adTest
  yrsDoneList <- list() # the years contributing to fireSizesList, kept in lockstep
  bailedEarly <- FALSE
  simsList <- list() # returnSims: the simulated fires of each batch of years
  pSumList <- list() #   ... and the spread probabilities of their pixels
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
        r = r, weighted = weighted, sizeWeightMean = sizeWeightMean,
        doSNLL_FSTest = doSNLL_FSTest,
        doMADTest = doMADTest, doADTest = doADTest,
        cells = cells,
        covCentre = covCentre,
        sizeLik = sizeLik, sizeLikDf = sizeLikDf, link = link,
        returnSims = returnSims, capSizes = capSizes, fireSpreadSD = fireSpreadSD,
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
      if (isTRUE(returnSims)) {
        simsList[[ii]] <- data.table::rbindlist(results$sims, use.names = TRUE)
        pSumList[[ii]] <- results$pSummary
        next
      }

      if (isTRUE(doADTest)) {
        fireSizesList[[ii]] <- unlist(results$allFireSizes)
        yrsDoneList[[ii]] <- yrs
      }

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

      if (isTRUE(doSNLL_FSTest)) {
        thresh <- round(thresh, 0)
        SNLL_FSTest <- round(sum(unlist(results$SNLL)), 0)
        failVal <- 1e6L
        numYrsDone <- length(results$SNLL_FS)
        ## lower is _more_ restrictive; too high takes too long. `pruneAbove` (default Inf, i.e. no
        ## effect) lets the caller tighten this with the worst value the current population would
        ## accept: block 2's SNLL is non-negative, so this block's value is a lower bound on the
        ## total, and a trial above `pruneAbove` is one every parent already beats.
        threshold <- min(thresh * numYrsDone, pruneAbove)
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
          bailedEarly <- TRUE
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

  if (isTRUE(returnSims)) {
    out <- data.table::rbindlist(simsList, use.names = TRUE)
    data.table::setattr(out, "spreadProb", combineSpreadProbSummaries(unlist(pSumList, recursive = FALSE)))
    return(out)
  }
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
  ## The Anderson-Darling test compares whole distributions, so it runs once, after
  ## all batches, on the simulated fire sizes pooled across batches -- and against the
  ## observed sizes from *those same years*. Previously it ran inside the loop under
  ## `ii == 2`, which pooled only the final batch's simulated sizes while comparing
  ## them to every year's observed sizes; the omitted years are the 2 with the largest
  ## area burned, so the observed sample kept an upper tail the simulated sample could
  ## not have.
  if (isTRUE(doADTest) && !isTRUE(bailedEarly)) {
    pooled <- pooledFireSizes(fireSizesList, yrsDoneList, historicalFiresAboveMin)
    adTest <- try(adStatistic(pooled$simulated, pooled$observed))
    if (is(adTest, "try-error")) {
      adTest <- 1e6L
    }
    adTest <- adTest * if (identical(adWeight, "auto")) {
      adWeightAuto(sum(vapply(historicalFiresAboveMin, nrow, integer(1))), sizeLik, weighted)
    } else {
      adWeight
    }
    objFunRes <- objFunRes + adTest
    if (verbose > 1) {
      print(paste0("  ", Sys.getpid(), " adTest:", adTest, "; "))
    }
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

#' Pool simulated and observed fire sizes over the same years
#'
#' The Anderson-Darling test in `objFunSpread()` compares whole distributions of
#' fire size, so both samples must be drawn from the *same* set of years. The years
#' are simulated in batches (largest-area years first, so a hopeless parameter set
#' can bail early), and this assembles the two samples from the per-batch
#' accumulators, skipping any batch that did not run.
#'
#' @param fireSizesList List with one element per batch, each a numeric vector of
#'   simulated fire sizes. Elements for batches that did not run are `NULL`.
#' @param yrsDoneList List with one element per batch, each the year names
#'   contributing to the matching element of `fireSizesList`.
#' @param historicalFiresAboveMin Named list of observed fires, one element per
#'   year, each with a `size` column.
#'
#' @return A list of two numeric vectors, `simulated` and `observed`, covering the
#'   same years.
#'
#' @keywords internal
#' @rdname pooledFireSizes
pooledFireSizes <- function(fireSizesList, yrsDoneList, historicalFiresAboveMin) {
  yrsDone <- unlist(yrsDoneList)
  list(
    simulated = unlist(fireSizesList),
    observed = unlist(purrr::transpose(historicalFiresAboveMin[yrsDone])$size)
  )
}

#' Two-sample Anderson-Darling statistic
#'
#' The statistic `kSamples::ad.test(x, y)$ad[1, 1]` reports (version 1, ties allowed), and nothing
#' else. `ad.test()` also standardises it, which takes a loop over every observation: 4 s for the
#' 44,000 simulated and 1,800 observed fires of one evaluation on ELF 4.3, quadratic in the number
#' of fires. The objective uses only the statistic.
#'
#' @param x,y Numeric vectors: the two samples.
#'
#' @return A single numeric.
#' @keywords internal
adStatistic <- function(x, y) {
  ## as ad.test() does: when every year bailed there are no simulated fires, and 0 / 0 below is NaN
  stopifnot(length(x) > 0, length(y) > 0)
  N <- as.numeric(length(x) + length(y)) # numeric: N * M below overflows integers past ~46,000 fires
  z <- sort(unique(c(x, y)))
  l <- tabulate(match(c(x, y), z), length(z)) # multiplicity of each distinct value in the pooled sample
  B <- cumsum(as.numeric(l))
  j <- seq_len(length(z) - 1L)
  sum(vapply(list(x, y), function(s) {
    M <- cumsum(as.numeric(tabulate(match(s, z), length(z))))
    sum((l / N * (N * M - length(s) * B)^2 / (B * (N - B)))[j]) / length(s)
  }, numeric(1)))
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
                        plot.it, verbose = 2, covCentre = NULL, sizeLik = "kde", sizeLikDf = 5,
                        sizeWeightMean = 1, link = NULL, returnSims = FALSE, capSizes = TRUE,
                        fireSpreadSD = 0) {
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
    doAssertions, logisticPars, covPars, maxFireSpread, lowerSpreadProb, covCentre = covCentre
  )

  set(shortAnnDT, NULL, "spreadProb",
      logisticAll(logisticPars, mat = as.matrix(shortAnnDT[, ..colsToUse]), covPars, lowerSpreadProb,
                  link = link))
  pSummary <- if (isTRUE(returnSims)) spreadProbSummary(shortAnnDT$spreadProb, ceiling = logisticPars[1])
  ## Initialised here, not inside the branch below, because the function returns
  ## it unconditionally. With no test selected -- which is how the module's
  ## `debug` mode calls this, passing tests = "" -- doFitting is FALSE, the branch
  ## is skipped, and `return(ret)` used to fail with "object 'ret' not found".
  ## Nothing appends to it before the branch, so this is the same list it always
  ## was; the only change is that it exists when no test asked for anything.
  ret <- list()
  doFitting <- any(c(doSNLL_FSTest, doMADTest, doADTest)) || isTRUE(returnSims)
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
    ## Taken from the spreadProb column, not by scanning `cells`. `cells` is landscape-length and
    ## zero everywhere except this year's pixels (6.0M cells vs a median 41k pixels on ELF 5.3.1),
    ## so `cells[cells > a | cells > b]` made four passes over the landscape to recover values that
    ## are already here. Same values: pixelID is unique within a year, the zeros never pass a
    ## non-negative threshold, and `x > a | x > b` is `x > min(a, b)`. Only quantile() and summary()
    ## read it, so order does not matter. Measured on ELF 5.3.1 with identical seeds: identical
    ## objective values, and an evaluation 1.1-1.4x faster.
    nonEdgeValues <- shortAnnDT$spreadProb[
      shortAnnDT$spreadProb > min(lowerSpreadProb * 1.025, logisticPars[1] * 0.99)]
    sdSP <- diff(quantile(nonEdgeValues, c(0.1, 0.9)))
    if (is.na(sdSP)) sdSP <- 0
    
    medSPRight <- medSP <= maxFireSpread & medSP >= lowerSpreadProb
    spreadOutEnough <- sdSP / medSP > 0.025
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
        ww <- if (isFALSE(weighted)) "unweighted" else "weighted"
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
      maxSizes <- if (isTRUE(capSizes)) {
        fireSenseUtils::multiplier(annualFires$size, minSize = minSize)
      } else {
        rep(Inf, NROW(annualFires))
      }
      # maxSizes <- annualFires$size * 2
      dups <- duplicated(annualFires$cells)
      if (any(dups)) {
        annualFires <- annualFires[which(!dups), ] #
        maxSizes <- maxSizes[!dups]
        loci <- annualFires$cells[!dups]
      }
      ## spread() runs on this year's bounding box, not the whole landscape. It allocates
      ## landscape-length state on every call (SpaDES.tools spread.R: `integer(ncells)`), Nreps
      ## times per fire year, and the year's pixels fill a small part of the landscape (a median
      ## 10-33% by bounding box on ELFs 5.3.1, 5.3.2, 13.1). See cropToCells() for why the result
      ## is the same. Indices are mapped back below, so nothing after spread() sees the crop.
      crop <- cropToCells(r, c(shortAnnDT$pixelID, loci))
      spreadProbCrop <- numeric(crop$ncell)
      spreadProbCrop[crop$toCrop(shortAnnDT$pixelID)] <- shortAnnDT$spreadProb
      # if (any(cells[loci] == 0)) {
      spreadProbCrop[crop$toCrop(loci)] <- 1
      # }
      if (SpaDES.core::anyPlotting(plot.it)) {
        ## the plotting code further down reads the landscape-length vector
        cells[shortAnnDT$pixelID] <- shortAnnDT$spreadProb
        cells[loci] <- 1
      }
      ## No system.time() here. It defaulted to gcFirst = TRUE, so every fire year
      ## forced a full garbage collection, and `st` was assigned and never read.
      ## On ELF 12.4 (1.7M cells, 18 fire years, a ~20-40 GB process) that was
      ## 109.96 s of a 134 s objective-function evaluation -- 82% of self time, with
      ## the spread() calls it was timing accounting for 6.4 s. It also explained why
      ## dropping Nreps from 25 to 5 barely moved the total (82.7 s to 77.6 s): the
      ## forced collection is once per year, outside the replicate loop.
      ## spreadCpp() instead of spread(): same rules (generations, one draw per
      ## burning-cell/unburned-neighbour pair against the neighbour's probability,
      ## one fire per cell, a maxSize that is never exceeded, NA unburnable) but a
      ## plain C++ loop, 2-3x faster on this call shape. It is NOT bit-compatible
      ## with spread(), so cached fits from before this change are not comparable.
      ## The landscape is already cropped above, which matters more for spreadCpp()
      ## than it did for spread().
      ## Per-fire random effect (fitFireSpreadSD in .objfunSpreadFit): each fire's buffer pixels get
      ## that fire's eps on the logit scale, drawn anew in each replicate. A pixel in two buffers takes
      ## the first fire's. Still one spreadCpp() call per replicate.
      if (fireSpreadSD > 0) {
        fb <- unique(annualFireBufferedDT[, list(pixelID, ids)], by = "pixelID")
        fireIds <- unique(annualFires$ids)
        pixFire <- match(fb$ids[match(shortAnnDT$pixelID, fb$pixelID)], fireIds)
        pixCrop <- crop$toCrop(shortAnnDT$pixelID)
        lp <- stats::qlogis(shortAnnDT$spreadProb)
        lociCrop <- crop$toCrop(loci)
      }
      spreadState <- lapply(seq_len(Nreps), function(i) {
        sp <- spreadProbCrop
        if (fireSpreadSD > 0) {
          eps <- stats::rnorm(length(fireIds), 0, fireSpreadSD)
          e <- eps[pixFire]
          e[is.na(e)] <- 0
          sp[pixCrop] <- stats::plogis(lp + e)
          sp[lociCrop] <- 1
        }
        SpaDES.tools::spreadCpp(
          landscape = crop$r,
          loci = crop$toCrop(loci),
          spreadProb = sp,
          maxSize = maxSizes
        )
      })
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
      ## back to landscape cell indices
      set(spreadState, NULL, "initialLocus", crop$toFull(spreadState$initialLocus))
      set(spreadState, NULL, "indices", crop$toFull(spreadState$indices))
      if (isTRUE(returnSims)) {
        simSize <- spreadState[, list(sim = .N), by = c("rep", "initialLocus")]
        f <- data.table::as.data.table(annualFires)[, list(ids, initialLocus = cells, size)]
        return(list(sims = data.table::data.table(yr = as.character(yr), simSize[f, on = "initialLocus"]),
                    pSummary = pSummary))
      }
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
        ## Observed and simulated sizes on the SAME (square-root) scale. With `x = size[1]` the
        ## observed size sat beyond sqrt(maxSize) for any fire above ~12 pixels, so its likelihood
        ## was the `minLik` floor whatever the parameters.
        emp <- emp[N > 1, list(size = size[1],
                               lik = if(.N > 1) {
                                 if (identical(sizeLik, "t")) {
                                   sizeLikT(size[1], N, sizeLikDf)
                                 } else {
                                   EnvStats::demp(x = sqrt(size[1]), obs = sqrt(N))
                                 }
                               } else {
                                 0
                               }), by = "ids"]
        # emp <- emp[N > 1, list(size = size[1],
        #                            lik = ifelse(.N > 1,
        #                                         EnvStats::demp(x = size[1], obs = sqrt(N)),
        #                                         0)), by = "ids"]
        # emp <- emp[N > 1, list(size = size[1], lik = EnvStats::demp(x = size[1], obs = N)), by = "ids"]
        ## A weight MULTIPLIES the log-likelihood. This was log(lik * log(size)), i.e.
        ## log(lik) + log(log(size)): an offset, under which a fire's likelihood moved the SNLL by the
        ## same amount whatever its size.
        set(emp, NULL, "lik", sizeWeight(emp$size, weighted) / sizeWeightMean * log(pmax(minLik, emp$lik)))
        
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
      if (isTRUE(doADTest)) {
        ## The adTest compares distributions, and the observed sample is single fires, so it gets
        ## every replicate's fire, not `fireSizes`: a per-fire mean over Nreps has a far shorter tail.
        ret <- append(ret, list(allFireSizes = spreadState[, .N, by = c("rep", "initialLocus")][["N"]]))
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
      if (isTRUE(returnSims)) {
        f <- data.table::as.data.table(annualFires)[, list(ids, initialLocus = cells, size)]
        return(list(sims = data.table::data.table(yr = as.character(yr), rep = NA_integer_, f,
                                                  sim = NA_integer_), pSummary = pSummary))
      }
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

#' Automatic weight for the Anderson-Darling term
#'
#' The multiplier that gives the `"adTest"` term and the fire-size SNLL about equal influence on the
#' objective, `c * sqrt(nFires)`. Both terms are sums over the fitted fires, so both grow with how
#' many there are; measured across six ELFs (138 to 1624 fires, 2,708 to 1,114,569 pixels burned,
#' 14 members of each fit's final population), the spread of the SNLL across members grew with
#' `nFires^0.5` to `nFires^1`, and neither term scaled with the ELF's area or its area burned. `c`
#' comes from the same measurements, one per `sizeLik` and `weighted` combination. Under it the AD
#' term's share of the influence stayed within 0.36-0.68 of the total, against 0.16-0.94 under the
#' fixed 50.
#'
#' The constants were measured on populations fitted with an older objective, and a spread estimated
#' from 14 members is coarse, so treat the result as the right order of magnitude rather than exact.
#'
#' @param nFires Number of fires being fitted (after `minFireSize`).
#' @param sizeLik As in [.objfunSpreadFit()].
#' @param weighted As in [.objfunSpreadFit()]; `TRUE` uses the unweighted constants.
#'
#' @return A single numeric.
#' @export
adWeightAuto <- function(nFires, sizeLik = "kde", weighted = FALSE) {
  cc <- c(kde = 1.897, kde_sqrt = 5.098, t = 0.282, t_sqrt = 0.804)
  key <- paste0(match.arg(sizeLik, c("kde", "t")), if (identical(weighted, "sqrt")) "_sqrt" else "")
  unname(cc[[key]]) * sqrt(nFires)
}

#' Student-t likelihood of an observed fire size given simulated sizes
#'
#' The density, at `sqrt(size)`, of a Student-t centred on the mean of `sqrt(N)` with the standard
#' deviation of `sqrt(N)` as its scale. The scale has a lower bound of 0.5, because replicates that
#' all stop at the same size (e.g. at `maxSize`) have a standard deviation of 0.
#'
#' @param size Observed fire size, in pixels.
#' @param N Simulated sizes of that fire, in pixels; at least two.
#' @param df Degrees of freedom.
#'
#' @return A single numeric: the likelihood.
#' @keywords internal
#' @importFrom stats dt sd
sizeLikT <- function(size, N, df) {
  sc <- max(sd(sqrt(N)), 0.5)
  dt((sqrt(size) - mean(sqrt(N))) / sc, df) / sc
}

#' Size weights for the fire-size likelihood
#'
#' @param size Observed fire sizes, in pixels.
#' @param weighted `FALSE` (weights of 1), `TRUE` or `"log"` (`log(size)`), or `"sqrt"` (`sqrt(size)`).
#'
#' @return A numeric vector the length of `size`.
#' @keywords internal
sizeWeight <- function(size, weighted) {
  if (isFALSE(weighted)) return(rep(1, length(size)))
  if (isTRUE(weighted)) weighted <- "log"
  switch(match.arg(weighted, c("log", "sqrt")), log = log(size), sqrt = sqrt(size))
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
#' @param covCentre Optional named list or vector: `covCentre[[col]]` is subtracted from column `col`
#'   after rescaling, mutual exclusivity and the range assertion. `NULL` (default) does no centring.
#'   Supply the fitted values at prediction.
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
                                      lowerSpreadProb, covCentre = NULL) {
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
    assertCovariateRange(shortAnnDT, colsToUse)
    if (logisticPars[1] > maxFireSpread) {
      warning(
        "The first parameter of the logistic is > ", maxFireSpread, ".",
        "The parameter should be lowered."
      )
    }
  }

  ## Centre last: after the mutual-exclusivity step, whose 0 means "this pixel does not have this
  ## cover" only on the uncentred scale, and after the range assertion, which centred (negative)
  ## values would fail. Fitting and prediction both pass through here, so both centre identically.
  if (!is.null(covCentre)) {
    for (cn in intersect(names(covCentre), colsToUse)) {
      set(shortAnnDT, NULL, cn, shortAnnDT[[cn]] - covCentre[[cn]])
    }
  }

  # covPars <- tail(x = par, n = parsModel)
  # logisticPars <- head(x = par, n = length(par) - parsModel)
  shortAnnDT
}

#' Crop a raster's cell index to the bounding box of some cells
#'
#' `SpaDES.tools::spread()` allocates landscape-length state on every call, so its cost grows with
#' the landscape even when the fires touch a few thousand cells. This gives it a smaller landscape:
#' the bounding box of `cells`, plus a margin, as an empty raster with functions that map cell
#' indices into it and back.
#'
#' A spread on the crop is the same as on the full landscape, random draws included, when every
#' cell with a non-zero `spreadProb` is in `cells`: a fire can then only occupy those cells, the
#' one-cell margin keeps all eight neighbours of each inside the crop (so each step draws for the
#' same neighbours), and both mappings are monotonic, so neighbours keep their order. The margin
#' stops at the landscape's own edge, where the full landscape has no neighbour either.
#'
#' @param r A `SpatRaster`; only its geometry is used.
#' @param cells Integer vector of cell indices of `r` that must be inside the crop.
#' @param margin Number of cells to add on each side.
#'
#' @return A list: `r`, the empty cropped `SpatRaster`; `ncell`, its number of cells; `toCrop()`
#'   and `toFull()`, functions mapping cell indices of `r` to the crop's and back.
#' @keywords internal
cropToCells <- function(r, cells, margin = 1L) {
  nc <- as.integer(terra::ncol(r))
  nr <- as.integer(terra::nrow(r))
  cells <- as.integer(cells)
  rows <- (cells - 1L) %/% nc # 0-based
  cols <- (cells - 1L) %% nc
  r0 <- max(min(rows) - margin, 0L)
  r1 <- min(max(rows) + margin, nr - 1L)
  c0 <- max(min(cols) - margin, 0L)
  c1 <- min(max(cols) + margin, nc - 1L)
  nrCrop <- r1 - r0 + 1L
  ncCrop <- c1 - c0 + 1L
  res <- terra::res(r)
  rCrop <- terra::rast(
    nrows = nrCrop, ncols = ncCrop, crs = terra::crs(r),
    xmin = terra::xmin(r) + c0 * res[1], xmax = terra::xmin(r) + (c1 + 1L) * res[1],
    ymin = terra::ymax(r) - (r1 + 1L) * res[2], ymax = terra::ymax(r) - r0 * res[2]
  )
  list(
    r = rCrop, ncell = nrCrop * ncCrop,
    toCrop = function(x) ((x - 1L) %/% nc - r0) * ncCrop + ((x - 1L) %% nc - c0) + 1L,
    toFull = function(x) ((x - 1L) %/% ncCrop + r0) * nc + ((x - 1L) %% ncCrop + c0) + 1L
  )
}
