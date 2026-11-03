yearSpreadSDTxt <- "yearSpreadSD"

utils::globalVariables(c(
  "..colsToKeep", "..colsToUse", ".N", "buffer", "burned", "burnedClass", "cut",
  "id", "ids", "indices", "initialLocus", "N", "numAvailPixels", "open", "pixelID", "prob",
  "simFireSize", "size", "spreadProb", "value"
))

#' Objective function for `fireSense_spreadFit` module
#'
#' @param par parameters: the logistic parameters (`maxAsymptote` first, then `upperTail1` if the link has
#'   one; without `hillSlope1` and `inflectionPoint1`, which are fixed at 1 -- see `fixLogisticPars()`), then the covariate coefficients (the intercept's first, when
#'   the formula has one -- see [spreadDesignCols()]), and optionally
#'   `yearSpreadSD` last (see `fitYearSpreadSD`).
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
#' @param maxFireSpread A value for `spreadProb` that is considered impossible to go above. Default
#'   [spreadProbCeiling], also the upper bound of `maxAsymptote`. `maxAsymptote` is the spread-probability
#'   ceiling in a typical year. The year random effect (`yearSpreadSD`) is added on the logit of the final
#'   spread probability, after the link, so in a given year p can exceed `maxAsymptote` or fall below
#'   `lowerSpreadProb`; that is intended. There is no absolute cap, because `spreadCpp` does not need one;
#'   `maxAsymptote` is bounded because runaway fires are slow to simulate and wasted if the parameters are wrong.
#'
#' @param minFireSize Integer. Minimum fire size (in pixels) to include when
#'   scoring simulated against historical fires; fires smaller than this are
#'   filtered out of the comparison. Default `1`.
#'
#' @param escapeSizeHa Numeric or `NULL`. Size (ha) a fire must reach to count as escaped. `NULL` (the
#'   default) keeps the historical behaviour. When set, only observed fires of at least
#'   [escapeSizePixels()] pixels are fitted, and every simulated fire starts as an escaped fire: it burns
#'   its first [escapeSizePixels()] cells regardless of `spreadProb` (`SpaDES.tools::spreadCpp(minSize =)`),
#'   then spreads with the fitted probabilities. So the spread model is fitted as "given the fire escaped".
#'
#' @param jumpTries,jumpMeanDist Used with `escapeSizeHa`: passed to `SpaDES.tools::spreadCpp()`, so a fire
#'   stuck below the escape size (no burnable neighbour left) can jump to nearby burnable land. `jumpTries`
#'   is how many times it may try, `jumpMeanDist` the mean jump distance in pixels. `0` (default) is off.
#'
#' @param yearAreaWeight Weight of the annual-area term: each fit year's observed area burned scored against
#'   that year's simulated totals, one per replicate, with the per-fire size likelihood (`yearAreaNLL()`).
#'   `0` (default) leaves it out; `"auto"` is (number of fitted fires) / (number of fit years), so the year
#'   view and the per-fire view carry equal total weight.
#'
#' @param areaDistWeight Weight of the area-weighted size-distribution term (`areaWeightedCvM()`), which
#'   compares simulated and observed fires by the share of area burned they make up. `0` (default) leaves
#'   it out; `"auto"` uses the AD term's weight, `adWeightAuto()`.
#'
#' @param returnTerms Logical. `TRUE` returns the objective with each term separately (diagnostics).
#'
#' @template mutuallyExclusive
#'
#' @param doAssertions Logical. If `TRUE`, the default, the function will test a few minor things
#'   for consistency. This should be set to `FALSE` for operational situations, as the assertions
#'   take some small amount of time.
#'
#' @param tests One or more of `"mad"`, `"adTest"`, `"pitAD"` or `"crps"`. Default: `"crps"`. (`"snll_fs"` and `"SNLL"` were removed; `"snll_fs"` stops with an error.) `"crps"` is a different objective made of proper scores, each smooth in the parameters: the per-fire CRPS on log size ([crpsFireSizes()], summed over fitted fires), plus the pooled `adTest`, plus `areaDist` ([areaWeightedCvM()]), plus the CRPS of the log total area burned over all fit years ([crpsTotalArea()]), plus the CRPS of the 10 largest fires ([crpsTopFires()]), plus, if its weight is not 0, the mean over fit years of the CRPS of each year's log area burned ([crpsAnnualArea()]), each times its weight in `crpsWeights`. It takes over the `adTest` and `areaDist` terms (their weights come from `crpsWeights`, not `adWeight` or `areaDistWeight`), has its own annual-area term, off by default (`crpsWeights$annualArea`; `yearAreaWeight` is not used with it). A cut replicate (burned an open buffer edge pixel) is scored in the per-fire term as at least the observed size, and in the other terms at the size it burned. Combined with the other tests it adds to them. `"pitAD"` is a rank-based Anderson-Darling term that does not depend on the fit buffers: each fitted fire's observed size ranked among its own replicates ([pitMidRank()]), scored against Uniform(0, 1) ([adUniform()]), weighted like `"adTest"` (`adWeight`).
#'
#' @param Nreps Integer. The number of replicates, per ignition, to run.
#'
#' @param crpsWeights Named list of the multipliers of the `"crps"` test's terms, `fire`, `adTest`, `areaDist`,
#'   `totalArea`, `top10` and `annualArea` (0: left out); names left out take the values of [crpsWeightsDefault()] (set by the 2026-10-28
#'   profile on ELFs 6.6.2, 14.3 and 11.3, "W3": 1, 0.34, 2.6, 164 and 16.5). `NULL` is all defaults. Not used
#'   without `"crps"` in `tests`.
#'
#' @param gateRefusal What a year that fails [spreadProbGateTest()] gets. `"minLik"` (default): it is not
#'   simulated and scores as if every fire stayed one pixel (`"crps"`), or the `minLik` floor (the
#'   annual-area term). `"simulate"`:
#'   the gates are not applied, the year is simulated and scored like any other, so a parameter set is never
#'   refused for its spread probabilities alone.
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
#' @param covCentre Optional named list or vector of values subtracted from the rescaled covariates
#'   of the same name. `NULL` (default) does no centring. Must be the same at fit and predict.
#'   [spreadCovCentre()] gives the centres to use with an intercept.
#' @param lanscape1stQuantileThresh A `spreadProb` value that represents a threshold for the
#'   1st quantile of the `spreadProbs` on the landscape; if that quantile is above this
#'   number, then the `.objFunSpredFit` will bail because it is "too burny" a landscape.
#'   Default 0.265, meaning if only 25%% of the pixels on the landscape are below
#'   this `spreadProb`, then it will bail.
#'
#' @param weighted Only picks the constant of `adWeight = "auto"` ([adWeightAuto()]): `"sqrt"` uses the
#'   one measured for size-weighted fits, anything else the unweighted one.
#'
#' @param sizeLik How the annual-area term (`yearAreaWeight`) takes a year's likelihood from its `Nreps`
#'   simulated totals, on the square-root scale. `"kde"` (default): their kernel density at the observed
#'   total ([EnvStats::demp()]); it is zero away from the simulated totals, so a year they never reach
#'   scores the `minLik` floor however far off it is. `"t"`: a Student-t centred on their mean with
#'   their standard deviation as scale; no floor, and a miss costs more the further it is. It also
#'   picks the constant of `adWeight = "auto"`.
#'
#' @param sizeLikDf Degrees of freedom of the Student-t when `sizeLik = "t"`. Default 5. Lower
#'   values forgive a missed fire more.
#'
#' @param adWeight What the `"adTest"` statistic is multiplied by. `"auto"` (default) uses
#'   [adWeightAuto()]. A number is used as given. Not used with `"crps"`, which has its own weight.
#'
#' @param link Passed to [logisticAll()]: `NULL` for the form implied by the logistic parameters, or
#'   `"logistic3pUpper"` for the 3-parameter logistic with Stukel's upper tail, whose fourth logistic
#'   parameter is `upperTail1`, or `"powerLaw"` for [powerLawLink()]: `maxAsymptote` is then `p_c`, fitted, and
#'   the year effect (`fitYearSpreadSD`) is added to the linear predictor `x`, not to the logit of `p`, so a year
#'   cannot take `p` past `p_c`. Pass it explicitly when fitting: `par` may arrive unnamed.
#' @param powerLawGamma The `"powerLaw"` link's `gamma` ([powerLawLink()]), fixed, not fitted; used only with
#'   `link = "powerLaw"`.
#' @param powerLawShape The `"powerLaw"` link's curvature exponent `shape` ([powerLawLink()]), fixed, not fitted;
#'   used only with `link = "powerLaw"`. Default 1.
#' @param climateRef `NULL` (default): the climate covariates are inside the logistic, under the
#'   `maxAsymptote` ceiling, as always. A named numeric or list ([spreadClimateRef()]): climate is taken
#'   outside the ceiling, `p = plogis(qlogis(p0) + beta_clim * (climate - climateRef))` with `p0` the
#'   logistic of the other covariates; see [spreadProbFromCovs()], which computes both links. The year
#'   effect (`fitYearSpreadSD`) is added on the same logit, after climate. The gates
#'   ([spreadProbGateTest()]) and the saturation summary are then applied to `p0`, the spread probability
#'   at the reference climate, not to `p`.
#'
#' @param returnSims If `TRUE`, return the simulated fires instead of the objective value: a
#'   `data.table` with one row per fire year x fire x replicate -- `yr`, `rep`, `initialLocus`, `ids`,
#'   the simulated size `sim` and the observed `size`, in pixels. A year the objective declines to
#'   simulate (its spread probabilities fail the sanity checks) keeps its fires with `rep` and `sim`
#'   `NA`, so the caller can count them. Every year is simulated. For diagnostics and validation; the value is not an objective.
#'   The result carries `attr(, "spreadProb")`: the spread probabilities of every fitted pixel-year,
#'   summarised by [spreadProbSummary()], for [linkSaturation()].
#' @param fitYearSpreadSD If `TRUE`, the last element of `par` is the sd of a per-year random effect on
#'   the logit of spread probability, a seasonal departure: in each replicate the year draws one
#'   `eps ~ N(0, sd^2)` and every pixel of that year spreads with `plogis(qlogis(p) + eps)`, so all of the
#'   year's fires are bigger or smaller together. Without it the simulated fire sizes were too alike:
#'   medians too large and the largest fires too small at once (2026-09-23). `sd = 0` is the model
#'   without it, exactly (no random numbers are drawn). `NULL` (default): `TRUE` when `par` is named and
#'   includes `yearSpreadSD`. DEoptim passes `par` unnamed, so `runDEoptim()` sets it explicitly.
#' @param returnBurned If `TRUE`, also record which pixels burned: implies `returnSims = TRUE`, and the
#'   result carries `attr(, "burned")`, a named list with one element per simulated fire year, each a
#'   `data.table` with the replicate `rep` and the `pixelID` of every pixel that burned in it (ignition
#'   cells included). A year the objective declines to simulate is `NULL`. For
#'   [spreadFitValidationData()]. `FALSE` (default) records nothing and changes nothing.
#' @param penaliseRunaways A simulated fire that burns an open edge pixel of any of its year's buffers
#'   (a buffer pixel next to flammable land outside every buffer, see [bufferOpen()]) is cut: it would
#'   have kept burning, so it is "at least this big". The per-fire `"crps"` term scores it as at least the
#'   observed size; the other size-based terms (`"adTest"`, `yearAreaWeight`, `areaDistWeight`, `"mad"`)
#'   score the size it burned. That is always so. `TRUE`
#'   (default) also stops such a fire at the first open pixel of its own buffer
#'   (`SpaDES.tools::spreadCpp(stopCells =)`), as nothing it burns after that is scored; `FALSE` lets it
#'   burn on. Not with `returnSims = TRUE`: the simulated sizes are returned as simulated. Fires are not
#'   capped at any size: spread is bounded only by the year's buffers. The open pixels are found once per fit
#'   (`runDEoptim()` adds `edge` and `open` columns to `fireBufferedListDT`, [addBufferEdge()]), or here,
#'   for a table without them.
#' @param runawayEdgeFrac,runawayEdgeMin Ignored. A replicate was a runaway when it burned at least
#'   `max(runawayEdgeMin, ceiling(runawayEdgeFrac * n))` of the `n` edge pixels of its own buffer, open or
#'   not; it is now cut at the first open pixel of any buffer of its year.
#' @param runawaySize,runawayBufferMultiple Ignored, with a warning unless `NULL`: they set the size a
#'   runaway counted as in the size-distribution and area terms. A cut replicate is scored at the size it
#'   reached.
#'
#' @param verbose If >= 2, then this will show more information about `spreadProb` fitting.
#'
#' @param lowerSpreadProb Numeric. Lower bound for `spreadProb`; if a candidate
#'   `spreadProb` falls at or below this, the call exits early to avoid wasted
#'   simulation work. Default [spreadProbFloor].
#'
#' @param mutuallyExclusive Named list of vectors describing groups of model
#'   terms that must not be active together (mutually exclusive in the
#'   formula). Default `list("youngAge" = c("class", "nf"))`.
#'
#' @param ... This is not used here, but allows for extraneous arguments to not break this function.
#'   The removed arguments `capSizes`, `penaliseCapHits` and `thresh` arrive here: they are ignored with a
#'   deprecation warning (fires are no longer capped at a size; use `penaliseRunaways`).
#'
#' @return
#' The objective, to be minimised: the sum of the terms `tests` asks for. With `returnTerms = TRUE`, a
#' named vector of the objective and each term.
#'
#' @export
#' @importFrom data.table := rbindlist set setDT setDTthreads setnames setorderv
#' @importFrom EnvStats demp
#' @importFrom graphics hist mtext
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
                             maxFireSpread = spreadProbCeiling,
                             lowerSpreadProb = spreadProbFloor,
                             minFireSize = 2,
                             escapeSizeHa = NULL,
                             jumpTries = 0,
                             jumpMeanDist = 0,
                             yearAreaWeight = 0,
                             areaDistWeight = 0,
                             returnTerms = FALSE,
                             tests = "crps",
                             Nreps = 10,
                             mutuallyExclusive = list("youngAge" = c("class", "nf")),
                             doAssertions = TRUE,
                             plot.it = FALSE, # TODO parameterize this line for plotting
                             objFunCoresInternal = 1,
                             lanscape1stQuantileThresh = 0.265,
                             covCentre = NULL,
                             weighted = TRUE,
                             sizeLik = "kde",
                             sizeLikDf = 5,
                             adWeight = "auto",
                             link = NULL,
                             climateRef = NULL,
                             powerLawGamma = spreadPowerLawGamma,
                             powerLawShape = spreadPowerLawShape,
                             returnSims = FALSE,
                             returnBurned = FALSE,
                             penaliseRunaways = TRUE,
                             runawayEdgeFrac = fireSenseRunawayEdgeFrac,
                             runawayEdgeMin = fireSenseRunawayEdgeMin,
                             runawaySize = NULL,
                             runawayBufferMultiple = NULL,
                             fitYearSpreadSD = NULL,
                             # bufferedRealHistoricalFiresList,
                             verbose = 2,
                             crpsWeights = NULL,
                             gateRefusal = c("minLik", "simulate"),
                             ...) { # fireSense_SpreadFitRaster
  # Optimization's objective function
  # lapply(historicalFires, setDT)

  data.table::setDTthreads(1)
  deprecatedArgs(c(...names(), names(which(c(
    runawaySize = !is.null(runawaySize), runawayBufferMultiple = !is.null(runawayBufferMultiple))))))
  checkObjectiveTests(tests)
  sizeLik <- match.arg(sizeLik, c("kde", "t"))
  gateRefusal <- match.arg(gateRefusal)
  if (isTRUE(returnBurned)) returnSims <- TRUE
  ## the per-year random effect's sd, fitted as the LAST element of `par`
  ## and hillSlope1 and inflectionPoint1 are fixed at 1, not fitted -- see fixLogisticPars() for why.
  parSplit <- splitSpreadPar(par, fitYearSpreadSD)
  par <- parSplit$par
  yearSpreadSD <- parSplit$yearSpreadSD

  doMADTest <- any(grepl("mad", tolower(tests)))
  doCrps <- any(tolower(tests) == "crps")
  crpsWeights <- utils::modifyList(crpsWeightsDefault(), as.list(crpsWeights))
  ## "crps" takes over these two terms, with its own weights
  doADTest <- any(grepl("adtest", tolower(tests))) || doCrps
  doPitADTest <- any(tolower(tests) == "pitad")
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
  ## the buffers' edge and open pixels: made once by runDEoptim(); here for a caller that did not
  fireBufferedListDT <- addBufferEdge(fireBufferedListDT, landscape)
  # lapply(annualDTx1000, setDT)
  lapply(nonAnnualDTx1000, setDT)
  # lapply(fireBufferedListDT, setDT)
  # dtThreadsOrig <- data.table::setDTthreads(1)
  if (is(formulaToFit, "formula")) {
    stop("formulaToFit must be provided as a charater string because it takes too much RAM otherwise.")
  }
  formulaToFit <- as.formula(formulaToFit, env = .GlobalEnv)
  ## the terms, preceded by the intercept when the formula has one (see spreadDesignCols())
  colsToUse <- spreadDesignCols(formulaToFit)
  # How many of the parameters belong to the model?
  parsModel <- length(colsToUse)
  ncells <- ncell(landscape)

  r <- rast(landscape)
  ## an escaped fire has at least `escapeMinPx` pixels: only such observed fires are fitted, and every
  ## simulated fire burns that many cells first (spreadCpp(minSize =)). NULL keeps the historical rules.
  escapeMinPx <- if (!is.null(escapeSizeHa)) escapeSizePixels(escapeSizeHa, landscape)
  minFireSize <- effectiveMinFireSize(minFireSize, escapeMinPx)
  if (!is.null(escapeMinPx) && isTRUE(jumpTries > 0) &&
      !all(c("jumpTries", "jumpMeanDist") %in% names(formals(SpaDES.tools::spreadCpp))))
    stop("jumpTries needs a SpaDES.tools whose spreadCpp() has `jumpTries` and `jumpMeanDist`")
  if (doCrps) { # its own terms and weights: its annual-area term is `crpsWeights$annualArea`
    yearAreaWeight <- 0
    areaDistWeight <- crpsWeights$areaDist
    adWeight <- crpsWeights$adTest
  }
  doYearArea <- !identical(yearAreaWeight, 0) && !identical(yearAreaWeight, 0L)
  doAreaDist <- !identical(areaDistWeight, 0) && !identical(areaDistWeight, 0L)
  years <- as.character(names(annualDTx1000))
  names(years) <- years
  ## numeric, not integer: it receives spreadProb (double), and assigning doubles into an
  ## integer vector coerces the whole landscape-length vector on every fire year.
  cells <- numeric(ncells)
  # Nreps <- 10
  indexNonAnnual <- nonAnnualIndex(nonAnnualDTx1000)
  historicalFiresAboveMin <- firesAboveMinSize(historicalFires, minFireSize)

  ## the first block is the two years with the most burned area, the second block is the rest
  largest <- firstBlockYears(historicalFiresAboveMin)
  smallest <- setdiff(names(historicalFiresAboveMin), largest)
  lrgSmallFireYears <- list(large = largest, small = smallest)
  objFunRes <- 0
  objFunResList <- list() # "mad": each batch's value and number of fires
  fireSizesList <- list() # simulated fire sizes, pooled across batches for the adTest
  yrsDoneList <- list() # the years contributing to fireSizesList, kept in lockstep
  pitList <- list() # each fitted fire's PIT (pitMidRank()), pooled across batches for the "pitAD" term
  crpsFireList <- list() # "crps": each batch's per-fire CRPS, summed over its fires
  crpsSizesList <- list() #   ... and every simulated fire of its years (year, rep, N)
  cutList <- list() #   ... and how many of its simulated fires were cut at an open buffer edge, of how many
  yearAreaList <- list() # yearAreaWeight: each year's simulated annual totals, one per replicate
  simsList <- list() # returnSims: the simulated fires of each batch of years
  pSumList <- list() #   ... and the spread probabilities of their pixels
  burnedList <- list() # returnBurned: per year, the pixels burned in each replicate
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
        r = r,
        doMADTest = doMADTest, doADTest = doADTest || doAreaDist, doPitADTest = doPitADTest,
        cells = cells,
        covCentre = covCentre,
        link = link, climateRef = climateRef,
        powerLawGamma = powerLawGamma, powerLawShape = powerLawShape,
        returnSims = returnSims, returnBurned = returnBurned, stopCutFires = isTRUE(penaliseRunaways),
        yearSpreadSD = yearSpreadSD,
        escapeMinPx = escapeMinPx, jumpTries = jumpTries, jumpMeanDist = jumpMeanDist,
        doYearArea = doYearArea, doCrps = doCrps, gateRefusal = gateRefusal,
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
        if (isTRUE(returnBurned)) burnedList[yrs] <- results$burned
        next
      }

      if (isTRUE(doADTest) || isTRUE(doAreaDist)) {
        fireSizesList[[ii]] <- unlist(results$allFireSizes)
        yrsDoneList[[ii]] <- yrs
      }
      if (isTRUE(doPitADTest)) {
        pitList[[ii]] <- unlist(results$pit)
      }
      if (isTRUE(doCrps)) {
        crpsFireList[[ii]] <- sum(unlist(results$crpsFire))
        crpsSizesList[[ii]] <- data.table::rbindlist(results$crpsSizes)
        cutList[[ii]] <- c(cut = sum(unlist(results$runaways)), sims = sum(unlist(results$nSims)))
      }
      if (isTRUE(doYearArea)) {
        yearAreaList[yrs] <- results$annualAreaByRep # a year the objective refused to simulate stays absent
      }

      if (isTRUE(doMADTest)) {
        a <- purrr::map2(historicalFiresAboveMin[yrs], results$fireSizes, function(x, y) {
          data.table(x, simFireSize = y)
        })

        a <- rbindlist(a)
        a[, dev := abs(size - simFireSize)]
        # a[, devLog := sqrt(abs(size - simFireSize))]
        mad <- round(mean(a$dev), 1)
        # mad <- round(mean(a$devLog), 1)
        objFunResList[ii] <- list(list(objFunRes = mad, nFires = NROW(a)))
      }
    }
  } # run through 2nd batch of smaller fires

  if (isTRUE(returnSims)) {
    out <- data.table::rbindlist(simsList, use.names = TRUE)
    data.table::setattr(out, "spreadProb", combineSpreadProbSummaries(unlist(pSumList, recursive = FALSE)))
    if (isTRUE(returnBurned)) data.table::setattr(out, "burned", burnedList)
    return(out)
  }
  bb <- purrr::transpose(objFunResList)
  bb <- purrr::map(bb, unlist)

  if (isTRUE(doMADTest)) {
    totalNFires <- sum(bb$nFires)
    objFunRes <- do.call(sum, purrr::map2(
      bb$objFunRes, bb$nFires,
      function(.x, .y) .x * .y
    )) / totalNFires
  }
  if (isTRUE(doCrps)) {
    ## Per-fire, total-area and largest-fire CRPS. All fit years are scored (a refused year counts as every
    ## fire staying one pixel). The total and the largest fires are pooled over the
    ## fit years, so they are found here, not per batch.
    fitYrs <- unlist(lrgSmallFireYears)
    obsSizes <- unlist(lapply(historicalFiresAboveMin[fitYrs], function(x) x$size))
    crpsSizes <- data.table::rbindlist(crpsSizesList)
    # the replicates of different years are one replicate: a rep's years add (total) or compete (largest)
    crpsFire <- sum(unlist(crpsFireList))
    crpsTotal <- crpsTotalArea(crpsSizes, obsSizes, Nreps)
    crpsTop <- crpsTopFires(crpsSizes, obsSizes, Nreps)
    objFunRes <- objFunRes + crpsWeights$fire * crpsFire + crpsWeights$totalArea * crpsTotal +
      crpsWeights$top10 * crpsTop
    crpsAnn <- NA_real_
    if (!identical(as.numeric(crpsWeights$annualArea), 0)) { # each year's observed area, from the same fires
      obsYear <- vapply(historicalFiresAboveMin[fitYrs], function(x) sum(x$size), numeric(1))
      crpsAnn <- crpsAnnualArea(crpsSizes, obsYear, Nreps)
      objFunRes <- objFunRes + crpsWeights$annualArea * crpsAnn
    }
    if (verbose > 1) {
      cuts <- do.call(rbind, cutList)
      cutShare <- sum(cuts[, "cut"]) / max(1, sum(cuts[, "sims"]))
      print(paste0("  ", Sys.getpid(), " crpsFire:", round(crpsFire, 2), "; crpsTotalArea:", round(crpsTotal, 4),
                   "; crpsTop10:", round(crpsTop, 3), "; crpsAnnualArea:", round(crpsAnn, 4),
                   "; cut share: ", round(cutShare, 3), "; "))
    }
  }
  ## The Anderson-Darling test compares whole distributions, so it runs once, after
  ## all batches, on the simulated fire sizes pooled across batches -- and against the
  ## observed sizes from *those same years*. Previously it ran inside the loop under
  ## `ii == 2`, which pooled only the final batch's simulated sizes while comparing
  ## them to every year's observed sizes; the omitted years are the 2 with the largest
  ## area burned, so the observed sample kept an upper tail the simulated sample could
  ## not have.
  if (isTRUE(doADTest)) {
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
  ## The rank-based AD term: each fitted fire's observed size ranked among its own replicates (pitMidRank()), and
  ## those ranks scored against Uniform(0, 1) (adUniform()). A replicate cut at an open buffer edge after it passed
  ## the observed size is larger whatever the buffer, so unlike the pooled adTest it does not depend on the buffer.
  pitPart <- NA_real_
  if (isTRUE(doPitADTest)) {
    u <- unlist(pitList)
    pitPart <- if (length(u)) adUniform(u, Nreps) else 1e6L
    pitPart <- pitPart * if (identical(adWeight, "auto")) {
      adWeightAuto(sum(vapply(historicalFiresAboveMin, nrow, integer(1))), sizeLik, weighted)
    } else {
      adWeight
    }
    objFunRes <- objFunRes + pitPart
    if (verbose > 1) {
      print(paste0("  ", Sys.getpid(), " pitAD:", pitPart, "; "))
    }
  }
  adPart <- if (isTRUE(doADTest)) adTest else NA_real_
  nFitFires <- sum(vapply(historicalFiresAboveMin, nrow, integer(1)))
  yearAreaTerm <- areaDistTerm <- NA_real_
  if (isTRUE(doYearArea)) {
    ## every fit year is scored; one with no simulated totals (refused as too burny / not spread out) gets minLik
    fitYrs <- unlist(lrgSmallFireYears)
    obsYear <- vapply(historicalFiresAboveMin[fitYrs], function(x) sum(x$size), numeric(1))
    yearAreaTerm <- yearAreaNLL(obsYear, yearAreaList, sizeLik = sizeLik, sizeLikDf = sizeLikDf) *
      if (identical(yearAreaWeight, "auto")) nFitFires / length(fitYrs) else yearAreaWeight
    objFunRes <- objFunRes + yearAreaTerm
  }
  if (isTRUE(doAreaDist)) {
    pooled <- pooledFireSizes(fireSizesList, yrsDoneList, historicalFiresAboveMin)
    areaDistTerm <- areaWeightedCvM(pooled$simulated, pooled$observed) *
      if (identical(areaDistWeight, "auto")) adWeightAuto(nFitFires, sizeLik, weighted) else areaDistWeight
    objFunRes <- objFunRes + areaDistTerm
  }
  if (verbose > 1 && (isTRUE(doYearArea) || isTRUE(doAreaDist)))
    print(paste0("  ", Sys.getpid(), " yearArea:", round(yearAreaTerm, 1), "; areaDist:", round(areaDistTerm, 1), "; "))
  if (length(objFunResList) > 1) {
    print(paste0(Sys.getpid(), "; FINISHED! ", Sys.time(), "; Objective Final: ", round(objFunRes, 0)))
  }
  if (isTRUE(returnTerms))
    return(c(objective = objFunRes, adTest = adPart, pitAD = pitPart, yearArea = yearAreaTerm,
             areaDist = areaDistTerm, crpsFire = if (isTRUE(doCrps)) crpsFire else NA_real_,
             crpsTotalArea = if (isTRUE(doCrps)) crpsTotal else NA_real_,
             crpsTop10 = if (isTRUE(doCrps)) crpsTop else NA_real_, nFires = nFitFires))
  return(objFunRes)
}

## The "SNLL_FS" and "SNLL" tests are gone. Called by the objective, and by runDEoptim() so a fit stops on the
## master, before it builds a cluster, not in every worker.
checkObjectiveTests <- function(tests) {
  if (any(grepl("snll", tolower(tests))))
    stop("The \"SNLL_FS\" and \"SNLL\" tests were removed; use tests = \"crps\".", call. = FALSE)
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

#' Mid-rank probability integral transform of one fire's observed size among its replicates
#'
#' The share of the fire's `Nreps` replicates smaller than its observed size, ties counted half. A replicate
#' with no row (it burned nothing) is smaller. A replicate cut at an open buffer edge has burned at least its
#' reached size: if that is at least the observed size it is larger (exact, whatever the buffer); if it is
#' smaller its order is unknown and it counts half.
#'
#' @param size Observed size of the fire (pixels).
#' @param sim Simulated sizes of its replicates that have a row (pixels).
#' @param cut Logical, as `sim`: was the replicate cut at an open buffer edge.
#' @param Nreps Number of replicates.
#' @return A number in `[0, 1]`.
#' @keywords internal
pitMidRank <- function(size, sim, cut, Nreps = length(sim)) {
  smaller <- !cut & sim < size
  tie <- !cut & sim == size
  unknown <- cut & sim < size
  (Nreps - length(sim) + sum(smaller) + 0.5 * (sum(tie) + sum(unknown))) / Nreps
}

#' Anderson-Darling statistic of PITs against Uniform(0, 1)
#'
#' @param u PITs ([pitMidRank()]), one per fire.
#' @param Nreps Number of replicates: `u` is kept `1 / (2 * Nreps)` from 0 and 1, the resolution of a
#'   mid-rank, so the statistic stays finite.
#' @return A single non-negative number; about 1 on average when `u` is uniform.
#' @keywords internal
adUniform <- function(u, Nreps) {
  eps <- 1 / (2 * Nreps)
  u <- sort(pmin(pmax(u, eps), 1 - eps))
  n <- length(u)
  i <- seq_len(n)
  -n - sum((2 * i - 1) * (log(u) + log(1 - rev(u)))) / n
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
                        r, Nreps, doMADTest, doADTest, doPitADTest = FALSE,
                        plot.it, verbose = 2, covCentre = NULL,
                        link = NULL, returnSims = FALSE, returnBurned = FALSE,
                        stopCutFires = FALSE, yearSpreadSD = 0, escapeMinPx = NULL, jumpTries = 0, jumpMeanDist = 0,
                        doYearArea = FALSE, climateRef = NULL, powerLawGamma = spreadPowerLawGamma,
                        powerLawShape = spreadPowerLawShape, doCrps = FALSE, gateRefusal = "minLik") {
  if (isTRUE(plot.it)) plot.it <- "screen"

  # needed because data.table objects were recovered from disk
  # Rescale to numerics and /1000
  # setDT(nonAnnDTx1000)
  # matrix multiplication
  parsList <- paramsSeparate(par, parsModel)
  logisticPars <- parsList[["logisticPars"]]
  covPars <- parsList[["covPars"]]

  shortAnnDT <- spreadProbCovsOfYear(
    yr, annDTx1000, nonAnnualDTx1000, indexNonAnnual, covMinMax, mutuallyExclusive, colsToUse,
    doAssertions, logisticPars, covPars, maxFireSpread, lowerSpreadProb, covCentre = covCentre
  )
  mat <- as.matrix(shortAnnDT[, ..colsToUse])
  set(shortAnnDT, NULL, "spreadProb", spreadProbFromCovs(mat, logisticPars, covPars, lowerSpreadProb, link,
                                                         climateRef = climateRef, powerLawGamma = powerLawGamma,
                                                         powerLawShape = powerLawShape))
  ## What the gates and the saturation summary look at. With climate outside the ceiling (climateRef) it is
  ## p0, the spread probability at the reference climate: the part the ceiling applies to. A dry year
  ## legitimately takes p past maxFireSpread under that link, and a gate on p would refuse exactly the
  ## years that carry the climate signal. Without climateRef it is spreadProb itself, as always.
  pGate <- if (is.null(climateRef)) shortAnnDT$spreadProb else
    as.numeric(spreadProbFromCovs(mat, logisticPars, covPars, lowerSpreadProb, link, climateRef = climateRef,
                                  climateAtRef = TRUE))
  pSummary <- if (isTRUE(returnSims)) spreadProbSummary(pGate, ceiling = logisticPars[1])
  ## Initialised here, not inside the branch below, because the function returns
  ## it unconditionally. With no test selected -- which is how the module's
  ## `debug` mode calls this, passing tests = "" -- doFitting is FALSE, the branch
  ## is skipped, and `return(ret)` used to fail with "object 'ret' not found".
  ## Nothing appends to it before the branch, so this is the same list it always
  ## was; the only change is that it exists when no test asked for anything.
  ret <- list()
  doFitting <- any(c(doMADTest, doADTest, doPitADTest, doCrps)) || isTRUE(returnSims)
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
    gate <- spreadProbGateTest(pGate, logisticPars[1], lowerSpreadProb, maxFireSpread,
                               lanscape1stQuantileThresh)
    nonEdgeValues <- gate$nonEdgeValues
    summ <- gate$summ
    medSPRight <- gate$medianOK
    spreadOutEnough <- !gate$notSpread
    lowSPLowEnough <- !gate$burny
    loci <- annualFires$cells

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
    
    if ((medSPRight && spreadOutEnough && lowSPLowEnough) || identical(gateRefusal, "simulate")) {
      if (verbose > 1) {
        print(paste0(
          " ",
          Sys.getpid(), ": ", yr, ", spreadProbs: ",
          paste(names(summ), round(summ, 3), collapse = ", ")
        ))
      }
      setDT(annualFireBufferedDT)
      ## the open edge pixels of the year's buffers (bufferOpen()): made by .objfunSpreadFit() or
      ## runDEoptim(); here, for a caller of this function alone, with every cell of `r` flammable
      if (!"open" %in% names(annualFireBufferedDT))
        annualFireBufferedDT <- addBufferEdge(list(annualFireBufferedDT), r)[[1]]
      openEdge <- annualFireBufferedDT[which(open), list(ids, pixelID)]
      ## the size (pixels) of each fire's buffer, burned and unburned
      bufferSizes <- annualFireBufferedDT[, list(numAvailPixels = .N), by = "ids"]
      if (doAssertions || SpaDES.core::anyPlotting(plot.it)) {
        tableOfBufferedMaps <- bufferSizes[annualFires, on = "ids"]
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
      dups <- duplicated(annualFires$cells)
      if (any(dups)) {
        annualFires <- annualFires[which(!dups), ] #
        loci <- annualFires$cells[!dups]
      }
      ## spread() runs on this year's bounding box, not the whole landscape. It allocates
      ## landscape-length state on every call (SpaDES.tools spread.R: `integer(ncells)`), Nreps
      ## times per fire year, and the year's pixels fill a small part of the landscape (a median
      ## 10-33% by bounding box on ELFs 5.3.1, 5.3.2, 13.1). See cropToCells() for why the result
      ## is the same. Indices are mapped back below, so nothing after spread() sees the crop.
      crop <- cropToCells(r, c(shortAnnDT$pixelID, loci))
      ## A cut replicate is scored the same however far it burns past the open edge, so with
      ## `stopCutFires` a fire stops at the first open pixel of its own buffer. (A fire that runs into a
      ## neighbour's buffer and through its open edge is cut too, below, but is not stopped: one stop list
      ## per fire, not per pair of fires.) Not when the simulated sizes are returned: they must be whole.
      ## An open pixel with no spreadProb row is outside the crop and never burns, so it is left out.
      stopArgs <- NULL
      if (isTRUE(stopCutFires) && !isTRUE(returnSims)) {
        stopOpen <- openEdge[ids %in% annualFires$ids & pixelID %in% shortAnnDT$pixelID]
        if (NROW(stopOpen))
          stopArgs <- list(stopCells = crop$toCrop(stopOpen$pixelID),
                           stopEvent = match(stopOpen$ids, annualFires$ids), stopAt = 1L)
      }
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
      ## one fire per cell, NA unburnable) but a
      ## plain C++ loop, 2-3x faster on this call shape. It is NOT bit-compatible
      ## with spread(), so cached fits from before this change are not comparable.
      ## The landscape is already cropped above, which matters more for spreadCpp()
      ## than it did for spread().
      ## Per-year random effect (fitYearSpreadSD in .objfunSpreadFit): a seasonal departure. Every
      ## pixel of this year shifts by the same eps on the logit scale, drawn anew in each replicate,
      ## so all of a year's fires are bigger or smaller together; ignition cells stay 1. Still one
      ## spreadCpp() call per replicate.
      ## The year effect is added on the logit of p, except with the powerLaw link, where it is added to the
      ## linear predictor x and p follows from it (p cannot pass p_c, so no year percolates).
      if (yearSpreadSD > 0) {
        pixCrop <- crop$toCrop(shortAnnDT$pixelID)
        powerLaw <- identical(link, "powerLaw")
        lp <- if (powerLaw) as.vector(mat %*% covPars) else stats::qlogis(shortAnnDT$spreadProb)
        toP <- if (powerLaw) function(z) powerLawLink(z, logisticPars[1L], powerLawGamma, lowerSpreadProb, powerLawShape) else stats::plogis
        lociCrop <- crop$toCrop(loci)
      }
      spreadState <- lapply(seq_len(Nreps), function(i) {
        ## each replicate copies the year's spreadProb vector (~70 MB on a large ELF) and R's heap growth lets the dead
        ## copies pile up (peak 5.5 GB on a 2.25 GB worker); a young-generation gc each replicate keeps that down at
        ## no measurable time cost (a full gc here is >15 times slower)
        gc(full = FALSE)
        sp <- spreadProbCrop
        if (yearSpreadSD > 0) {
          sp[pixCrop] <- toP(lp + stats::rnorm(1, 0, yearSpreadSD))
          sp[lociCrop] <- 1
        }
        if (is.null(escapeMinPx)) {
          SpaDES.tools::spreadCpp(
            landscape = crop$r,
            loci = crop$toCrop(loci),
            spreadProb = sp,
            stopCells = stopArgs$stopCells, stopEvent = stopArgs$stopEvent, stopAt = stopArgs$stopAt
          )
        } else { ## an escaped fire burns its first escapeMinPx cells whatever spreadProb
          SpaDES.tools::spreadCpp(
            landscape = crop$r,
            loci = crop$toCrop(loci),
            spreadProb = sp,
            minSize = escapeMinPx,
            jumpTries = jumpTries,
            jumpMeanDist = jumpMeanDist,
            stopCells = stopArgs$stopCells, stopEvent = stopArgs$stopEvent, stopAt = stopArgs$stopAt
          )
        }
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
                    pSummary = pSummary,
                    burned = if (isTRUE(returnBurned)) spreadState[, list(rep, pixelID = indices)]))
      }
      ## One row per simulated fire. A fire that burned an open edge pixel of any of the year's buffers
      ## is cut: it would have kept burning, so its size is "at least this big". The per-fire term gives
      ## it no density at the observed size; the other size-based terms score the size it burned.
      set(spreadState, NULL, "open", spreadState$indices %in% openEdge$pixelID) # one lookup, not one per fire
      simSizes <- spreadState[, list(N = .N, cut = any(open)), by = c("rep", "initialLocus")]
      ret <- append(ret, list(runaways = sum(simSizes$cut), nSims = NROW(simSizes)))
      if (isTRUE(doPitADTest)) {
        ## one PIT per fitted fire, in annualFires order; a fire with no simulated replicate (no row) gets u = 1
        pe <- simSizes[annualFires, on = c("initialLocus" = "cells")]
        ret <- append(ret, list(pit = pe[, list(u = pitMidRank(size[1], N[!is.na(N)], cut[!is.na(N)], Nreps)),
                                         by = "ids"][["u"]]))
      }
      if (isTRUE(doMADTest) || isTRUE(doADTest)) {
        fireSizes <- round(simSizes[, list(N = sum(N)), by = "initialLocus"][["N"]] / Nreps, 0) # Here tabulate() is equivalent to table() but faster
        ret <- append(ret, list(fireSizes = fireSizes))
      }
      if (isTRUE(doADTest)) {
        ## The adTest compares distributions, and the observed sample is single fires, so it gets
        ## every replicate's fire, not `fireSizes`: a per-fire mean over Nreps has a far shorter tail.
        ## fires of 1 pixel never spread; the per-fire likelihood already leaves them out (N > 1 above)
        ret <- append(ret, list(allFireSizes = simSizes[N > 1][["N"]]))
      }
      if (isTRUE(doCrps)) {
        ## per fitted fire, its replicates' sizes against its observed size; a replicate with no row is left out
        pe <- simSizes[annualFires, on = c("initialLocus" = "cells")]
        ret <- append(ret, list(
          crpsFire = pe[, list(s = crpsFireSizes(size[1], N[!is.na(N)], cut[!is.na(N)])), by = "ids"][, sum(s)],
          crpsSizes = simSizes[, list(year = as.character(yr), rep, N)]))
      }
      if (isTRUE(doYearArea)) {
        ## this year's simulated total per replicate (fires of more than 1 pixel); a replicate with none burns 0
        a <- simSizes[N > 1][, list(A = sum(N)), by = "rep"]
        tot <- numeric(Nreps); tot[a$rep] <- a$A
        ret <- append(ret, list(annualAreaByRep = tot))
      }
      if (SpaDES.core::anyPlotting(plot.it)) { # THIS IS PLOTTING STUFF
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
      }
    } else {
      if (isTRUE(returnSims)) {
        f <- data.table::as.data.table(annualFires)[, list(ids, initialLocus = cells, size)]
        return(list(sims = data.table::data.table(yr = as.character(yr), rep = NA_integer_, f,
                                                  sim = NA_integer_), pSummary = pSummary, burned = NULL))
      }
      if (isTRUE(doCrps)) {
        ## a refused year scores as if every replicate left every fire at one pixel
        ret <- append(ret, list(crpsFire = sum(log(annualFires$size)),
                                crpsSizes = data.table::data.table(year = as.character(yr),
                                                                   rep = rep(seq_len(Nreps), each = NROW(annualFires)),
                                                                   N = 1L)))
      }
      # stop("encountered error with spreadProb - contact module developers")
      # Ian added this stop. Unclear what is supposed to happen. Object ret doesn't exist
    }
  }

  return(ret)
}

#' Automatic weight for the Anderson-Darling term
#'
#' The multiplier of the `"adTest"` term, `c * sqrt(nFires)`: it grows with the number of fitted fires,
#' as the (retired) fire-size likelihood the term was balanced against did. `c` was measured across six
#' ELFs (138 to 1624 fires), one value per `sizeLik` and `weighted` combination, on populations fitted
#' with an older objective, so treat the result as the right order of magnitude rather than exact.
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
#' all stop at the same size (e.g. every one burned out at a single pixel) have a standard deviation of 0.
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
#'   as integers, named by their data year. The table used for `yr` is that of the latest data year
#'   strictly before `yr` ([dataYearBefore()]); the function stops if a pixel of `annDTx1000` is
#'   not in it.
#'
#' @param indexNonAnnual Unused, kept for callers: the table of `yr` is chosen from the names of
#'   `nonAnnualDTx1000`, which must be the data years (see [dataYearBefore()]).
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
#'   annual or nonAnnual, that will be used. It may start with `spreadInterceptTxt` (see
#'   [spreadDesignCols()]): that column is made here, as 1s, after rescaling and centring.
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
    ## the table of the latest data year strictly before `yr`, as fireSense_dataPrepFit built it
    tableYears <- as.integer(names(nonAnnualDTx1000))
    whKeep <- match(dataYearBefore(tableYears, yr), tableYears)
    shortAnnDTx1000 <- nonAnnualDTx1000[[whKeep]][annDTx1000, on = "pixelID"]
    covsInTable <- intersect(spreadCovCols(colsToUse), names(nonAnnualDTx1000[[whKeep]]))
    nNA <- sum(!stats::complete.cases(shortAnnDTx1000[, covsInTable, with = FALSE]))
    if (nNA > 0L) {
      stop("fire year ", yr, " has ", nNA, " pixel(s) without covariates in vegetation table \"",
           names(nonAnnualDTx1000)[whKeep], "\": the table and the fire data disagree about which ",
           "vegetation year a fire year uses")
    }
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
  covCols <- spreadCovCols(colsToUse) # the intercept, if any, is not a covariate
  for (cn2 in covCols) {
    set(shortAnnDTx1000, NULL, cn2, shortAnnDTx1000[[cn2]] / 1000)
  }
  
  # rename because it is no longer x1000
  shortAnnDT <- shortAnnDTx1000

  if (doAssertions) {
    assertCovariateRange(shortAnnDT, covCols)
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
    for (cn in intersect(names(covCentre), covCols)) {
      set(shortAnnDT, NULL, cn, shortAnnDT[[cn]] - covCentre[[cn]])
    }
  }

  ## The intercept is a constant 1, added last so it is never rescaled, centred or asserted on.
  if (spreadInterceptTxt %in% colsToUse) set(shortAnnDT, NULL, spreadInterceptTxt, 1)

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

#' Pixels a fire must burn to count as escaped
#'
#' The smallest whole number of pixels whose area is at least `escapeSizeHa`. Shared by the spread
#' fit ([.objfunSpreadFit()]) and fireSense's burning, so "escaped" means the same in both.
#'
#' @param escapeSizeHa Size (ha) a fire must reach to count as escaped.
#' @param landscape A `SpatRaster` (or anything with [terra::res()]) in metres.
#' @return An integer, at least 1.
#' @export
#' @examples
#' r <- terra::rast(nrows = 10, ncols = 10, xmin = 0, xmax = 2400, ymin = 0, ymax = 2400)
#' escapeSizePixels(50, r) # 240 m pixels are 5.76 ha: 9 pixels
escapeSizePixels <- function(escapeSizeHa, landscape) {
  pixHa <- prod(terra::res(landscape)) / 1e4
  max(1L, as.integer(ceiling(escapeSizeHa / pixHa - 1e-9)))
}

#' Negative log-likelihood of each year's observed area burned
#'
#' Scores each fit year's observed total area burned against that year's simulated totals, one per
#' replicate, with the size likelihood the per-fire term uses (`sizeLik`): a Student-t on the square-root
#' scale (`sizeLikT()`) or the empirical kernel density of `sqrt(sims)`. A year with fewer than 2 simulated
#' totals (the objective refused to simulate it) gets `minLik`, as its fires do in the per-fire term.
#'
#' @param obsYear Named numeric: observed total area burned per fit year, in pixels.
#' @param simYear Named list, per year: the simulated totals, one per replicate.
#' @param sizeLik,sizeLikDf As for `.objfunSpreadFit()`.
#' @param minLik Floor on each year's likelihood.
#'
#' @return A single numeric: the summed negative log-likelihood over the years of `obsYear`.
#' @keywords internal
yearAreaNLL <- function(obsYear, simYear, sizeLik = "t", sizeLikDf = 5, minLik = 1e-29) {
  lik <- vapply(names(obsYear), function(y) {
    N <- simYear[[y]]
    if (length(N) < 2) return(minLik)
    if (identical(sizeLik, "t")) sizeLikT(obsYear[[y]], N, sizeLikDf)
    else EnvStats::demp(x = sqrt(obsYear[[y]]), obs = sqrt(N))
  }, numeric(1))
  -sum(log(pmax(minLik, lik)))
}

#' Area-weighted Cramér-von Mises distance between two fire-size samples
#'
#' Compares the share of total area burned that comes from fires up to each size, `F_w(x) =
#' sum(size[size <= x]) / sum(size)`, for the simulated and the observed fires. The squared difference is
#' averaged over the observed fires, each weighted by its share of the observed area, and multiplied by the
#' number of observed fires, so it grows with the sample like the AD statistic. It depends only on the order
#' of sizes, so it is the same on the log scale. Large fires dominate it, as they do area burned.
#'
#' @param sim,obs Numeric vectors of fire sizes (pixels).
#'
#' @return A single non-negative numeric: 0 when the two samples have the same area shares, at most
#'   `length(obs)`, which is also what it returns when there are no simulated fires.
#' @keywords internal
areaWeightedCvM <- function(sim, obs) {
  stopifnot(length(obs) > 0)
  if (!length(sim)) return(length(obs))
  Fw <- function(s, x) {
    s <- sort(s)
    cs <- c(0, cumsum(as.numeric(s)) / sum(s))
    cs[findInterval(x, s) + 1L]
  }
  x <- sort(obs)
  w <- x / sum(x)
  length(obs) * sum(w * (Fw(sim, x) - Fw(obs, x))^2)
}
