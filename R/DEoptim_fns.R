utils::globalVariables(c(
  ".BY", ".SD", "bestValue", "value", "iter", "lower95", "upper95", "var",
  "dif", "variable", "pred"
))

#' Run `DEoptim`
#'
#' Provides a wrapper around [DEoptim::DEoptim], setting up the multiple cluster connections.
#' This will only work if ssh keys are preconfigured on all machines (if using multiple machines).
#'
#' @param landscape A `SpatRaster` which has the correct metadata associated with
#'   the `pixelID` and cells of other objects in this function call.
#'
#' @param annualDTx1000 A list of data.table objects. Each list element will be from 1
#'   year, and it must be the same length as `fireBufferedListDT` and `historicalFires`.
#'   All covariates must be integers, and must be `1000x` their actual values.
#'
#' @param nonAnnualDTx1000 A list of data.table objects. Each list element must be named
#'   with a concatenated sequence of names from `names(annualDTx1000)`,
#'   e.g., `1991_1992_1993`.
#'   It should contain all the years in `names(annualDTx1000)`.
#'   All covariates must be integers, and must be `1000x` their actual values.
#'
#' @param fireBufferedListDT A list of data.table objects. It must be same length as
#'   `annualDTx1000`, with same names. Each element is a `data.table` with columns:
#'   `buff`...TODO: INCOMPLETE
#'
#' @param historicalFires A named list of `data.frame`s (one per year, names
#'   matching those of `annualDTx1000`), each with columns `cells` (pixel
#'   indices of fire ignitions) and `size` (fire size in pixels). Used as the
#'   observed reference against which simulated fires are scored.
#'
#' @param itermax Maximum number of iterations for the [DEoptim::DEoptim] algorithm.
#'   Passed to [DEoptim::DEoptim.control].
#'
#' @param initialpop Optional. A matrix or vector specifying the initial population
#'   for [DEoptim::DEoptim]. If `NULL`, [DEoptim::DEoptim] generates one.
#'   Passed to [DEoptim::DEoptim.control].
#'
#' @param NP Optional. The number of population members (individuals) in [DEoptim::DEoptim].
#'   If `NULL`, [DEoptim::DEoptim] sets a default. [DEoptim::DEoptim.control].
#'
#' @param trace Integer or Logical. Controls the level of tracing information
#'   printed by [DEoptim::DEoptim] during optimization. Passed to [DEoptim::DEoptim.control].
#'
#' @param strategy Integer `[1,10]`. Defines the [DEoptim::DEoptim] strategy variant to use.
#'   Passed to [DEoptim::DEoptim.control].
#'
#' @param cores A numeric (for running on localhost only) or a character vector of
#'   machine names (including possibly "localhost"), where
#'   the length of the vector indicates how many cores should be used on that machine.
#'
#' @param libPath A character string indicating an R package library directory.
#'   This location must exist on each machine, though the function will make sure it
#'   does internally.
#'
#' @param logPath A character string indicating what file to write logs to. This
#'   `dirname(logPath)` must exist on each machine, though the function will make sure it
#'   does internally.
#'
#' @param doObjFunAssertions logical indicating whether to do assertions.
#'
#' @param paths list of paths containing the `cachePath` to store cache.
#'    Should likely be `cachePath(sim)`. See [SpaDES.core::paths].
#'
#' @param iterStep Integer. Must be less than `itermax`. This will cause [DEoptim::DEoptim] to run
#'   the `itermax` iterations in `ceiling(itermax / iterStep)` steps. At the end of
#'   each step, this function will plot, optionally, the parameter histograms (if
#'   `visualizeDEoptim` is `TRUE`)
#'
#' @param lower Numeric vector. Lower bounds for the parameters being optimized.
#'   Passed to [DEoptim::DEoptim].
#'
#'   If named and including `yearSpreadSD` (which must be last), the per-year random effect (a seasonal departure) of
#'   [.objfunSpreadFit()] (`fitYearSpreadSD`) is fitted, in the fit and in the re-score.
#'
#' @param upper Numeric vector. Upper bounds for the parameters being optimized.
#'   Passed to [DEoptim::DEoptim].
#'
#' @template mutuallyExclusive
#'
#' @param formulaToFit Passed to [DEoptim::DEoptim]
#'
#' @param objFunCoresInternal Integer. The number of cores to use for potential
#'   parallelization *within* a single call to the objective function ([.objfunSpreadFit()]).
#'   This is distinct from the parallelization managed by [DEoptim::DEoptim] across population members.
#'
#' @param covMinMax,tests,maxFireSpread,Nreps,.verbose Passed to [.objfunSpreadFit()].
#'
#' @param covCentre Passed to [.objfunSpreadFit()] in the fit, the re-score and the profile and
#'   simulation diagnostics: values subtracted from the rescaled covariates, from [spreadCovCentre()]
#'   when the formula has an intercept. `NULL` (default) centres nothing, and then the call, and so its
#'   cache key, is what it was before this argument existed.
#'
#' @param thresh Threshold multiplier used in SNLL fire size (`"snll_fs"`) test. Default 550.
#'
#' @param visualizeDEoptim Logical. If `TRUE`, then histograms will be made of [DEoptim::DEoptim] outputs.
#'
#' @param plotEvery Integer. Generations between DEoptim progress figures; the final figures are
#'   always drawn. Passed to `clusters::DEoptimIterative()`. It is not part of the fit's cache key,
#'   so changing it does not refit. Default 25.
#'
#' @param .plotSize List specifying plot `height` and `width`, in pixels.
#'
#' @param rep Integer. An identifier for the replication number of this optimization run.
#'   Used in cache tags and plot filenames. Default 1L.
#'
#' @param .plots Character string. Specifies the plot destination device (e.g., "screen", "png", "pdf").
#'   Passed to internal plotting functions (likely via [SpaDES.core::Plots()]).
#'   Default "screen".
#'
#' @param runName Character string used to label this run. Forwarded to
#'   `clusters::DEoptimIterative()` and used as a suffix for the cache `.functionName`
#'   so that runs with different `runName` values get distinct cache entries.
#'   Default `""` (no suffix).
#'
#' @param nCoresNeeded Integer. How many workers to request for the DEoptim cluster; defaults to
#'   about 10 per estimated parameter, `10 * length(lower)`. DEoptim's `NP` is set to the number of
#'   workers the cluster actually gets, so a smaller allocation means a smaller population.
#'
#' @param rescoreReps Integer. After the fit, each member of the final population is evaluated this
#'   many more times (full evaluations, no early stop) and the result is attached to the returned
#'   object as `attr(DE, "finalRescore")`, for [bestByReplicatedMean()]. `0` skips it.
#' @param profileReps Integer. If above 0, [profileCoefficients()] is run around the best member
#'   of the re-score, each point evaluated this many times, and attached as `attr(DE, "profile")`.
#'   Needs `rescoreReps > 0`. Costs about `6 * nCoefficients * profileReps` evaluations.
#' @param simulateMembers Integer. If above 0, that many best members of the re-score simulate the
#'   observed fires without the size cap ([simulateFireSizes()]), attached as `attr(DE, "fitSims")`
#'   for [scoreFireSizes()] and [linkSaturation()]. Needs `rescoreReps > 0`.
#' @param link Passed to [.objfunSpreadFit()] in the fit and the re-score: `NULL` or
#'   `"logistic3pUpper"`. Needed because DEoptim may hand the objective an unnamed `par`.
#' @param jumpTries,jumpMeanDist,yearAreaWeight,areaDistWeight Passed to [.objfunSpreadFit()], in the fit
#'   and the final-population re-score. Off by default (`0`).
#' @param penaliseRunaways Passed to [.objfunSpreadFit()] in the fit and the re-score: `TRUE` (default)
#'   scores a simulated fire that reaches the edge of its own buffer as a runaway, not as a fire of
#'   the size it reached. The edge ring of every buffer is computed here, once ([addBufferEdge()]).
#' @param runawayEdgeFrac,runawayEdgeMin Passed to [.objfunSpreadFit()] in the fit and the re-score: a
#'   replicate is a runaway only if it burns at least `max(runawayEdgeMin, ceiling(runawayEdgeFrac * n))`
#'   of the `n` edge-ring pixels of its fire (at most `n`). Defaults [fireSenseRunawayEdgeFrac] and [fireSenseRunawayEdgeMin].
#' @param penaliseCapHits Deprecated and ignored; use `penaliseRunaways`.
#' @param escapeSizeHa Passed to [.objfunSpreadFit()] in the fit and the re-score: the size (ha) a
#'   fire must reach to count as escaped. `NULL` keeps the historical rules.
#' @param sizeLik,sizeLikDf,weighted,adWeight Passed to [.objfunSpreadFit()], in the fit AND in the
#'   final-population re-score, so both evaluate the same objective. The defaults are that function's.
#'   Before these existed, every fit used those defaults whatever the caller wanted.
#' @param .c Numeric in `(0, 1]`. DEoptim's `c`, the speed of crossover adaptation (used by
#'   `strategy = 6`). Passed to [DEoptim::DEoptim.control()]; a `c` in `DEoptimControl` wins.
#' @param DEoptimControl Named list of further [DEoptim::DEoptim.control()] settings (for example
#'   `CR`, `F`, `p`, `reltol`), passed through [clusters::clusterSetup()] to DEoptim. `NP` is the
#'   number of workers the cluster gets.
#'
#' @return The result of the `clusters::DEoptimIterative()` call. This is typically a list where
#' each element contains the [DEoptim::DEoptim] object state after a block of `iterStep` iterations.
#' The final element represents the state after `itermax` iterations or upon early stopping.
#'
#' @export
#' @importFrom clusters clusterSetup
#' @importFrom crayon blurred
#' @importFrom data.table rbindlist as.data.table set
#' @importFrom parallel clusterExport clusterEvalQ stopCluster
#' @importFrom parallelly makeClusterPSOCK
#' @importFrom reproducible Cache checkPath
#' @importFrom SpaDES.core P
#' @importFrom RhpcBLASctl blas_get_num_procs blas_set_num_threads omp_get_max_threads omp_set_num_threads
#' @importFrom utils install.packages installed.packages packageVersion
runDEoptim <- function(landscape,
                       annualDTx1000,
                       nonAnnualDTx1000,
                       fireBufferedListDT,
                       historicalFires,
                       itermax,
                       initialpop = NULL,
                       NP = NULL,
                       trace,
                       strategy,
                       cores = NULL,
                       paths,
                       libPath = .libPaths()[1],
                       logPath = tempfile(sprintf(
                         "runDEoptim_%s_",
                         format(Sys.time(), "%Y-%m-%d_%H%M%S")
                       ), fileext = ".log"),
                       doObjFunAssertions = getOption("fireSenseUtils.assertions", TRUE),
                       iterStep = 25,
                       lower,
                       upper,
                       mutuallyExclusive,
                       formulaToFit,
                       objFunCoresInternal,
                       covMinMax = covMinMax,
                       tests = c("SNLL", "adTest"),
                       maxFireSpread,
                       Nreps,
                       thresh = 550,
                       .c = 0.5,
                       DEoptimControl = list(),
                       .verbose,
                       visualizeDEoptim = logPath,
                       .plots = "screen",
                       plotEvery = 25L,
                       .plotSize = list(height = 1600, width = 2000),
                       rep = 1L,
                       runName = "",
                       nCoresNeeded = 10L * length(lower),
                       rescoreReps = 10L,
                       sizeLik = "kde",
                       sizeLikDf = 5,
                       weighted = TRUE,
                       adWeight = "auto",
                       link = NULL,
                       escapeSizeHa = NULL,
                       jumpTries = 0,
                       jumpMeanDist = 0,
                       yearAreaWeight = 0,
                       areaDistWeight = 0,
                       penaliseRunaways = TRUE,
                       runawayEdgeFrac = fireSenseRunawayEdgeFrac,
                       runawayEdgeMin = fireSenseRunawayEdgeMin,
                       profileReps = 0L,
                       simulateMembers = 0L,
                       penaliseCapHits = NULL,
                       covCentre = NULL) {
  if (!is.null(penaliseCapHits)) deprecatedCapArgs("penaliseCapHits")
  ## the edge ring of each fire's buffer, once, before the tables go to the workers
  fireBufferedListDT <- addBufferEdge(fireBufferedListDT, landscape)
  if (isTRUE(is.na(cores))) cores <- NULL
  origBlas <- blas_get_num_procs()
  if (origBlas > 1) {
    blas_set_num_threads(1)
    on.exit(blas_set_num_threads(origBlas), add = TRUE)
  }
  origOmp <- omp_get_max_threads()
  ## NA where R has no OpenMP (macOS CRAN builds)
  if (isTRUE(origOmp > 1)) {
    omp_set_num_threads(1)
    on.exit(omp_set_num_threads(origOmp), add = TRUE)
  }

  ####################################################################
  #  Cluster
  ####################################################################
  objsNeeded <- list(
    "landscape",
    "annualDTx1000",
    "nonAnnualDTx1000",
    "fireBufferedListDT",
    "historicalFires",
    "mutuallyExclusive"
  )

  neededPkgs <- c("magrittr", "raster", "data.table", "SpaDES.core",
                  "SpaDES.tools", "fireSenseUtils", "sf", "plyr",# "mirai",
                  "munsell")

  control <- clusters::clusterSetup(
    messagePrefix = as.character(rep), # .runName,
    strategy = strategy, itermax = itermax,
    cores = cores, # logPath = file.path(dataPath(sim)),
    ## about 10 workers per estimated parameter; clusterSetup() sets NP to the workers it gets
    nCoresNeeded = nCoresNeeded,
    libPath = libPath[1], NP = NP,
    logPath = logPath,
    objsNeeded = objsNeeded,
    pkgsNeeded = neededPkgs, envir = environment(),
    ## every DEoptim setting the caller gave reaches DEoptim; .c is DEoptim's c
    controlArgs = utils::modifyList(list(c = .c), as.list(DEoptimControl))
  )
  cl <- control$cluster # This is to test whether it is actually closed
  
  #####################################################################
  # DEOptim call
  #####################################################################
  ## name the parameters from `lower`. termsInDEoptim() called every non-formula parameter a "logit"
  ## term, so yearSpreadSD was reported as a third logistic parameter.
  covTerms <- spreadDesignCols(formulaToFit) # includes the intercept, if the formula has one
  logisticTerms <- setdiff(names(lower), c(covTerms, yearSpreadSDTxt))
  message("Fitting ", length(lower), " parameters: logistic: ", paste(logisticTerms, collapse = ", "),
          "; covariates: ", paste(covTerms, collapse = ", "),
          if (yearSpreadSDTxt %in% names(lower)) paste0("; year effect: ", yearSpreadSDTxt))
  message("objectiveFunction threshold SNLL to run all years after first 2 years: ", thresh)
  ## the per-year random effect is fitted when the bounds include it (its sd is the LAST parameter);
  ## DEoptim passes `par` unnamed, so the objective is told explicitly
  fitYearSpreadSD <- yearSpreadSDTxt %in% names(lower)
  if (fitYearSpreadSD && !identical(names(lower)[length(lower)], yearSpreadSDTxt))
    stop("runDEoptim: `", yearSpreadSDTxt, "` must be the last element of `lower` and `upper`")

  # aaaa <<- 1; on.exit(rm(aaaa, envir = .GlobalEnv))
  DE <- Cache(
    clusters::DEoptimIterative(
      fn = fireSenseUtils::.objfunSpreadFit,
      # DE <- Cache(
      #   DEoptimIterative(
      itermax = itermax,
      lower = lower,
      upper = upper,
      ## only what was set here (NP from the built cluster, strategy, ...); DEoptimIterative()
      ## fills the rest, so passing a complete DEoptim.control() would override its defaults
      control = control,
      formulaToFit = formulaToFit,
      covMinMax = covMinMax,
      covCentre = covCentre,
      # tests = c("mad", "SNLL_FS"),
      tests = tests,
      figurePath = visualizeDEoptim,
      maxFireSpread = maxFireSpread,
      objFunCoresInternal = objFunCoresInternal,
      Nreps = Nreps,
      .verbose = .verbose,
      mutuallyExclusive = mutuallyExclusive,
      doAssertions = doObjFunAssertions,
      # visualizeDEoptim = visualizeDEoptim,
      .plots = .plots,
      plotEvery = plotEvery,
      .plotSize = .plotSize,
      iterStep = iterStep,
      thresh = thresh,
      sizeLik = sizeLik,
      sizeLikDf = sizeLikDf,
      weighted = weighted,
      adWeight = adWeight,
      link = link,
      fitYearSpreadSD = fitYearSpreadSD,
      escapeSizeHa = escapeSizeHa,
      jumpTries = jumpTries,
      jumpMeanDist = jumpMeanDist,
      yearAreaWeight = yearAreaWeight,
      areaDistWeight = areaDistWeight,
      penaliseRunaways = penaliseRunaways,
      runawayEdgeFrac = runawayEdgeFrac,
      runawayEdgeMin = runawayEdgeMin,
      rep = rep,
      runName = runName),
    cachePath = paths$cachePath,
    ## how often progress figures are drawn does not change the fit
    omitArgs = c(".verbose", "plotEvery", omitNullArgs(covCentre = covCentre)),
    ## The data the objective runs on reaches the workers through clusterSetup(objsNeeded), not as an
    ## argument above, so it is not in this key by itself: two held-out folds of one ELF shared their
    ## whole fit (2026-09-29). clusterSetup() digested those objects once; use that digest.
    .cacheExtra = clusters::shippedObjectsDigest(control)
    , .functionName = paste0("DEoptimIterative_", runName)
    # , cacheId = "8448b6a37b54361b"
  ) # iteration 201 to 300

  ## The fit moves workers (rebalancing, dead-worker rebuilds) in its own copy of the cluster and
  ## stops the old ones, so `cl` may hold closed connections now: "invalid connection" in the
  ## rescore (ELF 5.2.1 fold 2, 2026-10-03). Use the cluster the workers are on now.
  cl <- clusters::currentCluster(cl)

  ## Re-score the final population `rescoreReps` times each, on the same workers, with the early stop
  ## off (thresh = Inf) so every replicate is a full evaluation. DEoptim's own values are partly luck;
  ## the caller picks its best members from these means (bestByReplicatedMean()).
  finalPop <- if (length(DE)) DE[[length(DE)]]$member$pop
  rescoreArgs <- list(formulaToFit = formulaToFit, covMinMax = covMinMax, tests = tests,
                      maxFireSpread = maxFireSpread, objFunCoresInternal = objFunCoresInternal,
                      Nreps = Nreps, mutuallyExclusive = mutuallyExclusive, doAssertions = FALSE,
                      thresh = Inf, verbose = 0, sizeLik = sizeLik, sizeLikDf = sizeLikDf,
                      weighted = weighted, adWeight = adWeight, link = link,
                      fitYearSpreadSD = fitYearSpreadSD, escapeSizeHa = escapeSizeHa,
                      jumpTries = jumpTries, jumpMeanDist = jumpMeanDist,
                      yearAreaWeight = yearAreaWeight, areaDistWeight = areaDistWeight,
                      penaliseRunaways = penaliseRunaways,
                      runawayEdgeFrac = runawayEdgeFrac, runawayEdgeMin = runawayEdgeMin)
  rescoreArgs$covCentre <- covCentre # NULL adds nothing, so the cached re-score keeps its key
  if (isTRUE(rescoreReps > 0) && !is.null(finalPop)) {
    colnames(finalPop) <- names(lower)
    attr(DE, "finalRescore") <- Cache(
      rescorePopulation,
      pop = finalPop, fn = fireSenseUtils::.objfunSpreadFit, reps = as.integer(rescoreReps), cl = cl,
      fnArgs = rescoreArgs,
      omitArgs = "cl",
      cachePath = paths$cachePath,
      .functionName = paste0("rescoreFinalPopulation_", runName))
  }

  ## Diagnostics that need the workers, while they still hold the data (see ?fitDiagnostics)
  if (!is.null(attr(DE, "finalRescore")) && isTRUE(profileReps > 0 || simulateMembers > 0)) {
    top <- bestByReplicatedMean(finalPop, attr(DE, "finalRescore"), n = max(1L, simulateMembers))
    if (isTRUE(profileReps > 0))
      attr(DE, "profile") <- Cache(
        profileCoefficients,
        best = unlist(top$params[1]), pop = finalPop, fn = fireSenseUtils::.objfunSpreadFit,
        reps = as.integer(profileReps), cl = cl, fnArgs = rescoreArgs,
        omitArgs = "cl", cachePath = paths$cachePath,
        .functionName = paste0("profileBestMember_", runName))
    if (isTRUE(simulateMembers > 0))
      attr(DE, "fitSims") <- Cache(
        simulateFireSizes,
        pop = top$params, fn = fireSenseUtils::.objfunSpreadFit, cl = cl, fnArgs = rescoreArgs,
        omitArgs = "cl", cachePath = paths$cachePath,
        .functionName = paste0("simulateBestMembers_", runName))
  }
  DE
}

#' Make histograms of `DEoptim` object `pars`
#'
#' @param DE An object from a [DEoptim::DEoptim] call
#' @param cachePath A `cacheRepo` to pass to `showCache` and `loadFromCache` if `DE` is missing.
#' @param titles titles of plots
#' @param lower lower limit on x axis
#' @param upper upper limit on x axis
#'
#' @export
#' @importFrom data.table as.data.table
#' @importFrom graphics hist par
#' @importFrom reproducible loadFromCache showCache
#' @importFrom utils tail
#' @importFrom ggplot2 coord_cartesian ggtitle xlab theme_minimal
visualizeDE <- function(DE, cachePath, titles, lower, upper) {
  if (missing(DE)) {
    if (missing(cachePath)) {
      stop("Must provide either DE or cachePath")
    }
    message("DE not supplied; visualizing the most recent added to Cache")
    sc <- showCache(userTags = "DEoptim")
    cacheID <- tail(sc$cacheId, 1)
    DE <- reproducible::loadFromCache(cachePath, cacheId = cacheID)
  }
  if (is(DE, "list")) {
    DE <- tail(DE, 1)[[1]]
  }

  cc <- as.data.table(DE$member$pop)
  setnames(cc, titles)
  suppressWarnings(bb <- melt(cc))
  ff <- lapply(titles, function(p) {
    ggplot(bb[variable == p], aes(value)) +
      geom_histogram(bins = 15) + coord_cartesian(xlim = c(lower[p],upper[p])) +
      ggtitle(p) + xlab(NULL) +
      theme_minimal()
  })
  invisible(cowplot::plot_grid(plotlist = ff))
}

#' `termsInDEoptim`
#'
#' `termsInDEoptim()` is deprecated: it counted every parameter not in the formula as a "logit" term, so the per-year
#' random effect `yearSpreadSD` was reported as a third logistic parameter. `runDEoptim()` now
#' names the parameters with `names(lower)`.
#'
#' @param fireSense_spreadFormula The formula to be submitted to [DEoptim::DEoptim()],
#'                                from e.g., `sim$fireSense_spreadFormula`.
#'
#' @param thresh The threshold for accepting fits; e.g., from `mod$thresh`.
#'
#' @param numParams The number of parameters (TODO: improve description)
#'
#' @export
#' @rdname runDEoptim
termsInDEoptim <- function(fireSense_spreadFormula, thresh, numParams) {
  .Deprecated(msg = paste0("fireSenseUtils::termsInDEoptim() is deprecated; runDEoptim() names ",
                           "the parameters with names(lower)"))
  termsInForm <- attr(terms(as.formula(fireSense_spreadFormula, env = .GlobalEnv)), "term.labels")
  logitNumParams <- numParams - length(termsInForm)
  message("Using a ", logitNumParams, " parameter logistic equation")
  message("  There will be ", logitNumParams, " logit terms & ", numParams, " terms in all:")
  message("  ", paste(c(paste0("logit", seq(logitNumParams)), termsInForm), collapse = ", "))
  message("  objectiveFunction threshold SNLL to run all years after first 2 years: ", thresh)
  c(paste0("logit", seq(logitNumParams)), termsInForm)
}
