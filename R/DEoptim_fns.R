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
#' @param historicalFires TODO: DESCRIPTION NEEDED
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
#' @param thresh Threshold multiplier used in SNLL fire size (`"snll_fs"`) test. Default 550.
#'
#' @param visualizeDEoptim Logical. If `TRUE`, then histograms will be made of [DEoptim::DEoptim] outputs.
#'
#' @param .plotSize List specifying plot `height` and `width`, in pixels.
#'
#' @return The result of the [DEoptimIterative()] call. This is typically a list where
#' each element contains the [DEoptim::DEoptim] object state after a block of `iterStep` iterations.
#' The final element represents the state after `itermax` iterations or upon early stopping.
#'
#' @export
#' @importFrom clusters clusterSetup
#' @importFrom crayon blurred
#' @importFrom data.table rbindlist as.data.table set
#' @importFrom parallel clusterExport clusterEvalQ stopCluster
#' @importFrom parallelly makeClusterPSOCK
#' @importFrom qs qread qsave
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
                         "fireSense_SpreadFit_%s_",
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
                       .verbose,
                       visualizeDEoptim = logPath,
                       .plots = "screen",
                       .plotSize = list(height = 1600, width = 2000),
                       rep = 1L) {
  if (isTRUE(is.na(cores))) cores <- NULL
  origBlas <- blas_get_num_procs()
  if (origBlas > 1) {
    blas_set_num_threads(1)
    on.exit(blas_set_num_threads(origBlas), add = TRUE)
  }
  origOmp <- omp_get_max_threads()
  if (origOmp > 1) {
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

  neededPkgs <- c("kSamples", "magrittr", "raster", "data.table", "SpaDES.core",
                  "SpaDES.tools", "fireSenseUtils", "sf", # "mirai",
                  "munsell")

  control <- clusters::clusterSetup(
    messagePrefix = as.character(rep), # .runName,
    strategy = strategy, itermax = itermax,
    cores = cores, # logPath = file.path(dataPath(sim)),
    libPath = libPath[1], NP = NP,
    logPath = logPath,
    objsNeeded = objsNeeded,
    pkgsNeeded = neededPkgs, envir = environment()
  )
  cl <- control$cluster # This is to test whether it is actually closed
  # on.exit(parallel::stopCluster(control$cluster), add = TRUE)

  # control <- list(itermax = itermax, trace = trace, strategy = strategy)

  # if (!is.null(initialpop)) {
  #   control$initialpop <- initialpop
  # }
  #
  # if (!is.null(NP)) {
  #   control$NP <- NP
  # }
  # # if (is.null(cores)) cores <- "localhost"
  # if (!is.null(cores)) {
  #   logPath <- file.path(
  #     logPath,
  #     paste0(
  #       "fireSense_SpreadFit_", format(Sys.time(), "%Y-%m-%d_%H%M%S"),
  #       "_pid", Sys.getpid(), ".log"
  #     )
  #   )
  #   message(paste0(
  #     "Starting parallel model fitting for ",
  #     "fireSense_SpreadFit. Log: ", logPath
  #   ))
  #
  #   # Make sure logPath can be written in the workers -- need to create the dir
  #
  #   if (is.numeric(cores)) cores <- rep("localhost", cores)
  #
  #   ## Make cluster with just one worker per machine --> don't need to do these steps
  #   #     multiple times per machine, if not all 'localhost'
  #   revtunnel <- FALSE
  #   neededPkgs <- c("kSamples", "magrittr", "raster", "data.table",
  #                   "SpaDES.tools", "fireSenseUtils", "sf")
  #
  #   if (!identical("localhost", unique(cores))) {
  #     repos <- c("https://predictiveecology.r-universe.dev", getOption("repos"))
  #
  #     aa <- Require::pkgDep(unique(c("dqrng", "PredictiveEcology/SpaDES.tools@development",
  #                                    "PredictiveEcology/fireSenseUtils@development", "qs",
  #                                    "RCurl", neededPkgs)), recursive = TRUE)
  #     pkgsNeeded <- unique(Require::extractPkgName(unname(unlist(aa))))
  #
  #     revtunnel <- ifelse(all(cores == "localhost"), FALSE, TRUE)
  #
  #     coresUnique <- setdiff(unique(cores), "localhost")
  #     message("copying packages to: ", paste(coresUnique, collapse = ", "))
  #
  #     # RscriptPath = "/usr/local/bin/Rscript"
  #
  #     # st <- system.time(
  #     #   cl <- mirai::make_cluster(
  #     #     n = length(coresUnique),
  #     #     url = "tcp://localhost:5563",
  #     #     remote = mirai::ssh_config(
  #     #       remotes = paste0("ssh://", coresUnique),
  #     #       tunnel = TRUE,
  #     #       timeout = 1,
  #     #       rscript = RscriptPath
  #     #     )
  #     #   )
  #     # )
  #     st <- system.time({
  #       cl <- parallelly::makeClusterPSOCK(coresUnique, revtunnel = revtunnel, rscript_libs = libPath
  #                                          # , rscript = c("nice", RscriptPath)
  #       )
  #     })
  #     clusterExport(cl, list("libPath", "logPath", "repos", "pkgsNeeded"),
  #                   envir = environment())
  #
  #     # Missing `dqrng` and `sitmo`
  #     Require::Install(pkgsNeeded, libPaths = libPath)
  #
  #     parallel::clusterEvalQ(cl, {
  #       # If this is first time that packages need to be installed for this user on this machine
  #       #   there won't be a folder present that is writable
  #       if (!dir.exists(libPath)) {
  #         dir.create(libPath, recursive = TRUE)
  #       }
  #     })
  #
  #     message("Setting up packages on the cluster...")
  #     out <- lapply(setdiff(unique(cores), "localhost"), function(ip) {
  #       rsync <- Sys.which("rsync")
  #       if (!nzchar(rsync)) stop()
  #       system(paste0(rsync, " -aruv --update ", paste(file.path(libPath, pkgsNeeded), collapse = " "),
  #                     " ", ip, ":", libPath))
  #     })
  #
  #     parallel::clusterEvalQ(cl, {
  #       # If this is first time that packages need to be installed for this user on this machine
  #       #   there won't be a folder present that is writable
  #       if (tryCatch(packageVersion("Require") < "1.0.1", error = function() TRUE))
  #         install.packages("Require", lib = libPath)
  #       library(Require, lib.loc = libPath)
  #       dir.create(dirname(logPath), recursive = TRUE)
  #       out <- Require::Install(pkgsNeeded, libPaths = libPath)
  #     })
  #
  #     GDALversions <- parallel::clusterEvalQ(cl, {
  #       .libPaths(libPath)
  #       return(sf::sf_extSoftVersion()["GDAL"])
  #     })
  #
  #     stopifnot(length(unique(sf::sf_extSoftVersion()["GDAL"], GDALversions)) == 1)
  #     parallel::stopCluster(cl)
  #
  #     ## Now make full cluster with one worker per core listed in "cores"
  #     message("Starting ", paste(paste(names(table(cores))), "x", table(cores),
  #                                collapse = ", "), " clusters")
  #     message("Starting main parallel cluster ...")
  #     # sshCores <- paste0("ssh//", grep('localhost', cores, invert = TRUE, value = TRUE))
  #     # nonsshCores <- grep('localhost', cores, value = TRUE)
  #     # coresForMirai <- c(nonsshCores, sshCores)
  #
  #     st <- system.time({
  #       # cl <- mirai::make_cluster(
  #       #   length(coresForMirai),
  #       #   # url = "tcp://localhost:5555",
  #       #   remote = ssh_config(
  #       #     remotes = coresForMirai,
  #       #     # tunnel = TRUE,
  #       #     timeout = 1,
  #       #     rscript = RscriptPath
  #       #   )
  #       # )
  #
  #       cl <- parallelly::makeClusterPSOCK(cores,
  #                                          revtunnel = revtunnel,
  #                                          outfile = logPath, rscript_libs = libPath
  #                                          # , rscript = c("nice", RscriptPath)
  #       )
  #     })
  #
  #     on.exit(stopCluster(cl))
  #     message(
  #       "it took ", round(st[3], 2), "s to start ",
  #       paste(paste(names(table(cores))), "x", table(cores), collapse = ", "), " threads"
  #     )
  #     message("Moving objects to each node in cluster")
  #
  #     stMoveObjects <- try({
  #       system.time({
  #         objsToCopy <- mget(unlist(objsNeeded))
  #         objsToCopy <- lapply(objsToCopy, FUN = function(x) {
  #           if (inherits(x, "SpatRaster")) {
  #             x <- terra::wrap(x)
  #           } else {
  #             x
  #           }
  #           x
  #         })
  #         filenameForTransfer <- normalizePath(tempfile(fileext = ".qs"), mustWork = FALSE, winslash = "/")
  #         dir.create(dirname(filenameForTransfer), recursive = TRUE, showWarnings = FALSE) # during development, this was deleted accidentally
  #         qs::qsave(objsToCopy, file = filenameForTransfer)
  #         stExport <- system.time({
  #           outExp <- clusterExport(cl, varlist = "filenameForTransfer", envir = environment())
  #         })
  #         out11 <- clusterEvalQ(cl, {
  #           dir.create(dirname(filenameForTransfer), recursive = TRUE, showWarnings = FALSE)
  #         })
  #         out <- lapply(setdiff(unique(cores), "localhost"), function(ip) {
  #           rsync <- Sys.which("rsync")
  #           st1 <- system.time(system(paste0(rsync, " -av ",
  #                                            filenameForTransfer, " ", ip, ":",
  #                                            filenameForTransfer)))
  #         })
  #         out <- clusterEvalQ(cl, {
  #           out <- qs::qread(file = filenameForTransfer)
  #           out <- lapply(out, FUN = function(x) {
  #             if (inherits(x, "PackedSpatRaster")) {
  #               x <- terra::unwrap(x)
  #             } else {
  #               x
  #             }
  #             x
  #           })
  #           list2env(out, envir = .GlobalEnv)
  #         })
  #         # Delete the file
  #         out <- clusterEvalQ(cl, {
  #           if (dir.exists(dirname(filenameForTransfer))) {
  #             try(unlink(dirname(filenameForTransfer), recursive = TRUE), silent = TRUE)
  #           }
  #         })
  #       })
  #     })
  #
  #     if (is(stMoveObjects, "try-error")) {
  #       message("The attempt to move objects to cluster using rsync and qs failed; trying clusterExport")
  #       stMoveObjects <- system.time(clusterExport(cl, objsNeeded, envir = environment()))
  #       list2env(mget(unlist(objsNeeded), envir = environment()), envir = .GlobalEnv)
  #     }
  #     message("it took ", round(stMoveObjects[3], 2), "s to move objects to nodes")
  #     message("loading packages in cluster nodes")
  #
  #     clusterExport(cl, "neededPkgs", envir = environment())
  #     stPackages <- system.time(parallel::clusterEvalQ(
  #       cl,
  #       {
  #         for (i in neededPkgs) {
  #           library(i, character.only = TRUE)
  #         }
  #         message("loading ", i, " at ", Sys.time())
  #       }
  #     ))
  #     message("it took ", round(stPackages[3], 2), "s to load packages")
  #
  #     control$cluster <- cl
  #   } else {
  #     list2env(mget(unlist(objsNeeded), envir = environment()), envir = .GlobalEnv)
  #   }
  # }
  #####################################################################
  # DEOptim call
  #####################################################################
  termsInDEoptim(formulaToFit, thresh, length(lower))

  DE <- Cache(
    clusters:::DEoptimIterative2(
      fn = fireSenseUtils::.objfunSpreadFit,
      # DE <- Cache(
      #   DEoptimIterative(
      itermax = itermax,
      lower = lower,
      upper = upper,
      control = do.call("DEoptim.control", control),
      formulaToFit = formulaToFit,
      covMinMax = covMinMax,
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
      .c = .c,
      .plotSize = .plotSize,
      iterStep = iterStep,
      thresh = thresh,
      rep = rep),
    cachePath = paths$cachePath,
    omitArgs = c(".verbose")
    #,
    # cacheId = "cd495b412420ad4a"
  ) # iteration 201 to 300
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
    # geom_smooth(se = TRUE) +
    # geom_ribbon(aes(ymin = lower95, ymax = upper95)) +
    # facet_wrap(facets = "variable", scales = "free")
  })
  invisible(ggpubr::ggarrange(plotlist = ff))
  #
  #
  #
  #   aa <- as.data.table(t(DE$member$pop))
  #   melt(aa, id.vars = )
  #   aa[, pars := paste0("par", 1:NROW(aa))]
  #   dim1 <- floor(sqrt(NROW(aa)))
  #   dim2 <- NROW(aa) / dim1
  #   par(mfrow = c(dim1, dim2))
  #   aa[, hist(t(.SD), main = as.character(.BY))[[2]], by = pars]
}

#' Iterative `DEoptim` Runner with Caching and Visualization
#'
#' Internal function called by `runDEoptim`. Runs [DEoptim::DEoptim] in steps,
#' caching results and optionally visualizing progress after each step.
#'
#' @param rep Integer. An identifier for the replication number of this optimization run.
#'   Used in cache tags and plot filenames. Default 1L.
#'
#' @param .plots Character string. Specifies the plot destination device (e.g., "screen", "png", "pdf").
#'   Passed to internal plotting functions (likely via [SpaDES.core::Plots()]).
#'   Default "screen".
#'
#' @template mutuallyExclusive
#'
#' @param .c TODO: DESCRIPTION NEEDED
#'
#' @param figPath directory where figures will be saved, if relevant
#'
#' @param rep TODO: DESCRIPTION NEEDED
#'
#' @param control passed to [DEoptim::DEoptim.control]
#'
#' @param cachePath A `cacheRepo` (see [reproducible::Cache()]).
#'
#' @export
#' @importFrom crayon green
#' @importFrom data.table rbindlist setDTthreads
#' @importFrom DEoptim DEoptim DEoptim.control
#' @importFrom grDevices dev.off png
#' @importFrom ggplot2 geom_abline
#' @importFrom quickPlot isRstudioServer
#' @importFrom reproducible Cache isUpdated messageDF
#' @importFrom stats dnorm rnorm
#' @importFrom utils tail
#' @rdname runDEoptim
DEoptimIterative <- function(itermax, lower, upper,
                             control, formulaToFit, covMinMax,
                             tests = c("SNLL", "adTest"),
                             objFunCoresInternal,
                             maxFireSpread,
                             Nreps,
                             visualizeDEoptim,
                             figPath,
                             cachePath,
                             mutuallyExclusive,
                             doObjFunAssertions = getOption("fireSenseUtils.assertions", TRUE),
                             iterStep = 25,
                             thresh = 550,
                             .c = 0.5,
                             .verbose,
                             .plots = "screen",
                             .plotSize = list(height = 1600, width = 2000),
                             rep = 1L) {
  data.table::setDTthreads(1)
  message("starting DEoptimIFterative at ", Sys.time())
  x1 <- rnorm(1e2, 1, 2) # this is for debugging below
  DE <- list()

  itersToDo <- seq_len(ceiling(itermax / iterStep))
  cacheIds <- rep(NA_real_, itermax)

  if (FALSE) {
    sc <- showCache(Function = "DEoptimForCache")
    sc <- sc[tagKey == "function"]
    sc[, time := as.POSIXct(createdDate)]
    setorder(sc, time)
    cacheIdsFromCache <- unique(sc$cacheId)#[-(1:3)]
    cacheIds[seq_along(cacheIdsFromCache)] <- cacheIdsFromCache
    cacheIds <- na.omit(cacheIds)
    itersToDo <- min( (sum(!is.na(cacheIds))), itermax):itermax
  }

  for (iter in itersToDo) {
    control$itermax <- pmin(iterStep, itermax - iterStep * (iter - 1))
    control$storepopfrom <- control$itermax + 1
    control$reltol <- 0.1
    control$c <- .c

    controlArgs <- do.call("DEoptim.control", control)
    controlForCache <- controlArgs[c(
      "VTR", "strategy", "NP", "CR", "F", "bs", "trace",
      "initialpop", "p", "c", "reltol",
      "packages", "parVar", "foreachArgs"
    )]

    if (TRUE) {
      DE[[iter]] <- Cache(
        DEoptimForCache(
          fireSenseUtils::.objfunSpreadFit,
          lower = lower,
          upper = upper,
          control = controlArgs,
          formulaToFit = formulaToFit,
          covMinMax = covMinMax,
          tests = tests,
          maxFireSpread = maxFireSpread,
          mutuallyExclusive = mutuallyExclusive,
          doAssertions = doObjFunAssertions,
          Nreps = Nreps,
          plot.it = FALSE,
          controlForCache = controlForCache,
          objFunCoresInternal = objFunCoresInternal,
          thresh = thresh),
        cacheId = cacheIds[iter],
        .functionName = paste0("DEoptimForCache_", rep),
        verbose = .verbose,
        omitArgs = c("verbose", "control")
      )
      if (!isUpdated(DE[[iter]]))
        message(paste(round(unname(DE[[iter]]$optim$bestmem), 4), collapse = " "))
      message(crayon::green("Iteration ", iter, " done!"))
    } else {
      # This is for testing --> it is fast
      fn <- function(par, x) {
        -sum(dnorm(log = TRUE, x, mean = par[1], sd = par[2]))
      }

      st1 <- system.time(DE[[iter]] <- Cache(DEoptimForCache,
                                             fn,
                                             lower = lower,
                                             upper = upper,
                                             mutuallyExclusive = mutuallyExclusive,
                                             controlForCache = controlForCache,
                                             control = control,
                                             omitArgs = c("verbose", "control"),
                                             x = x1
      ))
    }

    control$initialpop <- DE[[iter]]$member$pop

    rng <- 25; # do 25 iteration steps, i.e., 1:100, 25:125
    dataRunToUse <- 150 # this will do the lm on this many items
    numSegments <- (length(DE) - dataRunToUse) / rng + 1# (length(DE) - dataRunToUse + 1) / rng
    pvals <- c(0,0)


    # Do these here because we need them for both sections below
    dfForGGplotSimple <- DEoptimToDataFrame(DE)
    gg1 <- ggPlotFnSimple(dfForGGplotSimple)

    if (numSegments > 1) {
      isNewSegment <- numSegments %% 1 == 0
      if (isNewSegment) {
        pvals <- numeric(floor(numSegments))
        iters <- list()
        s <- list()
        l <- list()
        segmentSeq <- seq_len(floor(numSegments))
        # if (!exists("dfForGGplotSimple", inherits = FALSE))
        for (i in segmentSeq) {
          col <- "black"
          if (i == tail(segmentSeq, 2)[1]) col <- "blue"
          if (i == tail(segmentSeq, 1)[1]) col <- "red"
          iters[[i]] <- seq_len(dataRunToUse) + (i-1) * rng;
          a <- data.table(iter = seq_along(DE), val = sapply(DE, function(x) x$member$bestvalit))
          l[[i]] <- lm(val ~ iter, data = a[iters[[i]]]);
          s[[i]] <- summary(l[[i]]);
          pvals[i] <- round(s[[i]]$coefficients[2, 4], 4)

          newdat <- data.table(iter = iters[[i]])
          set(newdat, NULL, "pred", predict(l[[i]], newdata = newdat))
          int <- s[[i]]$coefficients[1, 1]
          slop <- s[[i]]$coefficients[2, 1]
          # gg1 <- gg1 + geom_line(data = newdat,
          #                           aes(x = iter, y = pred), #, xend = tail(iter, 1), yend = tail(pred, 1)),
          #                           col = col)
          gg1 <- gg1 + geom_abline(intercept = int, slope = slop,
                                   #                        aes(x = iter, y = pred), #, xend = tail(iter, 1), yend = tail(pred, 1)),
                                   col = col)
        }
        pvalDT <- data.table(dataRange = sapply(segmentSeq, function(x) paste(range(iters[[x]]), collapse = ":")),
                             pvals = pvals)
        # Plots(gg1, types = .plots,
              # filename = ggDEoptimFilename(visualizeDEoptim, rep, text = "objFun/"))
        messageDF(pvalDT, colour = "yellow")
      }
    }
    if (!isFALSE(visualizeDEoptim) && (isUpdated(DE[[iter]]))) { # i.e., should be a path
      terms <- suppressMessages(termsInDEoptim(formulaToFit, thresh, length(lower)))
      nVars <- NCOL(DE[[iter]]$member$pop)
      if (length(terms) != nVars )
        terms <- c(terms, paste0("V", seq(nVars - length(terms))))
      dfForGGplot <- visualizeDEoptimLines(DE, terms = terms)
      dfForGGplotAllPoints <- visualizeDEoptimLines(DE, terms = terms, allPoints = TRUE)
      dfForGGplotSimple <- DEoptimToDataFrame(DE)


      withCallingHandlers({
        Plots(gg1, types = .plots,
              filename = ggDEoptimFilename(visualizeDEoptim, rep, text = "objFun/"))
        #Plots(dfForGGplotSimple, ggPlotFnSimple, types = .plots,
        #      filename = ggDEoptimFilename(visualizeDEoptim, rep, text = "objFun/"))
        Plots(dfForGGplotAllPoints, ggPlotFnMeansAllPoints, types = .plots,
              filename = ggDEoptimFilename(visualizeDEoptim, rep, text = "lines_mean_AllPoints/"));
        Plots(dfForGGplot, ggPlotFnMeans, types = .plots,
              filename = ggDEoptimFilename(visualizeDEoptim, rep, text = "lines_mean/"))
        Plots(dfForGGplot, ggPlotFnDif, types = .plots, ,
              filename = ggDEoptimFilename(visualizeDEoptim, rep, text = "lines_dif/"))
        Plots(dfForGGplot, ggPlotFnVars, types = .plots, ,
              filename = ggDEoptimFilename(visualizeDEoptim, rep, text = "lines_variance/"))
        Plots(fn = visualizeDE, DE = DE[[iter]], cachePath = cachePath,
              titles = terms, lower = lower, upper = upper, types = .plots,
              filename = ggDEoptimFilename(visualizeDEoptim, rep = rep, iter = iter, text = "hists/", time = TRUE))
      }, message = function(m) {
        if (any(grepl("geom_smooth|SavingSaved", m$message)))
          invokeRestart("muffleMessage")
      })
      reproducible::messageColoured(colour = "green",
                                    "5 Figures saved to: ", dirname(ggDEoptimFilename("~", 1, text = "")),
                                    verbose = .verbose)

    }


    # Break out if the last N segments are "non-significant slope at p == 0.1 i.e., conservative
    if (all(tail(pvals, 2) > 0.1) && length(DE) > 349) {
      break
    }
  }

  DE
}

#' @importFrom DEoptim DEoptim
DEoptimForCache <- function(...) {
  dots <- list(...)
  dots["controlForCache"] <- NULL
  do.call(DEoptim, dots)
}

#' `termsInDEoptim`
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
  termsInForm <- attr(terms(as.formula(fireSense_spreadFormula, env = .GlobalEnv)), "term.labels")
  logitNumParams <- numParams - length(termsInForm)
  message("Using a ", logitNumParams, " parameter logistic equation")
  message("  There will be ", logitNumParams, " logit terms & ", numParams, " terms in all:")
  message("  ", paste(c(paste0("logit", seq(logitNumParams)), termsInForm), collapse = ", "))
  message("  objectiveFunction threshold SNLL to run all years after first 2 years: ", thresh)
  c(paste0("logit", seq(logitNumParams)), termsInForm)
}

#' @importFrom stats setNames
DEoptimToDataFrame <- function(d, item = "bestvalit") {
  b <- lapply(d, function(dr) as.data.frame(dr$member[[item]]) |> setNames("bestValue"))
  b <- rbindlist(b, idcol = "iter")
  b
}

visualizeDEoptimLines <- function(d, terms, allPoints = FALSE) {
  iter <- length(d)
  se <- seq(iter)
  # this commented code will use "all the population
  if (isTRUE(allPoints)) {
    b <- lapply(d, function(dr) as.data.frame(dr$member$pop))
    b <- rbindlist(b, idcol = "iter")#
    setnames(b, old = grep("^V", colnames(b), value = TRUE),  terms)

  } else {
    b <- do.call(rbind, lapply(d, function(dr) colMeans(dr$member$pop))) |> as.data.table()
    setnames(b, terms)
    blower <- do.call(rbind, lapply(d, function(dr)
      sapply(seq(NCOL(dr$member$pop)), function(x) quantile(dr$member$pop[, x], 0.025)))) |>
      as.data.table()
    bupper <- do.call(rbind, lapply(d, function(dr)
      sapply(seq(NCOL(dr$member$pop)), function(x) quantile(dr$member$pop[, x], 0.975)))) |>
      as.data.table()
    bvar <- do.call(rbind, lapply(d, function(dr)
      sapply(seq(NCOL(dr$member$pop)), function(x) var(dr$member$pop[, x])))) |>
      as.data.table()
    setnames(blower, names(b))
    setnames(bupper, names(b))
    setnames(bvar, names(b))
    b[, iter := se]
    blower[, iter := se]
    bupper[, iter := se]
    bvar[, iter := se]
    blower <- melt(blower, id.vars = "iter")
    setnames(blower, old = "value", new = "lower95")
    bupper <- melt(bupper, id.vars = "iter")
    setnames(bupper, old = "value", new = "upper95")
    bvar <- melt(bvar, id.vars = "iter")
    setnames(bvar, old = "value", new = "var")


  }

  b <- melt(b, id.vars = "iter")
  if (isTRUE(allPoints)) {
    bmerged <- b
  } else {
    ons <- c("iter", "variable")
    bmerged <- b[blower, on = ons][bupper, on = ons][bvar, on = ons]
    bmerged[, dif := upper95 - lower95]
  }
  bmerged[]
}


ggPlotFnMeans <- function(bmerged) {
  ggplot(bmerged, aes(iter, value)) +
    geom_point() +
    geom_smooth(se = TRUE) +
    # geom_ribbon(aes(ymin = lower95, ymax = upper95)) +
    facet_wrap(facets = "variable", scales = "free")
}

#' @importFrom ggplot2 geom_point geom_smooth
ggPlotFnSimple <- function(bmerged) {
  ggplot(bmerged, aes(iter, bestValue)) +
    geom_point() +
    geom_smooth(se = TRUE)
}

#' @importFrom ggplot2 geom_point geom_smooth
ggPlotFnDif <- function(bmerged) {
  ggplot(bmerged, aes(iter, dif)) +
    geom_point() +
    geom_smooth(se = TRUE) +
    # geom_ribbon(aes(ymin = lower95, ymax = upper95)) +
    facet_wrap(facets = "variable", scales = "free")
}

#' @importFrom ggplot2 geom_point geom_smooth
ggPlotFnVars <- function(bmerged) {
  ggplot(bmerged, aes(iter, var)) +
  geom_point() +
  geom_smooth(se = TRUE) +
  # geom_ribbon(aes(ymin = lower95, ymax = upper95)) +
  facet_wrap(facets = "variable", scales = "free")
}

#' @importFrom ggplot2 geom_smooth geom_jitter ggplot
ggPlotFnMeansAllPoints <- function(b) {
  ggplot(b, aes(iter, value)) +
  # geom_point() +
  geom_jitter(size = 0.05, width = 0.2, col = "grey") +
  geom_smooth(se = TRUE) +
  # geom_ribbon(aes(ymin = lower95, ymax = upper95)) +
  facet_wrap(facets = "variable", scales = "free")
}


#' @importFrom reproducible paddedFloatToChar
ggDEoptimFilename <- function(visualizeDEoptim, rep, iter = NULL, text = "DE_hists_", time = FALSE) {
  file.path(visualizeDEoptim,
            "fireSense_SpreadFit",
            paste0(text, "rep", paddedFloatToChar(rep, padL = 3),
                   ifelse(is.null(iter), "", paste0("_iter", iter)), "_", Sys.getpid(),
                   ifelse(isTRUE(time), paste0("_", as.character(round(Sys.time(), 0))), ""), ".png"))
}
