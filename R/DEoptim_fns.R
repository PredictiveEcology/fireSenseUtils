utils::globalVariables(c(
  ".BY", ".SD", "pars"
))

#' Wrapper around \code{DEoptim} call
#'
#' Does the multiple cluster connections. This will only work if
#' ssh keys are correctly made between machines (if using multiple machines).
#'
#' @param landscape A \code{RasterLayer} which has the correct metadata associated with
#'   the \code{pixelID} and cells of other objects in this function call
#' @param annualDTx1000 A list of data.table objects. Each list element will be from 1
#'   year, and it must be the same length as \code{fireBufferedListDT} and \code{historicalFires}
#' @param nonAnnualDTx1000 A list of data.table objects. Each list element must be named
#'   with a concatenated sequence of names from \code{names(annualDTx1000)},
#'   e.g., \code{1991_1992_1993}.
#'   It should contain all the years in \code{names(annualDTx1000)}.
#' @param fireBufferedListDT A list of data.table objects. It must be same length as
#'   \code{annualDTx1000}, with same names. Each element is a \code{data.table} with columns:
#'   \code{buff}...TODO: INCOMPLETE
#' @param historicalFires DESCRIPTION NEEDED
#' @param itermax Passed to \code{DEoptim.control}
#' @param initialpop DESCRIPTION NEEDED
#' @param NP DESCRIPTION NEEDED
#' @param trace Passed to \code{DEoptim.control}
#' @param strategy Passed to \code{DEoptim.control}
#' @param cores A numeric (for running on localhost only) or a character vector of
#'   machine names (including possibly "localhost"), where
#'   the length of the vector indicates how many cores should be used on that machine.
#' @param logPath A character string indicating what file to write logs to. This
#'   \code{dirname(logPath)} must exist on each machine, though the function will make sure it
#'   does internally.
#' @param cachePath The \code{cachePath} to store cache in. Should likely be \code{cachePath(sim)}
#' @param iterStep Integer. Must be less than \code{itermax}. This will cause \code{DEoptim} to run
#'   the \code{itermax} iterations in \code{ceiling(itermax / iterStep)} steps. At the end of
#'   each step, this function will plot, optionally, the parameter histograms (if
#'   \code{visualizeDEoptim} is \code{TRUE})
#' @param lower Passed to \code{DEoptim}
#' @param upper Passed to \code{DEoptim}
#' @param FS_formula Passed to \code{DEoptim}
#' @param objFunCoresInternal DESCRIPTION NEEDED
#' @param covMinMax Passed to \code{fireSenseUtils::.objfun}
#' @param tests Passed to \code{fireSenseUtils::.objfun}
#' @param maxFireSpread Passed to \code{fireSenseUtils::.objfun}
#' @param Nreps Passed to \code{fireSenseUtils::.objfun}
#' @param .verbose Passed to \code{fireSenseUtils::.objfun}
#' @param visualizeDEoptim Logical. If \code{TRUE}, then histograms will be made of
#'   \code{DEoptim} outputs
#'
#' @return DESCRIPTION NEEDED
#'
#' @export
#' @importFrom crayon blurred
#' @importFrom data.table rbindlist as.data.table set
#' @importFrom future makeClusterPSOCK
#' @importFrom parallel clusterExport clusterEvalQ stopCluster
#' @importFrom reproducible Cache checkPath
runDEoptim <- function(landscape,
                       annualDTx1000,
                       nonAnnualDTx1000,
                       fireBufferedListDT,
                       historicalFires,
                       itermax,
                       initialpop,
                       NP,
                       trace,
                       strategy,
                       cores,
                       logPath,
                       cachePath,
                       iterStep = 25,
                       lower,
                       upper,
                       FS_formula,
                       objFunCoresInternal,
                       covMinMax = covMinMax,
                       tests,
                       maxFireSpread,
                       Nreps,
                       .verbose,
                       visualizeDEoptim) {
  ####################################################################
  #  Cluster
  ####################################################################
  control <- list(itermax = itermax,
                  trace = trace,
                  strategy = strategy)#,
  if (!is.null(initialpop))
    control$initialpop <- initialpop
  if (!is.null(NP))
    control$NP <- NP
  objsNeeded <- list("landscape",
                     "annualDTx1000",
                     "nonAnnualDTx1000",
                     "fireBufferedListDT",
                     "historicalFires")
  if (!is.null(cores)) {
    message("Starting ", paste(paste(names(table(cores))), "x", table(cores),
                               collapse = ", "), " clusters")
    logPath <- file.path(logPath,
                         paste0("fireSense_SpreadFit_log", Sys.getpid()))
    message(crayon::blurred(paste0("Starting parallel model fitting for ",
                                   "fireSense_SpreadFit. Log: ", logPath)))

    # Make sure logPath can be written in the workers -- need to create the dir

    ## Make cluster with just one worker per machine --> don't need to do these steps
    #  multiple times per machine
    browser()
    if (is.numeric(cores)) cores <- rep("localhost", cores)
    revtunnel <- if (all(cores == "localhost")) FALSE else TRUE
    browser()
    st <- system.time({
      cl <- future::makeClusterPSOCK(unique(cores), revtunnel = revtunnel)
    })
    clusterExport(cl, list("logPath"), envir = environment())

    parallel::clusterEvalQ(
      cl, {
        reproducible::checkPath(dirname(logPath), create = TRUE)
        devtools::install_github("PredictiveEcology/fireSenseUtils@development", upgrade = FALSE)
      }
    )
    stopCluster(cl)
    st <- system.time({
      cl <- future::makeClusterPSOCK(cores, revtunnel = revtunnel, outfile = logPath)
    })

    on.exit(stopCluster(cl))
    message("it took ", round(st[3],2), "s to start ",
            paste(paste(names(table(cores))), "x", table(cores),
                  collapse = ", "), " threads")
    clusterExport(cl, objsNeeded, envir = environment())
    parallel::clusterEvalQ(
      cl, {
        for (i in c("kSamples", "magrittr", "raster", "data.table",
                    "SpaDES.tools", "fireSenseUtils"))
          library(i, character.only = TRUE)
      }
    )
    control$cluster <- cl

  } else {
    list2env(mget(unlist(objsNeeded), envir = environment()), envir = .GlobalEnv)
  }

  #####################################################################
  # DEOptim call
  #####################################################################
  DE <- Cache(DEoptimIterative, itermax = itermax, lower = lower,
              upper = upper,
              control = do.call("DEoptim.control", control),
              FS_formula = FS_formula,
              covMinMax = covMinMax,
              # tests = c("mad", "SNLL_FS"),
              tests = c("SNLL_FS"),
              maxFireSpread = maxFireSpread,
              objFunCoresInternal = objFunCoresInternal,
              Nreps = Nreps,
              .verbose = .verbose,
              visualizeDEoptim = visualizeDEoptim,
              cachePath = cachePath,
              iterStep = iterStep,
              omitArgs = c("verbose"))#,
              #cacheId = "cd495b412420ad4a") # iteration 201 to 300
  DE
}

#' Make histograms of \code{DEoptim} object \code{pars}
#'
#' @param DE An object from a \code{DEoptim} call
#' @param cachePath A \code{cacheRepo} to pass to \code{showCache} and
#'        \code{loadFromCache} if \code{DE} is missing.
#'
#' @export
#' @importFrom data.table as.data.table
#' @importFrom graphics hist par
#' @importFrom reproducible loadFromCache showCache
#' @importFrom utils tail
visualizeDE <- function(DE, cachePath) {
  if (missing(DE)) {
    if (missing(cachePath))
      stop("Must provide either DE or cachePath")
    message("DE not supplied; visulizing the most recent added to Cache")
    sc <- showCache(userTags = "DEoptim")
    cacheID <- tail(sc$cacheId, 1)
    DE <- reproducible::loadFromCache(cachePath, cacheId = cacheID)
  }
  if (is(DE, "list")) {
    DE <- tail(DE, 1)[[1]]
  }

  aa <- as.data.table(t(DE$member$pop))
  aa[, pars := paste0("par", 1:NROW(aa))];
  dim1 <- floor(sqrt(NROW(aa)))
  dim2 <- NROW(aa) / dim1
  par(mfrow = c(dim1, dim2));
  aa[, hist(t(.SD), main = as.character(.BY))[[2]], by = pars]
}

#' @param control DESCRIPTION NEEDED
#'
#' @export
#' @importFrom data.table rbindlist setDTthreads
#' @importFrom DEoptim DEoptim DEoptim.control
#' @importFrom grDevices dev.off png
#' @importFrom quickPlot isRstudioServer
#' @importFrom reproducible Cache
#' @importFrom stats dnorm rnorm
#' @importFrom utils tail
#' @rdname runDEoptim
DEoptimIterative <- function(itermax,
                             lower,
                             upper,
                             control,
                             FS_formula,
                             covMinMax,
                             tests,
                             objFunCoresInternal,
                             maxFireSpread,
                             Nreps,
                             visualizeDEoptim,
                             cachePath,
                             iterStep = 25,
                             .verbose) {
  data.table::setDTthreads(1)
  x1 <- rnorm(1e2, 1, 2) # this is for debugging below
  DE <- list()
  for (iter in seq_len(ceiling(itermax / iterStep))) {
    control$itermax <- pmin(iterStep, itermax - iterStep * (iter - 1))
    control$storepopfrom <- control$itermax + 1

    if (TRUE) {
      controlArgs <- do.call(DEoptim.control, control)
      controlForCache <- controlArgs[c("VTR", "strategy", "NP", "CR", "F", "bs", "trace",
                                       "initialpop", "p", "c", "reltol",
                                       "packages", "parVar", "foreachArgs")]
      st1 <- system.time(DE[[iter]] <- Cache(DEoptimForCache,
        fireSenseUtils::.objfun,
        lower = lower,
        upper = upper,
        control = controlArgs,
        FS_formula = FS_formula,
        covMinMax = covMinMax,
        # tests = c("mad", "SNLL_FS"),
        tests = c("SNLL_FS"),
        maxFireSpread = maxFireSpread,
        Nreps = Nreps,
        controlForCache = controlForCache,
        objFunCoresInternal = objFunCoresInternal,
        verbose = .verbose,
        omitArgs = c("verbose", "control")
      ))
    } else {
      # This is for testing --> it is fast
      fn <- function(par, x) {
        -sum(dnorm(log = TRUE, x, mean = par[1], sd = par[2]))
      }
      controlArgs <- do.call("DEoptim.control", control)
      controlForCache <- controlArgs[c("VTR", "strategy", "NP", "CR", "F", "bs", "trace",
                                       "initialpop", "p", "c", "reltol",
                                       "packages", "parVar", "foreachArgs")]
      st1 <- system.time(DE[[iter]] <- Cache(DEoptimForCache,
        fn,
        lower = lower,
        upper = upper,
        controlForCache = controlForCache,
        control = control,
        omitArgs = c("verbose", "control"),
        x = x1
      ))
    }

    control$initialpop <- DE[[iter]]$member$pop
    if (isTRUE(visualizeDEoptim)) {
      if (!isRstudioServer()) {
        png(filename = paste0("DE_pars", as.character(Sys.time()), "_", Sys.getpid(), ".png"),
            width = 800, height = 1000)
      }
      visualizeDE(DE[[iter]], cachePath)
      if (!isRstudioServer()) {
        dev.off()
      }
    }
  }

  DE1 <- tail(DE, 1)[[1]]
  if (iter > 1) {
    bestvals <- which.min(unlist(lapply(DE, function(x) x$optim$bestval)))
    DE1$optim$bestmem <- DE[[bestvals]]$optim$bestmem
    DE1$optim$bestval <- DE[[bestvals]]$optim$bestval
    DE1$optim$iter <- sum(unlist(lapply(DE, function(x) x$optim$iter)))
    DE1$member$bestmemit <- as.matrix(rbindlist(lapply(DE, function(x) as.data.table(x$member$bestmemit))))
    DE1$member$bestvalit <- rbindlist(lapply(DE, function(x) as.data.table(x$member$bestvalit)))[[1]]
  }
  # DE1$member <- as.matrix(rbindlist(lapply(DE, function(x) as.data.table(x$member$bestmemit))))

  DE
}

#' @importFrom DEoptim DEoptim
DEoptimForCache <- function(...) {
  dots <- list(...)
  dots["controlForCache"] <- NULL
  do.call(DEoptim, dots)
}
