#' Wrapper around DEoptim call
#'
#' Does the multiple cluster connections. This will only work if 
#' ssh keys are correctly made between machines (if using multiple machines).
#'
#' @param landscape A RasterLayer which has the correct metadata associated with
#'   the pixelID and cells of other objects in this function call
#' @param annualDTx1000 A list of data.table objects. Each list element will be from 1
#'   year, and it must be the same length as \code{fireBufferedListDT} and \code{hhistoricalFires}
#' @param nonAnnualDTx1000 A list of data.table objects. Each list element must be named 
#'   with a concatenated sequence of names from \code{names(annualDTx1000)}, e.g., \code{1991_1992_1993}.
#'   It should contain all the years in \code{names(annualDTx1000)}
#' @param fireBufferedListDT A list of data.table objects. It must be same length as
#'   \code{annualDTx1000}, with same names. Each element is a data.table with columns:
#'   \code{buff}...
#' @param itermax Passed to \code{DEoptim.control}
#' @param trace Passed to \code{DEoptim.control}
#' @param strategy Passed to \code{DEoptim.control}
#' @param cores A numeric (for running on localhost only) or a character vector of 
#'   machine names (including possibly "localhost"), where
#'   the length of the vector indicates how many cores should be used on that machine.
#' @param logPath A character string indicating what file to write logs to. This 
#'   \code{dirname(logPath)} must exist on each machine, though the function will make sure it
#'   does internally
#' @param cachePath The cachePath to store cache in. Should likely be \code{cachePath(sim)}
#' @param lower Passed to \code{DEoptim}
#' @param upper Passed to \code{DEoptim}
#' @param formula Passed to \code{DEoptim}
#' @param covMinMax Passed to \code{fireSenseUtils::.objfun}
#' @param tests Passed to \code{fireSenseUtils::.objfun}
#' @param maxFireSpread Passed to \code{fireSenseUtils::.objfun}
#' @param Nreps Passed to \code{fireSenseUtils::.objfun}
#' @param .verbose Passed to \code{fireSenseUtils::.objfun}
#' @param visualizeDEoptim Logical. If \code{TRUE}, then histograms will be made of
#'   DEoptim outputs
#'
#' @return
#'
#' @export
#' @importFrom data.table rbindlist as.data.table set
#' @importFrom crayon blurred
#' @importFrom parallel clusterExport clusterEvalQ
#' @importFrom reproducible checkPath 
#' @importFrom future makeClusterPSOCK
runDEoptim <- function(landscape,
                       annualDTx1000,
                       nonAnnualDTx1000,
                       fireBufferedListDT,
                       historicalFires,
                       itermax,
                       trace,
                       strategy,
                       cores,
                       logPath,
                       cachePath,
                       lower,
                       upper,
                       formula,
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
  objsNeeded <- list("landscape",
                     "annualDTx1000",
                     "nonAnnualDTx1000",
                     "fireBufferedListDT",
                     "historicalFires")
  if (!is.null(cores)) {
    message("Starting ", paste(paste(unique(cores)), "x", table(cores),
                               collapse = ", "), " clusters")
    logPath <- file.path(logPath,
                         paste0("fireSense_SpreadFit_log", Sys.getpid()))
    message(crayon::blurred(paste0("Starting parallel model fitting for ",
                                   "fireSense_SpreadFit. Log: ", logPath)))
    
    # Make sure logPath can be written in the workers -- need to create the dir
    
    ## Make cluster with just one worker per machine --> don't need to do these steps
    #  multiple times per machine
    if (is.numeric(cores)) cores <- rep("localhost", cores)
    revtunnel <- if (all(cores == "localhost")) FALSE else TRUE
    st <- system.time(cl <- future::makeClusterPSOCK(unique(cores), revtunnel = revtunnel))
    clusterExport(cl, list("logPath"), envir = environment())
    
    parallel::clusterEvalQ(
      cl, {
        reproducible::checkPath(dirname(logPath), create = TRUE)
        devtools::install_github("PredictiveEcology/fireSenseUtils@development")
      }
    )
    stopCluster(cl)
    
    
    st <- system.time(cl <- future::makeClusterPSOCK(cores, revtunnel = revtunnel, outfile = logPath))
    
    on.exit(stopCluster(cl))
    message("it took ", round(st[3],2), "s to start ",
            paste(paste(unique(cores)), "x", table(cores),
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
    browser()
    
  } else {
    list2env(mget(unlist(objsNeeded), envir = environment()), envir = .GlobalEnv)
  }
  #####################################################################
  # DEOptim call
  #####################################################################
  browser()
  DE <- Cache(DEoptimIterative, itermax = itermax, lower = lower,
              upper = upper,
              control = do.call("DEoptim.control", control),
              formula = formula,
              covMinMax = covMinMax,
              # tests = c("mad", "SNLL_FS"),
              tests = c("SNLL_FS"),
              maxFireSpread = maxFireSpread,
              Nreps = Nreps,
              .verbose = .verbose,
              omitArgs = c("verbose"))#,
              #cacheId = "cd495b412420ad4a") # iteration 201 to 300
  DE
}

#' Make histograms of DEoptim object
#' 
#' @export
#' @param DE An object from a \code{DEoptim} call
#' @param cachePath A \code{cacheRepo} to pass to \code{showCache} and
#'   \code{loadFromCache} if \code{DE} is missing.
visualizeDE <- function(DE, cachePath) {
  cacheID <- if (missing(DE)) {
    sc <- showCache(userTags = "DEoptim")
    tail(sc$cacheId, 1)
  } else {
    gsub("cacheId:", "", grep("cacheId", attr(DE, "tags"), value = TRUE))
  }
  DE <- reproducible::loadFromCache(cachePath, cacheId = cacheID)
  aa <- as.data.table(t(DE$member$pop))
  aa[, pars := paste0("par", 1:NROW(aa))];
  dim1 <- floor(sqrt(NROW(aa)))
  dim2 <- NROW(aa) / dim1
  par(mfrow = c(dim1,dim2)); 
  aa[, hist(t(.SD), main = as.character(.BY))[[2]], by = pars]
}

#' @inheritParams runDEoptim
#' @export
#' @rdname runDEoptim
DEoptimIterative <- function(itermax, lower,
                             upper,
                             control,
                             formula,
                             covMinMax,
                             tests,
                             maxFireSpread,
                             Nreps,
                             .verbose) {
  data.table::setDTthreads(1)
  iterStep <- 10
  for (iter in seq_len(itermax / iterStep) * iterStep) {
    control$itermax <- iterStep
    browser()
    st1 <- system.time(DE <- Cache(DEoptim,
                                   fireSenseUtils::.objfun,
                                   lower = lower,
                                   upper = upper,
                                   control = do.call("DEoptim.control", control),
                                   formula = formula,
                                   covMinMax = covMinMax,
                                   # tests = c("mad", "SNLL_FS"),
                                   tests = c("SNLL_FS"),
                                   maxFireSpread = maxFireSpread,
                                   Nreps = Nreps,
                                   verbose = .verbose,
                                   omitArgs = c("verbose")
    ))
    control$initialpop <- DE$member$pop
    if (isTRUE(visualizeDEoptim)) {
      visualizeDE(DE, cachePath)
    }
  }
  
}