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
#'   year, and it must be the same length as \code{fireBufferedListDT} and \code{historicalFires}.
#'   All covariates must be integers, and must be 1000x their actual values.
#' @param nonAnnualDTx1000 A list of data.table objects. Each list element must be named
#'   with a concatenated sequence of names from \code{names(annualDTx1000)},
#'   e.g., \code{1991_1992_1993}.
#'   It should contain all the years in \code{names(annualDTx1000)}.
#'   All covariates must be integers, and must be 1000x their actual values.
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
#' @param libPath A character string indicating an R package library directory. This
#'   location must exist on each machine, though the function will make sure it
#'   does internally.
#' @param logPath A character string indicating what file to write logs to. This
#'   \code{dirname(logPath)} must exist on each machine, though the function will make sure it
#'   does internally.
#' @param doObjFunAssertions logical indicating whether to do assertions.
#' @param cachePath The \code{cachePath} to store cache in. Should likely be \code{cachePath(sim)}
#' @param iterStep Integer. Must be less than \code{itermax}. This will cause \code{DEoptim} to run
#'   the \code{itermax} iterations in \code{ceiling(itermax / iterStep)} steps. At the end of
#'   each step, this function will plot, optionally, the parameter histograms (if
#'   \code{visualizeDEoptim} is \code{TRUE})
#' @param lower Passed to \code{DEoptim}
#' @param upper Passed to \code{DEoptim}
#' @template mutuallyExclusive
#' @param FS_formula Passed to \code{DEoptim}
#' @param objFunCoresInternal DESCRIPTION NEEDED
#' @param covMinMax Passed to \code{fireSenseUtils::.objfunSpreadFit}
#' @param tests Passed to \code{fireSenseUtils::.objfunSpreadFit}
#' @param maxFireSpread Passed to \code{fireSenseUtils::.objfunSpreadFit}
#' @param Nreps Passed to \code{fireSenseUtils::.objfunSpreadFit}
#' @param thresh Threshold multiplier used in SNLL fire size (SNLL_FS) test. Default 550.
#' @param .verbose Passed to \code{fireSenseUtils::.objfunSpreadFit}
#' @param visualizeDEoptim Logical. If \code{TRUE}, then histograms will be made of
#'   \code{DEoptim} outputs.
#' @param .plotSize List specifying plot \code{height} and \code{width}, in pixels.
#'
#' @return DESCRIPTION NEEDED
#'
#' @export
#' @importFrom crayon blurred
#' @importFrom data.table rbindlist as.data.table set
#' @importFrom parallel clusterExport clusterEvalQ stopCluster
#' @importFrom parallelly makeClusterPSOCK
#' @importFrom qs qread qsave
#' @importFrom reproducible Cache checkPath
#' @importFrom RhpcBLASctl blas_get_num_procs blas_set_num_threads omp_get_max_threads omp_set_num_threads
#' @importFrom utils install.packages packageVersion
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
                       libPath = .libPaths()[1],
                       logPath = tempfile(sprintf("fireSense_SpreadFit_%s_",
                                          format(Sys.time(), "%Y-%m-%d_%H%M%S")), fileext = ".log"),
                       doObjFunAssertions = getOption("fireSenseUtils.assertions", TRUE),
                       cachePath,
                       iterStep = 25,
                       lower,
                       upper,
                       mutuallyExclusive,
                       FS_formula,
                       objFunCoresInternal,
                       covMinMax = covMinMax,
                       tests = c("SNLL", "adTest"),
                       maxFireSpread,
                       Nreps,
                       thresh = 550,
                       .verbose,
                       visualizeDEoptim,
                       .plotSize = list(height = 1600, width = 2000)) {
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
  control <- list(itermax = itermax, trace = trace, strategy = strategy)

  if (!is.null(initialpop)) {
    control$initialpop <- initialpop
  }

  if (!is.null(NP)) {
    control$NP <- NP
  }

  objsNeeded <- list(
    "landscape",
    "annualDTx1000",
    "nonAnnualDTx1000",
    "fireBufferedListDT",
    "historicalFires",
    "mutuallyExclusive"
  )

  if (!is.null(cores)) {
    logPath <- file.path(
      logPath,
      paste0(
        "fireSense_SpreadFit_", format(Sys.time(), "%Y-%m-%d_%H%M%S"),
        "_pid", Sys.getpid(), ".log"
      )
    )
    message(crayon::blurred(paste0(
      "Starting parallel model fitting for ",
      "fireSense_SpreadFit. Log: ", logPath
    )))

    # Make sure logPath can be written in the workers -- need to create the dir

    if (is.numeric(cores)) cores <- rep("localhost", cores)

    ## Make cluster with just one worker per machine --> don't need to do these steps
    #     multiple times per machine, if not all 'localhost'
    revtunnel <- FALSE
    if (!identical("localhost", unique(cores))) {
      revtunnel <- ifelse(all(cores == "localhost"), FALSE, TRUE)

      coresUnique <- setdiff(unique(cores), "localhost")
      message(
        "Making sure packages with sufficient versions installed and loaded on: ",
        paste(coresUnique, collapse = ", ")
      )
      st <- system.time({
        cl <- parallelly::makeClusterPSOCK(coresUnique, revtunnel = revtunnel, rscript_libs = libPath)
      })
      packageVersionFSU <- packageVersion("fireSenseUtils")
      packageVersionST <- packageVersion("SpaDES.tools")
      clusterExport(cl, list("libPath", "logPath", "packageVersionFSU", "packageVersionST"), envir = environment())

      parallel::clusterEvalQ(
        cl,
        {
          ## Use the binary packages for install if Ubuntu & Linux
          if (Sys.info()["sysname"] == "Linux" && grepl("Ubuntu", utils::osVersion)) {
            .os.version <- strsplit(system("lsb_release -c", intern = TRUE), ":\t")[[1]][[2]]
            .user.agent <- paste0(
              "R/", getRversion(), " R (",
              paste(getRversion(), R.version["platform"], R.version["arch"], R.version["os"]),
              ")"
            )
            optsNew <- list(
              "repos" = c(CRAN = paste0("https://packagemanager.rstudio.com/all/__linux__/",
                                        .os.version, "/latest")),
              "HTTPUserAgent" = .user.agent
            )
            opts <- options(optsNew)
            on.exit(options(opts), add = TRUE)
          }

          # If this is first time that packages need to be installed for this user on this machine
          #   there won't be a folder present that is writable
          if (!dir.exists(libPath)) {
            dir.create(libPath, recursive = TRUE)

            if (!dir.exists(libPath)) {
              stop("libPath directory creation failed.\n",
                   "Try creating on each machine manually, using e.g.,\n",
                   "  mkdir -p ", libPath)
            }
          }

          if (!"Require" %in% rownames(installed.packages())) {
            remotes::install_github("PredictiveEcology/Require@development")
          } else if (packageVersion("Require") < "0.1.0.9003") {
            remotes::install_github("PredictiveEcology/Require@development")
          }

          logPath <- Require::checkPath(dirname(logPath), create = TRUE)

          message(Sys.info()[["nodename"]])

          # Use Require with minimum version number as the mechanism for updating; remotes is
          #    too crazy with installing same package multiple times as recursive packages
          #    are dealt with
          Require::Require("PredictiveEcology/SpaDES.install@development")
          SpaDES.install::installSourcePackages() ## should be "rerun" proof, i.e., won't reinstall

          # This will install the versions of SpaDES.tools and fireSenseUtils that are on the main machine
          Require::Require(
            c(
              "dqrng",
              paste0("PredictiveEcology/SpaDES.tools@development (>=", packageVersionST, ")"),
              paste0("PredictiveEcology/fireSenseUtils@development (>=", packageVersionFSU, ")")
            ),
            upgrade = FALSE
          )
        }
      )
      parallel::stopCluster(cl)
    }

    ## Now make full cluster with one worker per core listed in "cores"
    message("Starting ", paste(paste(names(table(cores))), "x", table(cores),
      collapse = ", "
    ), " clusters")
    message("Starting main parallel cluster ...")
    st <- system.time({
      cl <- parallelly::makeClusterPSOCK(cores,
        revtunnel = revtunnel,
        outfile = logPath, rscript_libs = libPath
      )
    })

    on.exit(stopCluster(cl))
    message(
      "it took ", round(st[3], 2), "s to start ",
      paste(paste(names(table(cores))), "x", table(cores), collapse = ", "), " threads"
    )
    message("Moving objects to each node in cluster")

    stMoveObjects <- try({
      system.time({
        objsToCopy <- mget(unlist(objsNeeded))
        filenameForTransfer <- tempfile(fileext = ".qs")
        Require::checkPath(dirname(filenameForTransfer), create = TRUE) # during development, this was deleted accidentally
        qs::qsave(objsToCopy, file = filenameForTransfer)
        stExport <- system.time({
          outExp <- clusterExport(cl, varlist = "filenameForTransfer", envir = environment())
        })
        out11 <- clusterEvalQ(cl, {
          Require::checkPath(dirname(filenameForTransfer), create = TRUE)
        })
        out <- lapply(setdiff(unique(cores), "localhost"), function(ip) {
          st1 <- system.time(system(paste0("rsync -a ", filenameForTransfer, " ", ip, ":", filenameForTransfer)))
        })
        out <- clusterEvalQ(cl, {
          out <- qs::qread(file = filenameForTransfer)
          list2env(out, envir = .GlobalEnv)
        })
        # Delete the file
        out <- clusterEvalQ(cl, {
          if (dir.exists(dirname(filenameForTransfer))) {
            try(unlink(dirname(filenameForTransfer), recursive = TRUE), silent = TRUE)
          }
        })
      })
    })

    if (is(stMoveObjects, "try-error")) {
      message("The attempt to move objects to cluster using rsync and qs failed; trying clusterExport")
      stMoveObjects <- system.time(clusterExport(cl, objsNeeded, envir = environment()))
      list2env(mget(unlist(objsNeeded), envir = environment()), envir = .GlobalEnv)
    }
    message("it took ", round(stMoveObjects[3], 2), "s to move objects to nodes")
    message("loading packages in cluster nodes")
    stPackages <- system.time(parallel::clusterEvalQ(
      cl,
      {
        for (i in c(
          "kSamples", "magrittr", "raster", "data.table",
          "SpaDES.tools", "fireSenseUtils"
        )) {
          library(i, character.only = TRUE)
        }
        message("loading ", i, " at ", Sys.time())
      }
    ))
    message("it took ", round(stPackages[3], 2), "s to load packages")

    control$cluster <- cl
  } else {
    list2env(mget(unlist(objsNeeded), envir = environment()), envir = .GlobalEnv)
  }

  #####################################################################
  # DEOptim call
  #####################################################################
  DE <- # Cache(
    DEoptimIterative(
      itermax = itermax, lower = lower,
      upper = upper,
      control = do.call("DEoptim.control", control),
      FS_formula = FS_formula,
      covMinMax = covMinMax,
      # tests = c("mad", "SNLL_FS"),
      tests = tests,
      maxFireSpread = maxFireSpread,
      objFunCoresInternal = objFunCoresInternal,
      Nreps = Nreps,
      .verbose = .verbose,
      mutuallyExclusive = mutuallyExclusive,
      doObjFunAssertions = doObjFunAssertions,
      visualizeDEoptim = visualizeDEoptim,
      .plotSize = .plotSize,
      # cachePath = cachePath,
      iterStep = iterStep,
      thresh = thresh,
      # omitArgs = c("verbose")
    ) # ,
  # cacheId = "cd495b412420ad4a") # iteration 201 to 300
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

  aa <- as.data.table(t(DE$member$pop))
  aa[, pars := paste0("par", 1:NROW(aa))]
  dim1 <- floor(sqrt(NROW(aa)))
  dim2 <- NROW(aa) / dim1
  par(mfrow = c(dim1, dim2))
  aa[, hist(t(.SD), main = as.character(.BY))[[2]], by = pars]
}

#' @param control DESCRIPTION NEEDED
#' @template mutuallyExclusive
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
                             tests = c("SNLL", "adTest"),
                             objFunCoresInternal,
                             maxFireSpread,
                             Nreps,
                             visualizeDEoptim,
                             cachePath,
                             mutuallyExclusive,
                             doObjFunAssertions = getOption("fireSenseUtils.assertions", TRUE),
                             iterStep = 25,
                             thresh = 550,
                             .verbose,
                             .plotSize = list(height = 1600, width = 2000)) {
  data.table::setDTthreads(1)
  message("starting DEoptimIterative at ", Sys.time())
  x1 <- rnorm(1e2, 1, 2) # this is for debugging below
  DE <- list()
  for (iter in seq_len(ceiling(itermax / iterStep))) {
    control$itermax <- pmin(iterStep, itermax - iterStep * (iter - 1))
    control$storepopfrom <- control$itermax + 1

    controlArgs <- do.call("DEoptim.control", control)
    controlForCache <- controlArgs[c(
      "VTR", "strategy", "NP", "CR", "F", "bs", "trace",
      "initialpop", "p", "c", "reltol",
      "packages", "parVar", "foreachArgs"
    )]

    if (TRUE) {
      st1 <- system.time(DE[[iter]] <- # Cache(
        DEoptimForCache(
          fireSenseUtils::.objfunSpreadFit,
          lower = lower,
          upper = upper,
          control = controlArgs,
          FS_formula = FS_formula,
          covMinMax = covMinMax,
          tests = tests,
          maxFireSpread = maxFireSpread,
          mutuallyExclusive = mutuallyExclusive,
          doAssertions = doObjFunAssertions,
          Nreps = Nreps,
          controlForCache = controlForCache,
          objFunCoresInternal = objFunCoresInternal,
          thresh = thresh,
          verbose = .verbose # ,
          # omitArgs = c("verbose", "control")
        ))
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
    if (isTRUE(visualizeDEoptim)) {
      if (!isRstudioServer()) {
        png(
          filename = paste0("DE_pars", as.character(Sys.time()), "_", Sys.getpid(), ".png"),
          width = .plotSize$width, height = .plotSize$height
        )
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
