#' Threads this R session may use for multi-threaded libraries
#'
#' The thread budget for work that should be multi-threaded, i.e. training `xgboost` models.
#' (`xgboost` otherwise uses every core: restoring one set of saved boosters took 110 s of CPU in
#' 2.3 s on an 80-core host, with many R workers at once. See [setXgbThreads()] for the default.)
#'
#' The count is, in order:
#' 1. `getOption("fireSense.nThreads")` if set (e.g. `options(fireSense.nThreads = 4)`);
#' 2. otherwise `max(1, floor(cores / R processes on this host))`, where the processes are those
#'    running R's `exec/R` binary (workers, cluster nodes, this session). A lone session gets
#'    every core; one of 26 workers on 80 cores gets 3. Where the processes cannot be counted
#'    (not Linux), the session is assumed to be alone.
#'
#' @param cores Number of cores on the host.
#' @return Integer, at least 1.
#' @export
nThreads <- function(cores = parallel::detectCores()) {
  n <- getOption("fireSense.nThreads")
  if (!is.null(n)) return(max(1L, as.integer(n[1])))
  cores <- if (is.na(cores)) 1L else cores
  max(1L, as.integer(floor(cores / nRProcesses())))
}

nRProcesses <- function() {
  ps <- tryCatch(suppressWarnings(system2("ps", c("-eo", "args="), stdout = TRUE, stderr = FALSE)),
                 error = function(e) character())
  max(1L, sum(grepl("exec/R( |$)", ps)))
}

#' Set the number of threads `xgboost` uses in this session
#'
#' Sets `xgboost`'s global `nthread` (`xgboost::xgb.set.config()`). It covers restoring saved
#' boosters (`readRDS()`, a cache restore) and `predict()`, which do not take their own
#' `nthread`. The default is **1**, applied when fireSenseUtils is loaded (a booster is restored
#' before any module code runs when a cached event is read back). To change it, set
#' `options(fireSense.nThreads = N)` before loading the package, or call `setXgbThreads(N)`.
#' Code that trains boosters passes `nthread = nThreads()` itself. Other packages that use
#' `xgboost` should set their own `nthread`; this does not manage them.
#'
#' @param n Number of threads; default `getOption("fireSense.nThreads", 1)`.
#' @return `n`, invisibly.
#' @export
setXgbThreads <- function(n = getOption("fireSense.nThreads", 1L)) {
  xgboost::xgb.set.config(nthread = max(1L, as.integer(n[1])))
  invisible(n)
}
