#' Threads for training multi-threaded libraries
#'
#' The thread budget for work that should be multi-threaded, i.e. training `xgboost` models.
#' (`xgboost` otherwise uses every core: restoring one set of saved boosters took 110 s of CPU in
#' 2.3 s on an 80-core host, with many R workers at once. See [setXgbThreads()] for the default.)
#'
#' The count is, in order:
#' 1. `getOption("fireSense.nThreads")` if set (e.g. `options(fireSense.nThreads = 4)`);
#' 2. otherwise `min(8, max(2, floor(cores - load)))`, where `load` is the host's 1-minute load
#'    average (`/proc/loadavg`). Training therefore gets 2 to 8 threads, based on the free cores.
#'    Where the load cannot be read (not Linux), it is taken as 0, i.e. `min(8, cores)`.
#'
#' @param cores Number of cores on the host.
#' @return Integer, at least 1.
#' @export
nThreads <- function(cores = parallel::detectCores()) {
  n <- getOption("fireSense.nThreads")
  if (!is.null(n)) return(max(1L, as.integer(n[1])))
  cores <- if (is.na(cores)) 1L else cores
  load <- .loadAvg1()
  if (is.na(load)) load <- 0
  as.integer(min(8, max(2, floor(cores - load))))
}

## 1-minute load average of this host; NA where /proc/loadavg cannot be read
.loadAvg1 <- function() {
  x <- tryCatch(suppressWarnings(readLines("/proc/loadavg", n = 1L)), error = function(e) NA_character_)
  if (length(x) != 1L) return(NA_real_)
  suppressWarnings(as.numeric(strsplit(x, " ", fixed = TRUE)[[1]][1]))
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
