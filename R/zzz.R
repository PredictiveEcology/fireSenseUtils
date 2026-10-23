.onAttach <- function(libname, pkgname) {
  packageStartupMessage("Using fireSenseUtils version ", utils::packageVersion("fireSenseUtils"))
  invisible()
}

.onLoad <- function(libname, pkgname) {
  if (requireNamespace("xgboost", quietly = TRUE))
    try(setXgbThreads(), silent = TRUE) # 1, or getOption("fireSense.nThreads")
  invisible()
}
