.onAttach <- function(libname, pkgname) {
  packageStartupMessage("Using fireSenseUtils version ", utils::packageVersion("fireSenseUtils"))
  invisible()
}
