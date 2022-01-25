#' Update name of layers in a climate raster stack
#'
#' @param annualDataStack \code{RasterStack}
#' @param desiredYears character
#'
#' @export
updateStackYearNames <- function(annualDataStack, desiredYears) {
  objectName <- as.character(substitute(annualDataStack))

  grepTest4DigitYear <- "[[:digit:]]{4,4}$"
  namingConventionTxt <- paste(
    "does not have names that include the 4 digit year.",
    "Please use that naming convention."
  )

  # Sanity check -- both objects must have 4 length year in them
  hasYearStack <- all(grepl(grepTest4DigitYear, names(annualDataStack)))
  if (!isTRUE(hasYearStack)) {
    stop("The raster stack, ", objectName, ", ", namingConventionTxt)
  }

  hasYearParam <- all(grepl(grepTest4DigitYear, desiredYears))
  if (!isTRUE(hasYearParam)) {
    stop("The desiredYears argument", namingConventionTxt)
  }

  annualDataStackInt <- as.integer(gsub(".*([[:digit:]]{4,4}).*", "\\1", names(annualDataStack)))
  stopifnot(identical(annualDataStackInt, desiredYears))

  whichannualDataStackToKeep <- which(annualDataStackInt %in% desiredYears)
  annualDataStackInt <- annualDataStackInt[whichannualDataStackToKeep]
  namesWithYearPrepended <- paste0("year", annualDataStackInt)

  names(annualDataStack) <- namesWithYearPrepended
  annualDataStack
}
