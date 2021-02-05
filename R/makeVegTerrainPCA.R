globalVariables(c(
  ".SD", ".SDcols", "..dontWant", "year"
))

#' Create PCA data for glm - for fitting and predicting
#'
#' @param dataForPCA a data.table of fireSense covariates for each pixelID,
#'  with additional cols pixelGroup and year if running for FS fit
#' @param PCAmodel if predicting, the PCA model from the fitting stage
#' @param dontAlter The columns that shouldn't be transformed with log (x + 1)
#' @param dontWant The columns that shouldn't be part of the PCA
#'
#' @return a list of the PCA and a data.table with components converted to integer
#'
#' @export
#' @importFrom data.table data.table set setDT
#' @importFrom LandR asInteger
#' @importFrom stats prcomp
#' @rdname makeVegTerrainPCA
makeVegTerrainPCA <- function(dataForPCA, PCAmodel = NULL,
                              dontAlter = c("TPI", "HLI", "nonForest_highFlam", "nonForest_lowFlam"),
                              dontWant = c("pixelGroup", "pixelID", "year", "youngAge")) {

  # Log the species data
  doModify <- setdiff(colnames(dataForPCA), c(dontWant, dontAlter))
  dataForPCA[, (doModify) := lapply(.SD, function(x) log(x + 1)), .SDcols = doModify]
  dd <- dataForPCA[, !(..dontWant)]

  if (is.null(PCAmodel)) {
    #year should be present in data
    vegTerrainPCA <- prcomp(dd, center = TRUE, scale. = TRUE, rank = 10)
  } else{
    vegTerrainPCA <- predict(object = PCAmodel, newdata = dd)
  }
  # store as Integer
  if (is.null(PCAmodel)){
    vegComponents <- as.data.table(vegTerrainPCA$x)
  } else {
    vegComponents <- as.data.table(vegTerrainPCA) #object is a matrix, not prcomp
    vegTerrainPCA <- NULL # predict only needs the data.table
  }
    #put age and pixelID back in
  set(vegComponents, NULL, "pixelID" ,dataForPCA$pixelID)
  set(vegComponents, NULL, "youngAge", dataForPCA$youngAge)

  #put year in if present (ie during fitting)
  if (!is.null(dataForPCA$year)) {
    vegComponents[, year := dataForPCA$year]
  }
  vegTerrainPCA$x <- NULL

  #put age back in

  return(list("vegComponents" = vegComponents, "vegTerrainPCA" = vegTerrainPCA))
}
