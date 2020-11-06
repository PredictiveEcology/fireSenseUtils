#' Create PCA data for glm - for fitting and predicting
#'
#' @param dataForPCA a data.table of fireSense covariates for each pixelID,
#'  with additional cols pixelGroup and year if running for FS fit
#' @param PCAmodel if predicting, the PCA model from the fitting stage
#'
#' @return a list of the PCA and a data.table with components converted to integer
#'
#' @importFrom data.table data.table
#' @importFrom LandR asInteger
#'
#' @export
#' @rdname makeVegTerrainPCA
makeVegTerrainPCA <- function(dataForPCA, PCAmodel = NULL) {

  if (is.null(PCAmodel)) {
    #year should be present in data
    vegTerrainPCA <- prcomp(dataForPCA[, .SD, .SDcols = !c("pixelGroup", 'pixelID', 'year')],
                            center = TRUE, scale. = TRUE, rank = 10)
  } else{
    #To be tested
    vegTerrainPCA <- prcomp(object = PCAmodel, newdata = dataForPCA[, .SD., .SDcols = !c("pixelGroup", 'pixelID')])
  }
  # store as Integer
  vegComponents <- as.data.table(vegTerrainPCA$x * 1000)
  vegComponents <- vegComponents[, lapply(.SD, asInteger), .SDcols = colnames(vegComponents)]
  vegComponents[, pixelID := dataForPCA$pixelID]

  if (!is.null(dataForPCA$year)) {
    vegComponents[, year := dataForPCA$year]
  }

  return(list('vegComponents' = vegComponents, 'vegTerrainPCA' = vegTerrainPCA))
  }
