globalVariables(c(
  ".SD", ".SDcols", "..dontWant", "year", "flammable"
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
  if (is.null(PCAmodel)) {
    vegComponents <- as.data.table(vegTerrainPCA$x)
  } else {
    vegComponents <- as.data.table(vegTerrainPCA) #object is a matrix, not prcomp
    vegTerrainPCA <- NULL # predict only needs the data.table
  }
    #put age and pixelID back in
  set(vegComponents, NULL, "pixelID" ,dataForPCA$pixelID)
  set(vegComponents, NULL, "youngAge", dataForPCA$youngAge)

  #put year in if present (i.e., during fitting)
  if (!is.null(dataForPCA$year)) {
    vegComponents[, year := dataForPCA$year]
  }
  vegTerrainPCA$x <- NULL ## remove data from PCA object

  #put age back in
  return(list("vegComponents" = vegComponents, "vegTerrainPCA" = vegTerrainPCA))
}

#' Create landcoverDT object for vegetation and terrain PCA
#'
#' @param rstLCC landcover raster
#' @template flammableRTM
#' @param forestedLCC vector of values representing forested landcover classes in \code{rstLCC}
#' @param nonForestedLCCGroups a named list of non-forested flammable landcover groups
#' @return a data.table with columns for pixelID and binary presence of landcover
#'
#' @export
#' @importFrom data.table data.table set setDT
#' @importFrom raster getValues ncell
#' @importFrom fastDummies dummy_cols
#' @importFrom magrittr %>%
#'
#' @rdname makeLandcoverDT
makeLandcoverDT <- function(rstLCC, flammableRTM, forestedLCC, nonForestedLCCGroups) {

  lcc <- data.table(pixelID = 1:ncell(rstLCC),
                    lcc = getValues(rstLCC),
                    flammable = getValues(flammableRTM)) %>%
    .[lcc != 0,] %>%
    .[!is.na(flammable),] %>%
    .[flammable == 1,] %>%
    .[, lcc := as.factor(lcc)]
  set(lcc, NULL, "flammable", NULL)
  lcc <- dummy_cols(.data = lcc,
                    select_columns = "lcc",
                    remove_selected_columns = TRUE,
                    remove_first_dummy = FALSE, #no need if all forest LC is removed anyway
                    ignore_na = TRUE)
  #
  forestedLCC <- paste0("lcc_", forestedLCC)
  forestedLCC <- colnames(lcc)[colnames(lcc) %in% forestedLCC]
  set(lcc, NULL, j = forestedLCC, NULL) #removes columns of forested LCC

  # Group dummy columns into similar landcovers
  lccColsPreGroup <- colnames(lcc)[!colnames(lcc) %in% c("pixelID")]
  setDT(lcc) #pre-allocate space for new columns

  for (i in names(nonForestedLCCGroups)) {
    classes <- paste0("lcc_", nonForestedLCCGroups[[i]])
    classes <- classes[classes %in% colnames(lcc)]
    set(lcc, NULL,  eval(i), rowSums(lcc[, .SD, .SDcols = classes]))
  }
  rm(classes)
  set(lcc, NULL, lccColsPreGroup, NULL)
  return(lcc)
}
