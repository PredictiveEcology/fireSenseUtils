globalVariables(c(
  ".SD", ".SDcols", "..dontWant", "year", "flammable"
))

#' Create landcoverDT object to classify and track non-forest lcc
#' @param rstLCC landcover raster
#' @template flammableRTM
#' @param forestedLCC vector of values representing forested landcover classes in `rstLCC`
#' @param nonForestedLCCGroups a named list of non-forested flammable landcover groups
#' @return a data.table with columns for pixelID and binary presence of landcover
#'
#' @export
#' @importFrom data.table data.table set setDT
#' @importFrom fastDummies dummy_cols
#' @importFrom sf %>%
#' @importFrom terra values ncell
#'
#' @rdname makeLandcoverDT
makeLandcoverDT <- function(rstLCC, flammableRTM, forestedLCC, nonForestedLCCGroups) {
  lcc <- data.table(
    pixelID = 1:ncell(rstLCC),
    lcc = values(rstLCC, mat = FALSE),
    flammable = values(flammableRTM, mat = FALSE)
  ) %>%
    .[lcc != 0, ] %>%
    .[!is.na(flammable), ] %>%
    .[flammable == 1, ] %>%
    .[, lcc := as.factor(lcc)]
  set(lcc, NULL, "flammable", NULL)
  lcc <- dummy_cols(
    .data = lcc,
    select_columns = "lcc",
    remove_selected_columns = TRUE,
    remove_first_dummy = FALSE, # no need if all forest LC is removed anyway
    ignore_na = TRUE
  )
  #
  forestedLCC <- paste0("lcc_", forestedLCC)
  forestedLCC <- colnames(lcc)[colnames(lcc) %in% forestedLCC]
  set(lcc, NULL, j = forestedLCC, NULL) # removes columns of forested LCC

  # Group dummy columns into similar landcovers
  lccColsPreGroup <- colnames(lcc)[!colnames(lcc) %in% c("pixelID")]
  setDT(lcc) # pre-allocate space for new columns

  for (i in names(nonForestedLCCGroups)) {
    classes <- paste0("lcc_", nonForestedLCCGroups[[i]])
    classes <- classes[classes %in% colnames(lcc)]
    set(lcc, NULL, eval(i), rowSums(lcc[, .SD, .SDcols = classes]))
  }
  rm(classes)
  set(lcc, NULL, lccColsPreGroup, NULL)
  return(lcc)
}
