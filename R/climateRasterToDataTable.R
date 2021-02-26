#' Converts stacks of climate rasters to data.table and optionally subsets to index
#'
#' @param historicalClimateRasters named list of raster stack(s)
#' @param Index optional list of `data.table`s named by `fireYear` containing fire buffer indices
#'
#' @return a `data.table` of variables for fireSense PCA
#'
#' @export
#' @importFrom data.table data.table setnames
#' @importFrom LandR asInteger
#' @rdname climateRasterToDataTable
climateRasterToDataTable <- function(historicalClimateRasters, Index = NULL) {
  climatePCAdat <- lapply(names(historicalClimateRasters),
                          FUN = function(var, stackList = historicalClimateRasters) {
                            annualStack <- stackList[[var]]
                            annualVars <- lapply(1:nlayers(annualStack), FUN = function(layer) {
                              layer <- annualStack[[layer]]
                              rasterDat <- data.table(pixelID = 1:ncell(layer),
                                                      value = asInteger(getValues(layer)),
                                                      year = names(layer))
                            })
                            annualVars <- rbindlist(annualVars)
                            annualVars <- na.omit(annualVars)

                            ## Index will be flammable cells only - should always be passed
                            if (!is.null(Index)) {
                              annualVars <- annualVars[pixelID %in% Index,]
                            }
                            setnames(annualVars, old = "value", new = var)
                            return(annualVars)
                          })
  return(climatePCAdat)
}
