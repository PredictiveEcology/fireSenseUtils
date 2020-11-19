#' converts stacks of climate rasters to data.table and optionally subsets to index
#'
#' @param historicalClimateRasters named list of raster stack(s)
#' @param Index optional list of data.tables named by fireYear containing fire buffer indices
#' @return a data.table of variables for fireSense PCA
#' @importFrom data.table setnames data.table
#' @export
#' @rdname climateRasterToDataTable
climateRasterToDataTable <- function(historicalClimateRasters, Index = NULL) {
  climatePCAdat <- lapply(names(historicalClimateRasters),
                          FUN = function(var, stackList = historicalClimateRasters) {
                            annualStack <- stackList[[var]]
                            annualVars <- lapply(1:nlayers(annualStack), FUN = function(layer){
                              layer <- annualStack[[layer]]
                              rasterDat <- data.table(pixelID = 1:ncell(layer),
                                                      value = getValues(layer),
                                                      year = names(layer))

                            })
                            annualVars <- rbindlist(annualVars)
                                                        #most climate variables can be integer if coming from climateNA
                            if (!class(annualVars$value) == 'integer') {
                              warning(paste(var, "is not in integer form. It will be converted to integer.",
                                            "Please supply climate data with integer type if this is problematic"))
                              set(annualVars, NULL, 'value', as.integer(annualVars$value))
                            }
                            #Index will be flammable cells only - should always be passed
                            if (!is.null(Index)) {
                              annualVars <- annualVars[pixelID %in% Index,]
                            }
                            setnames(annualVars, old = 'value', new = var)
                            return(annualVars)
                          })
  return(climatePCAdat)
}
