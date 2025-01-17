utils::globalVariables(c(
  "burned", "cell", "lcc", "pixelIndex", "speciesCode", "YEAR", "yearRange"
))

#' Calculate proportional burn of landcover and tree species
#'
#' @template pixelGroupMap
#'
#' @template cohortData
#'
#' @param rstLCC a landcover map
#'
#' @param nonflammableLCC nonflammable landcover in `rstLCC`
#'
#' @param fires a single `sf` or `SpatVector` object of fire polygons containing a `YEAR` column
#'
#' @param yearRange  the range of years represented by this landscape
#'
#' @return `data.table` with cell, biomass, tree species or lcc for non-forest, and year of fire
#'
#' @examples
#' # fuelClassPrep(
#' #   sim$pixelGroupMap2011, sim$cohortData2011, sim$rstLCC2011,
#' #   nonflammableLCC = P(sim)$nonflammableLCC,
#' #   fires = TODO,
#' #   yearRange = c(2012, 2020)
#' # )
#'
#' @export
#' @importFrom data.table as.data.table
#'
fuelClassPrep <- function(pixelGroupMap, cohortData, rstLCC, nonflammableLCC, fires, yearRange) {
  ## Filter fires by year range
  firesSubset <- fires[fires$YEAR >= yearRange[1] & fires$YEAR <= yearRange[2], ] |>
    rasterize(pixelGroupMap, field = "YEAR", fun = "min") |>
    as.data.frame(cells = TRUE) |> as.data.table()

  ## Add pixels to cohort data
  speciesData <- addPixels2CohortData(cohortData = cohortData, pixelGroupMap)

  ## Prepare landscape data excluding non-flammable land cover classes
  landscapeData <- as.data.table(as.data.frame(rstLCC, cells = TRUE)) |>
    setnames(c("cell", "lcc"))
  landscapeData <- landscapeData[!lcc %in% nonflammableLCC]

  ## Merge species data with landscape data
  landscapeData <- speciesData[landscapeData, on = c("pixelIndex" = "cell")]

  ## Merge fires data with landscape data
  landscapeData <- firesSubset[landscapeData, on = c("cell" = "pixelIndex")]
  landscapeData[, burned := ifelse(is.na(YEAR), 0, 1)]
  landscapeData[is.na(speciesCode), speciesCode := paste0("LCC", lcc)]
  landscapeData[, year := yearRange[1]]

  return(landscapeData)
}
