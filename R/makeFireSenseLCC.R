#' Prepare NTEMS Land Cover Data for Fire Simulation Models
#'
#' Processes NTEMS Land Cover Classification (LCC) data for a specific year,
#' adjusting non-flammable classes based on a flammability threshold during
#' reprojection to a target raster's resolution and extent.
#'
#' @details
#' This function performs several steps:
#' 1. Downloads and prepares NTEMS LCC data for the specified `neededYear` using
#'    `LandR::prepInputs_NTEMS_LCC_FAO`, cropping and masking to the target extent.
#' 2. Identifies the dominant *flammable* land cover class within each pixel of
#'    `to`. This is done by temporarily masking non-flammable classes
#'    (defined in `nonflammableLCC`) to NA in the source raster and then
#'    reprojecting using the 'mode' method.
#' 3. Calculates the proportion of *flammable* cover within each pixel of
#'    `to` by reprojecting a binary (1=flammable, 0=non-flammable)
#'    version of the source raster using the 'average' method.
#' 4. Creates the final raster: Pixels where the proportion of flammable cover
#'    (calculated in step 3) is below the `flammabilityThreshold` are assigned
#'    a value of 0 (representing non-flammable). Otherwise, pixels retain the
#'    dominant flammable LCC code identified in step 2.
#'
#' @param neededYear Numeric. The specific year required for the NTEMS Land
#'   Cover Classification (LCC) data.
#' @param nonflammableLCC Numeric vector. LCC codes representing non-flammable
#'   land cover types (e.g., water, rock, urban). Defaults to `c(20, 31, 32, 33)`.
#'   Pixels with these values are initially masked or converted for thresholding.
#' @param flammabilityThreshold Numeric. A value between 0 and 1. Target pixels
#'   where the proportion of underlying *flammable* source pixels is below this
#'   threshold will be reclassified to 0 (non-flammable) in the output raster.
#'   Defaults to 0.1.
#' @param writeTo Character string or `NULL`. Optional filename for the output
#'   LCC raster, relative to `destinationPath`. If `NULL` (default), the raster is
#'   not written to disk by this function.
#' @param destinationPath Character string. Directory path where source NTEMS data
#'   will be downloaded/loaded (via `prepInputs_NTEMS_LCC_FAO`) and where the
#'   output raster will be written if `writeTo` is specified.
#'
#' @return A list of two  `SpatRaster` objects:
#'   1) the processed land cover classification at the resolution, extent,
#'   and CRS of `to` with  values representing either the modal flammable
#'   LCC code or 0 if the pixel is below the `flammabilityThreshold`, and
#'   2) the proportion of each pixel that is flammable in the native resolution
#'
#' @importFrom terra project values `values<-` writeRaster classify
#' @importFrom LandR prepInputs_NTEMS_LCC_FAO
#' @importFrom reproducible postProcessTo
#' @importFrom reproducible .suffix
#' @export
#'
#' @examples
#' \dontrun{
#'   # Requires terra and potentially LandR installed
#'   # Need example to and maskTo objects
#'   library(terra)
#'
#'   # Create dummy to and maskTo
#'   ras <- rast(xmin = 0, xmax = 10, ymin = 0, ymax = 10, res = 1)
#'   values(ras) <- 1
#'   crs(ras) <- "EPSG:4326" # Example CRS
#'   sa <- ext(ras) |> as.polygons() |> vect()
#'   crs(sa) <- crs(ras)
#'
#'   # Define a destination path (replace with a real path)
#'   destPath <- tempdir()
#'
#'   # Run the function (will likely try to download data)
#'   # Note: prepInputs_NTEMS_LCC_FAO might require specific setup or data sources
#'   # This example might fail if data download/access isn't configured
#'   # fireLCC <- makeFireSenseLCC(neededYear = 2011, # Example year
#'   #                            to = ras,
#'   #                            maskTo = sa,
#'   #                            destinationPath = destPath)
#'   # plot(fireLCC)
#' }
makeFireSenseLCC <- function(neededYear, to, maskTo = NULL, # to, maskTo = NULL,
                             nonflammableLCC = c(20, 31, 32, 33),
                             flammabilityThreshold = 0.1, writeTo = NULL,
                             destinationPath) {
  # 1. Retrieve and prepare base NTEMS LCC data for the specified year
  #    - Crops to the extent of to
  #    - Masks to the maskTo polygon(s)
  message("Preparing base NTEMS LCC data...")

  maskToArg <- if (is.null(maskTo)) to else maskTo


  # if (exists("aaaa", envir = .GlobalEnv)) browser()
  rstLCC <- prepInputs_NTEMS_LCC_FAO(year = neededYear,
                                     disturbedCode = 240, # Optional: specify disturbed code if needed
                                     overwrite = TRUE,    # Consider parameterizing overwrite?
                                     destinationPath = destinationPath,
                                     cropTo = to,
                                     maskTo = maskToArg)
  # 2. Determine the dominant *flammable* LCC code at the target resolution
  #    - Mask non-flammable codes to NA in the source resolution LCC
  #    - Project to the target resolution using 'mode' aggregation. This finds
  #      the most frequent flammable LCC code among the source pixels
  #      contributing to each target pixel.
  message("Calculating dominant flammable LCC at target resolution...")
  allFlam <- terra::classify(rstLCC, matrix(c(nonflammableLCC, rep(NA, length(nonflammableLCC))), ncol = 2)) |>
    terra::project(y = to, method = "mode")

  # 3. Calculate the proportion of *flammable* cover at the target resolution
  #    - Create a binary raster at source resolution: 1 = flammable, 0 = non-flammable
  #    - Project this binary raster to the target resolution using 'average'
  #      aggregation. The result represents the proportion of flammable source
  #      pixels within each target pixel.
  message("Calculating proportion of flammable cover...")
  # Create binary raster: 0 where LCC is non-flammable, keep original LCC otherwise
  nonFlamMap <- terra::classify(rstLCC, matrix(c(nonflammableLCC, rep(0, length(nonflammableLCC))), ncol = 2))

  # Convert remaining LCC codes (flammable ones) to 1
  allVals <- freq(nonFlamMap)$value
  allVals <- setdiff(allVals, 0)
  # not memory safe:
  # nonFlamMap[nonFlamMap > 0] <- 1 # Now 1 = flammable, 0 = non-flammable
  if (length(allVals) > 1)
    nonFlamMap <- terra::subst(nonFlamMap, from = allVals, to = 1L)

  # Project using average to get proportion of flammable pixels
  flammableProp <- terra::project(nonFlamMap, to, method = "average")

  # 4. Apply flammability threshold
  #    - Where the proportion of flammable cover is less than the threshold,
  #      set the LCC value in the result raster ('allFlam') to 0 (non-flammable).
  message("Applying flammability threshold...")
  # not memory safe
  allFlam[flammableProp < flammabilityThreshold] <- 0

  # 5. Optionally write the result to disk
  if (!is.null(writeTo)) {
    outPath <- file.path(destinationPath, writeTo)
    propFlamPath <- .suffix(outPath, "_propFlam")
    message("Writing output raster to: ", outPath)
    # Ensure data type supports 0 if original LCC didn't (though usually integer, so fine)
    allFlam <- terra::writeRaster(allFlam, filename = outPath, overwrite = TRUE)
    flammableProp <- writeRaster(flammableProp, filename = propFlamPath, overwrite = TRUE)
  }

  output <- list("lcc" = allFlam, "flammableProp" = flammableProp)

  return(output)
}
