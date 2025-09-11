#' Get current state of FireSense Fit parameters as a map
#'
#' Download the FireSense parameter object, strip the list columns and
#' convert to a SpatVector or SpatRaster. If `plot = TRUE`, this will also
#' download a map that represents the forested lands of canada, to be plotted
#' on the plot device with the fireSense parameter ELFs that currently have
#' estimated parameters.
#'
#' @param rasterize Logical. If `TRUE`, then the function returns the rasterized
#'   map. `FALSE`, the default, is a `SpatVector`.
#' @param plot Logical. If `TRUE`, the default, then the map will be plotted to
#'   device.
#' @param res The resolution of the raster if `rasterize = TRUE`
#' @param ... Other parameters passed to `fireSenseCloudParameters`
#' @return The map, either `SpatRaster` or `SpatVector`
#' @export
#' @seealso [fireSenseCloudParameters()]
fireSenseCloudParametersMap <-
  function(rasterize = FALSE, plot = TRUE,
           res = 5000, ...) {
    oo <- fireSenseCloudParameters(...)
    if (isTRUE(plot)) {
      can <- scfmutils::prepInputsFireRegimePolys(type = "FRU")
      templ <- terra::ext(can) |>
        terra::rast(res = res, crs = can, vals = 1)
      canRas <- terra::rasterize(can, templ)
    }

    if (isTRUE(rasterize)) {
      b <- terra::rasterize(terra::vect(sf::st_as_sf(oo[, 5:6])), canRas)
      canRas[b>0] <- 2
      if (plot %in% TRUE)
        terra::plot(canRas, main = "FireSense Fit completed")
      oo <- canRas
    } else {
      oo <- sf::st_as_sf(oo[, 5:6])
      if (plot %in% TRUE) {
        oo <- terra::project(terra::vect(oo), can)
        # canV <- terra::vect(can)
        # canV <- terra::union(canV)
        canV2 <- terra::as.polygons(canRas)
        terra::plot(canV2, col = "transparent", main = "FireSense Fit completed")
        terra::plot(oo, add = TRUE, col = "yellow")
      }

    }
    oo
  }



#' Get current state of FireSense Fit parameters as a map
#'
#' Download the FireSense parameter object, strip the list columns and
#' convert to a SpatVector or SpatRaster.
#'
#' @param url Http url of the fireSense object with parameters. Default is
#'   correct on GoogleDrive. Must be a folder.
#' @param targetFile A filename to search on the folder on Google Drive. The
#'   default is correct for fireSense parameters.
#' @param destinationPath A path for the downloaded file.
#' @return The object and the `rds` file saved to destinationPath.
#' @export
#' @seealso [fireSenseCloudParametersMap()]
fireSenseCloudParameters <- function(url = "https://drive.google.com/file/d/1xQGAhBCRimYQC_GWA0lOctZU4VSlyojS/view?usp=drivesdk",
                                     targetFile = "fireSenseParams.rds",
                                     destinationPath = ".", useCache = TRUE) {
  prepInputs(targetFile = targetFile,
                   url = url,
                   destinationPath = destinationPath,
                   useCache = useCache, purge = 7, overwrite = TRUE)
}
