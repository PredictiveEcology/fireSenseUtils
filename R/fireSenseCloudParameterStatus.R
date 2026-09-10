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
        terra::plot(canV2, col = "transparent", main = "FireSense Fit completed; ELFs with labels")
        if ("destinationPath" %in% ...names())
          destinationPath <- list(...)$destinationPath
        else
          destinationPath <- "."
        ELFs <- makeELFs(destinationPath = destinationPath, singleSpatVector = TRUE)
        ELFs2 <- makeELFs(destinationPath = destinationPath, singleSpatVector = FALSE)

        ## makeELFs() returns a list (rasters + `poly`); the single SpatVector is the `poly` element.
        ## cf. ELFsInStudyArea(), which already does this. `singleSpatVector` no longer switches the
        ## return type -- the code that used it is commented out in makeELFs().
        ELFsPoly <- ELFs$poly

        # terra::plot(ELFsPoly, col = c("turquoise", "yellow")[ELFsPoly$buffer], alpha = 0.5)
        terra::plot(ELFsPoly)
        terra::plot(oo, add = TRUE, col = rgb(255, 255, 0, alpha = 127, maxColorValue = 255))
        centroids_points <- terra::centroids(oo)
        keep <- ELFsPoly$buffer %in% 2 & !ELFsPoly$ID %in% centroids_points[[polygonIDTxt]]
        terra::text(terra::centroids(ELFsPoly[ keep, ]), labels = ELFsPoly$ID[keep], cex = 0.7)
        terra::text(centroids_points, labels = centroids_points[[polygonIDTxt]], cex = 0.8, col = "blue")
      }

    }
    oo
  }

#' @export
#' @rdname makeELFs
#' @seealso [fireSenseCloudParameters()]
plotELFs <- function(destinationPath = ".") {
  ELFs <- makeELFs(destinationPath = destinationPath, singleSpatVector = TRUE) |>
    Cache()
  ELFs2 <- makeELFs(destinationPath = destinationPath, singleSpatVector = FALSE)|>
    Cache()

  ## makeELFs() returns a list (rasters + `poly`); the single SpatVector is the `poly` element.
  ## cf. ELFsInStudyArea(), which already does this. `singleSpatVector` no longer switches the
  ## return type -- the code that used it is commented out in makeELFs().
  ELFsPoly <- ELFs$poly

  # terra::plot(ELFsPoly, col = c("turquoise", "yellow")[ELFsPoly$buffer], alpha = 0.5)
  terra::plot(ELFsPoly)
  keep <- ELFsPoly$buffer %in% 2
  terra::text(terra::centroids(ELFsPoly[ keep, ]), labels = ELFsPoly$ID[keep], cex = 0.7)
  return(invisible(ELFs))
}


#' Get current state of FireSense Fit parameters as a map
#'
#' Download the FireSense parameter object, strip the list columns and
#' convert to a SpatVector or SpatRaster.
#'
#' The file is downloaded from Google Drive on every call, so the result is always
#' the current state of the shared parameter file, never a copy left on disk by an
#' earlier call.
#'
#' @param url Http url of the fireSense object with parameters on Google Drive:
#'   either the file itself (the default) or the folder that contains `targetFile`.
#' @param targetFile The name of the file to find when `url` is a folder. The
#'   default is correct for fireSense parameters.
#' @param destinationPath A path where a copy of the downloaded file is saved.
#' @param useCache Ignored. Kept so existing calls still work: the file is always
#'   downloaded fresh.
#' @return The object and the `rds` file saved to destinationPath.
#' @export
#' @seealso [fireSenseCloudParametersMap()]
fireSenseCloudParameters <- function(
    url = paste0(
      "https://drive.google.com/file/d/",
      "1oJX7V5wPSj49C6wt59yD9RBplfZ7Cp-f/view?usp=drive_link"
    ),
    targetFile = "fireSenseParams.rds",
    destinationPath = ".", useCache = TRUE) {

  # KNN = "https://drive.google.com/file/d/1xQGAhBCRimYQC_GWA0lOctZU4VSlyojS/view?usp=drivesdk"
  ## Downloaded directly rather than through prepInputs(): prepInputs() returns the copy
  ## already on disk (or linked from destinationPathShared) whenever it matches
  ## CHECKSUMS.txt, and neither `purge` (rebuilds CHECKSUMS.txt entries) nor `overwrite`
  ## (the written output) makes it download again -- so a changed file on Drive was
  ## never seen. The file is small, so a fresh download each call is cheap.
  remote <- googledrive::drive_get(googledrive::as_id(url))
  if (isTRUE(googledrive::is_folder(remote))) {
    inFolder <- googledrive::drive_ls(remote)
    remote <- inFolder[inFolder$name %in% targetFile, ]
    if (NROW(remote) != 1L)
      stop("Expected one file named '", targetFile, "' in ", url, "; found ", NROW(remote), ".")
  }
  tmp <- tempfile(fileext = ".rds")
  on.exit(unlink(tmp), add = TRUE)
  googledrive::drive_download(remote, path = tmp, overwrite = TRUE)
  out <- readRDS(tmp)

  dir.create(destinationPath, showWarnings = FALSE, recursive = TRUE)
  dest <- file.path(destinationPath, targetFile)
  unlink(dest) # a new file, not a write through a hard link into a shared copy
  file.copy(tmp, dest)
  out
}
