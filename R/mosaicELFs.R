#' Blend per-ELF rasters into one national raster
#'
#' Each ELF's raster covers the ELF's core and buffer, so neighbours overlap. Each ELF is given a weight that grows
#' from 0 at the outer edge of its footprint (core plus buffer) to the distance inward in metres, and each cell takes
#' the weighted mean of the ELFs that have a value there, so values change smoothly across the borders. The
#' footprint is the ELF's own map (`rasWhole` above 0), not the cells of its raster that have a value: holes in the
#' raster (non-burnable cells) make no edge. A cell with a value from one ELF takes it. A cell with values
#' only where the weights are 0 (the outermost cells of the ELF's footprint, with no neighbour) takes their mean. A cell
#' with no value stays `NA`. Each raster is projected to the national grid with `method`, so the ELFs may be in different
#' projections.
#'
#' @param files A named list with one element per product (e.g. `simulated`, `observed`): a named character vector
#'   of raster files, one per ELF, named by the ELF. Every product is blended with the same weights.
#' @param ELFs A list with `rasWhole`, the ELF maps (one layer per ELF on one grid; 0 is outside the ELF, 1 buffer, 2
#'   core), as [makeELFs()] and [mergeELFs()] give them. ELFs of `files` that it lacks are an error.
#' @param filename A character vector named as `files`: where to write each product. The output is a Float32
#'   cloud-optimized GeoTIFF (`gdal`), `NA` as no-data.
#' @param res,crs Cell size (m) and CRS of the national grid. The default CRS is the Canada Atlas Lambert
#'   (EPSG:3978), whose grid starts at multiples of `res`.
#' @param method Projection method for continuous values; see [terra::project()].
#' @param gdal Options of the COG driver.
#'
#' @return `filename`, invisibly.
#' @export
#' @importFrom terra classify crop distance ext extend ifel mosaic project rast sprc trim writeRaster
mosaicELFs <- function(files, ELFs, filename, res = 250, crs = "EPSG:3978", method = "bilinear",
                       gdal = c("COMPRESS=DEFLATE", "PREDICTOR=3", "OVERVIEW_RESAMPLING=AVERAGE")) {
  stopifnot(is.list(files), setequal(names(files), names(filename)))
  ids <- unique(unlist(lapply(files, names)))
  missingELFs <- setdiff(ids, names(ELFs$rasWhole))
  if (length(missingELFs)) stop("mosaicELFs(): no ELF map for ", paste(missingELFs, collapse = ", "))

  ## the footprint of each ELF, in the national CRS: trimmed to the ELF, with a ring of outside around it
  foot <- lapply(stats::setNames(ids, ids), function(id) {
    f <- terra::extend(terra::trim(terra::classify(ELFs$rasWhole[[id]], cbind(0, NA))), 2)
    terra::ifel(is.na(f), 0, 1)
  })
  grid <- .nationalGrid(foot, crs, res)
  ## weight inward from the footprint's outer edge, on the national grid (a cell outside has 0)
  weight <- lapply(foot, function(f) {
    g <- terra::crop(grid, terra::ext(terra::project(terra::as.polygons(terra::ext(f), crs = terra::crs(f)), crs)), snap = "out")
    inside <- terra::project(f, g, method = "near")
    terra::distance(terra::ifel(inside > 0, NA, 1))
  })

  for (product in names(files)) {
    parts <- lapply(stats::setNames(names(files[[product]]), names(files[[product]])), function(id) {
      w <- weight[[id]]
      v <- terra::project(terra::rast(files[[product]][[id]]), w, method = method)
      valued <- !is.na(v)
      list(num = terra::ifel(valued, v * w, 0), den = terra::ifel(valued, w, 0),
           sum = terra::ifel(valued, v, 0), n = terra::ifel(valued, 1, 0))
    })
    total <- lapply(c("num", "den", "sum", "n"), function(k)
      terra::mosaic(terra::sprc(lapply(parts, `[[`, k)), fun = "sum"))
    names(total) <- c("num", "den", "sum", "n")
    out <- terra::ifel(total$n == 0, NA, terra::ifel(total$den > 0, total$num / total$den, total$sum / total$n))
    names(out) <- product
    terra::writeRaster(out, filename[[product]], filetype = "COG", datatype = "FLT4S", gdal = gdal, overwrite = TRUE)
  }
  invisible(filename)
}

## The national grid of `res` cells (cell edges at multiples of `res`) covering the footprints `foot` (rasters in any CRS)
.nationalGrid <- function(foot, crs, res) {
  exts <- vapply(foot, function(f)
    as.vector(terra::ext(terra::project(terra::as.polygons(terra::ext(f), crs = terra::crs(f)), crs))), numeric(4))
  e <- c(floor(min(exts[1, ]) / res) * res, ceiling(max(exts[2, ]) / res) * res,
         floor(min(exts[3, ]) / res) * res, ceiling(max(exts[4, ]) / res) * res)
  terra::rast(terra::ext(e), resolution = res, crs = crs)
}
