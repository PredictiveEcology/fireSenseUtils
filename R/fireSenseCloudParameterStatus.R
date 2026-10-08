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

#' @param which Character. ELF names (e.g., `c("13.1", "4.1")`) whose cores `plotELFs()`
#'   fills with `fill`. Default `NULL` fills none.
#' @param fill The colour used to fill the `which` ELFs.
#' @param labels What `plotELFs()` prints on each ELF: `"code"` (the ELF name, e.g. `"13.1"`),
#'   `"name"` (the name of the ecozone the ELF lies in, e.g. `"Pacific Maritime"`; see
#'   [elfLabels()]) or `"none"`.
#' @param labelWhich `"all"` labels every ELF; `"highlighted"` labels only those in `which`.
#' @param buffers If `FALSE`, `plotELFs()` draws only the ELF cores, without the buffer rings.
#' @param axes `"m"` shows the map's own axes (metres); `"longlat"` keeps the map projection
#'   and draws a graticule labelled in degrees; `"none"` draws no axes and no box.
#' @export
#' @rdname makeELFs
#' @seealso [fireSenseCloudParameters()]
plotELFs <- function(destinationPath = ".", which = NULL, fill = "green",
                     labels = c("code", "name", "none"), labelWhich = c("all", "highlighted"),
                     buffers = TRUE, axes = c("m", "longlat", "none")) {
  labels <- match.arg(labels)
  labelWhich <- match.arg(labelWhich)
  axes <- match.arg(axes)
  ELFs <- makeELFs(destinationPath = destinationPath, singleSpatVector = TRUE) |>
    Cache()
  ELFs2 <- makeELFs(destinationPath = destinationPath, singleSpatVector = FALSE)|>
    Cache()

  ## makeELFs() returns a list (rasters + `poly`); the single SpatVector is the `poly` element.
  ## cf. ELFsInStudyArea(), which already does this. `singleSpatVector` no longer switches the
  ## return type -- the code that used it is commented out in makeELFs().
  ELFsPoly <- ELFs$poly
  drawn <- elfLayer(ELFsPoly, buffers = buffers)

  # terra::plot(ELFsPoly, col = c("turquoise", "yellow")[ELFsPoly$buffer], alpha = 0.5)
  ## a few labels on a map edge need room beside it to avoid the ELFs: widen the extent then
  args <- list(drawn)
  if (!identical(labels, "none") && identical(labelWhich, "highlighted")) {
    e <- terra::ext(drawn)
    args$ext <- e + 0.08 * (e$xmax - e$xmin) # terra pads each side by this amount
  }
  if (!identical(axes, "m")) args <- c(args, list(axes = FALSE, box = FALSE))
  if (identical(axes, "longlat")) args$mar <- c(3, 3, 1, 1)
  do.call(terra::plot, args)
  if (identical(axes, "longlat")) elfGraticule(drawn)
  keep <- ELFsPoly$buffer %in% 2
  if (length(which)) {
    missingELFs <- setdiff(which, ELFsPoly$ID)
    if (length(missingELFs))
      warning("ELFs not found: ", paste(missingELFs, collapse = ", "))
    terra::plot(ELFsPoly[keep & ELFsPoly$ID %in% which, ], col = fill, add = TRUE)
  }
  if (!identical(labels, "none")) {
    if (identical(labelWhich, "highlighted")) keep <- keep & ELFsPoly$ID %in% which
    if (any(keep)) {
      ## halo: a white outline keeps the text readable over the polygon lines
      cex <- if (identical(labelWhich, "highlighted")) 1 else 0.7 # few labels: print them larger
      cen <- terra::centroids(ELFsPoly[ keep, ])
      txt <- elfLabels(ELFsPoly$ID[keep], labels)
      xy <- terra::crds(cen)
      lab <- elfLabelXY(xy[, 1], xy[, 2], graphics::strwidth(txt, cex = cex),
                        graphics::strheight(txt, cex = cex), usr = graphics::par("usr"))
      moved <- lab$x != xy[, 1] | lab$y != xy[, 2]
      graphics::segments(xy[moved, 1], xy[moved, 2], lab$x[moved], lab$y[moved], col = "grey20", lwd = 1.5)
      ## halo: a white outline keeps the text readable over the polygon lines
      terra::text(terra::vect(cbind(lab$x, lab$y), crs = terra::crs(cen)), labels = txt,
                  cex = cex, halo = TRUE, hc = "white", hw = 0.2)
    }
  }
  return(invisible(ELFs))
}

## Ecozone (Canada's national ecological framework) names. An ELF code is
## "<ecozone>.<ecoprovince>[.<piece>]"; the ecoprovince shapefile that makeELFs() reads
## (ECOPROVINC) carries codes only, so the ecozone is the finest named level the ELF codes follow.
## Same names as the ZONE_NAME column of ecozone_shp.zip from sis.agr.gc.ca ("Boreal PLain" there).
ecozoneNames <- c(
  `1` = "Arctic Cordillera", `2` = "Northern Arctic", `3` = "Southern Arctic",
  `4` = "Taiga Plain", `5` = "Taiga Shield", `6` = "Boreal Shield",
  `7` = "Atlantic Maritime", `8` = "Mixedwood Plain", `9` = "Boreal Plain",
  `10` = "Prairie", `11` = "Taiga Cordillera", `12` = "Boreal Cordillera",
  `13` = "Pacific Maritime", `14` = "Montane Cordillera", `15` = "Hudson Plain"
)

#' Text to print on ELFs
#'
#' @param ids Character. ELF codes, e.g. `c("13.1", "6.2.1")`.
#' @param labels `"code"` returns `ids`; `"name"` returns the name of the ecozone each ELF
#'   lies in, with the code added in brackets where two of `ids` share a zone name, so that
#'   labels stay distinct; `"none"` returns `""`.
#' @return A character vector the length of `ids`.
#' @export
#' @examples
#' elfLabels(c("13.1", "5.4", "6.2.1", "6.2.2"), "name")
elfLabels <- function(ids, labels = c("code", "name", "none")) {
  labels <- match.arg(labels)
  ids <- as.character(ids)
  switch(labels,
         code = ids,
         none = rep("", length(ids)),
         name = {
           out <- unname(ecozoneNames[sub("\\..*$", "", ids)])
           out[is.na(out)] <- ids[is.na(out)]
           dup <- out %in% out[duplicated(out)]
           out[dup] <- paste0(out[dup], " (", ids[dup], ")")
           out
         })
}

## Where to centre each label so that boxes of width `w` and height `h` (user units) do not
## overlap one another and stay inside `usr` (the plot region, `c(xmin, xmax, ymin, ymax)`).
## A label stays on its point (`x`, `y`) unless that overlaps one already placed; it then moves
## to the nearest free spot in a ring around the point, so plotELFs() draws a leader line to it.
## A label with no free spot stays on its point.
elfLabelXY <- function(x, y, w, h, usr = c(-Inf, Inf, -Inf, Inf)) {
  ## boxes carry a margin of one line height, so neighbouring labels do not touch
  box <- function(cx, cy, i) c(cx - w[i] / 2 - h[i], cx + w[i] / 2 + h[i], cy - h[i], cy + h[i])
  free <- function(b, placed) {
    inside <- b[1] >= usr[1] && b[2] <= usr[2] && b[3] >= usr[3] && b[4] <= usr[4]
    inside && !any(vapply(placed, function(p) b[1] < p[2] && p[1] < b[2] && b[3] < p[4] && p[3] < b[4],
                          logical(1)))
  }
  angles <- seq(0, 2 * pi, length.out = 9)[-9]
  out <- list(x = x, y = y)
  placed <- list()
  for (i in seq_along(x)) {
    b <- box(x[i], y[i], i)
    if (!free(b, placed)) {
      for (d in c(3, 5, 7, 9) * h[i]) {
        cand <- lapply(angles, function(a) c(cos(a) * (w[i] / 2 + d), sin(a) * (h[i] / 2 + d)))
        ok <- Filter(function(o) free(box(x[i] + o[1], y[i] + o[2], i), placed), cand)
        if (length(ok)) {
          ## of the free spots at this distance, the one furthest from the other points
          away <- vapply(ok, function(o) min(c(Inf, sqrt((x[-i] - x[i] - o[1])^2 + (y[-i] - y[i] - o[2])^2))),
                         numeric(1))
          o <- ok[[which.max(away)]]
          out$x[i] <- x[i] + o[1]
          out$y[i] <- y[i] + o[2]
          b <- box(out$x[i], out$y[i], i)
          break
        }
      }
    }
    placed[[i]] <- b
  }
  out
}

## The polygons plotELFs() draws: with `buffers = FALSE`, only the cores (buffer == 2),
## as makeELFs() codes them (1 = buffer ring, 2 = ELF core).
elfLayer <- function(ELFsPoly, buffers = TRUE) {
  if (isTRUE(buffers)) ELFsPoly else ELFsPoly[ELFsPoly$buffer %in% 2, ]
}

## Graticule over the current map, in degrees, labelled where each line meets the
## left (parallels) or bottom (meridians) edge of `v`'s extent. The map keeps its projection.
elfGraticule <- function(v, step = 10) {
  e <- terra::ext(v)
  ll <- terra::ext(terra::project(terra::as.polygons(e, crs = terra::crs(v)), "EPSG:4326"))[]
  ## the multiples of `step` inside the range; none when the map spans less than one step
  grid <- function(lo, hi) {
    from <- ceiling(ll[[lo]] / step) * step
    to <- floor(ll[[hi]] / step) * step
    if (from > to) numeric(0) else seq(from, to, step)
  }
  ## one graticule line cut to the map extent
  line <- function(x, y) {
    terra::vect(cbind(x, y), "lines", crs = "EPSG:4326") |>
      terra::project(terra::crs(v)) |>
      terra::crop(e)
  }
  n <- 200
  lonRange <- seq(ll[["ymin"]] - 20, min(ll[["ymax"]] + 20, 89), length.out = n)
  latRange <- seq(ll[["xmin"]] - 20, ll[["xmax"]] + 20, length.out = n)
  ## a label sits at the line's end on the edge; lines that do not reach that edge get none
  edgeAt <- function(l, col, edgeValue, size) {
    xy <- terra::crds(l)
    if (!NROW(xy)) return(NA_real_)
    i <- which.min(xy[, col])
    if (xy[i, col] - edgeValue > 0.01 * size) NA_real_ else xy[i, 3 - col]
  }
  degLab <- function(d, pos, neg) paste0(abs(d), "\u00b0", if (d < 0) neg else pos)
  for (side in c(meridians = 1, parallels = 2)) {
    isLon <- side == 1
    vals <- if (isLon) grid("xmin", "xmax") else grid("ymin", "ymax")
    lines <- lapply(vals, function(d) if (isLon) line(rep(d, n), lonRange) else line(latRange, rep(d, n)))
    for (l in lines) terra::plot(l, add = TRUE, col = "grey60", lty = 3)
    at <- vapply(lines, function(l) {
      if (isLon) edgeAt(l, 2, e$ymin, e$ymax - e$ymin) else edgeAt(l, 1, e$xmin, e$xmax - e$xmin)
    }, numeric(1))
    ok <- !is.na(at)
    if (any(ok))
      graphics::axis(side, at = at[ok], pos = if (isLon) e$ymin else e$xmin, las = 1, cex.axis = 0.7, lwd = 0, lwd.ticks = 1,
                     labels = mapply(degLab, vals[ok], if (isLon) "E" else "N", if (isLon) "W" else "S"))
  }
  invisible(NULL)
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
