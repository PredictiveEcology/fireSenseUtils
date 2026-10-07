utils::globalVariables(c(
  ".N", "goalSize", "simSize", "pixels"
))

#' Create buffers around polygons based on area target for buffer
#'
#' @param poly `sf` polygons or a list of `sf` containing polygons to buffer.
#' @param rasterToMatch A `SpatRaster` with `res`, `origin`, `extent`,
#'   `crs` of desired outputted `pixelID` values.
#' @param areaMultiplier Either a scalar that will buffer `areaMultiplier * fireSize` or
#'   a function of `fireSize.` Default is 1. See [multiplier()] for an example.
#' @param verb Logical or numeric related to how much verbosity is printed. `FALSE` or
#'   `0` is none. `TRUE` or `1` is some. `2` is much more.
#' @param polyName Optional character string of the polygon layer name (not the individual polygons
#'   on a `sf` polygon object)
#' @param field Passed to `fasterize::fasterize`. If this is unique (such as polygon id),
#'   then each polygon will have its buffer calculated independently for each unique value
#'   in `field`
#' @param minSize The absolute minimum size of the buffer & non-buffer together. This will
#'   be imposed after `areaMultiplier`.
#' @param cores number of processor cores to use
#' @param ... passed to `fasterize::fasterize`
#'
#' @return
#' A `data.table` (or list of `data.table`s if `poly` was a list) with 2 columns:
#' `buffer` and `pixelID`. `buffer` is either `1` (the original polygon) or
#' `0` (in the buffer).
#'
#' @export
#' @rdname bufferToArea
bufferToArea <- function(poly, rasterToMatch, areaMultiplier,
                         verb = FALSE, polyName = NULL, field = NULL,
                         minSize = 500, cores = 1, ...) {
  UseMethod("bufferToArea")
}

#' @export
#' @importFrom parallelly availableCores
#' @importFrom purrr pmap
#'
#' @rdname bufferToArea
bufferToArea.list <- function(poly, rasterToMatch, areaMultiplier = 10,
                              verb = FALSE, polyName = NULL, field = NULL,
                              minSize = 500, cores = 1, ...) {
  if (is.null(polyName)) {
    polyName <- names(poly)
  }
  maxCores <- parallelly::availableCores(constraints = "connections", omit = 1)
  cores <- min(min(length(poly), cores), maxCores)
  ## Buffer pixels are sampled at random. A forked child seeds itself independently, so draw one seed per
  ## polygon set here and run each set under it: the same session seed then gives the same buffers,
  ## forked or not.
  seeds <- sample.int(.Machine$integer.max, length(poly))
  if (cores > 1 && !"tools:rstudio" %in% search()) { # no forked workers inside RStudio
    out <- parallel::mcMap(
      mc.cores = cores,
      poly = poly,
      polyName = polyName,
      .seed = seeds,
      MoreArgs = list(
        rasterToMatch = rasterToMatch, verb = verb,
        areaMultiplier = areaMultiplier, field = field, minSize = minSize,
        cores = 1,
        ...
      ),
      .singleThreadedGDAL(.withSeed(bufferToArea))
    )
  } else {
    out <- purrr::pmap(
      .l = list(poly = poly, polyName = polyName, .seed = seeds),
      rasterToMatch = rasterToMatch, verb = verb,
      areaMultiplier = areaMultiplier, field = field, minSize = minSize,
      cores = 1,
      ...,
      .f = .withSeed(bufferToArea)
    )
  }
  names(out) <- names(poly)
  out
}

## Wrap `FUN` so that it runs under the seed passed as `.seed`, restoring the caller's stream afterwards.
.withSeed <- function(FUN) {
  function(..., .seed) {
    withr::with_seed(.seed, FUN(...))
  }
}

## Wrap `FUN` for a forked child. GDAL keeps one process-wide worker-thread pool, created at the
## first multi-threaded raster write (terra passes `NUM_THREADS = terraOptions()$threads`). A child
## forked after that inherits the pool but none of its threads, so its own multi-threaded write
## waits for them forever with 0 CPU. Writing single-threaded in the child never uses the pool.
.singleThreadedGDAL <- function(FUN) {
  function(...) {
    terra::terraOptions(threads = 1)
    FUN(...)
  }
}

#' @export
#' @importFrom sf st_as_sf
#' @rdname bufferToArea
bufferToArea.SpatialPolygons <- function(poly, rasterToMatch, areaMultiplier = 10,
                                         verb = FALSE, polyName = NULL, field = NULL,
                                         minSize = 500, cores = 1, ...) {
  bufferToArea.sf(sf::st_as_sf(poly), rasterToMatch,
    areaMultiplier = areaMultiplier,
    verb = verb, polyName = polyName, field = field,
    minSize = minSize, cores = cores, ...
  )
}

#' @export
#' @importFrom data.table data.table rbindlist setorderv
#' @importFrom LandR asInteger
#' @importFrom terra rasterize values
#' @importFrom sf st_crs st_transform
#' @rdname bufferToArea
bufferToArea.sf <- function(poly, rasterToMatch, areaMultiplier = 10,
                            verb = FALSE, polyName = NULL, field = NULL,
                            minSize = 500, cores = 1, ...) {
  if (is.null(polyName)) polyName <- "Layer 1"
  if (as.integer(verb) >= 1) print(paste("Buffering polygons on", polyName))
  r <- rasterize(
    x = sf::st_transform(poly, sf::st_crs(rasterToMatch)),
    y = rasterToMatch, field = field, ...
  )

  emptyDT <- data.table(pixelID = integer(0), buffer = integer(0), ids = integer(0))

  if (all(is.na(r[]))) {
    return(emptyDT)
  }
  rvals <- values(r, mat = FALSE)
  loci <- which(!is.na(rvals))
  ids <- as.integer(rvals[loci])

  am <- if (is(areaMultiplier, "function")) {
    areaMultiplier
  } else {
    function(x) areaMultiplier * x
  }
  fireIds <- unique(ids)
  nFires <- length(fireIds)
  fire <- match(ids, fireIds) # fire (index into `fireIds`) of each burned cell, `loci`
  actualSize <- tabulate(fire, nFires)
  goalSize <- asInteger(pmax(minSize, vapply(actualSize, am, numeric(1))))

  ## A buffer grows one ring (the 8 neighbours of its newest cells) per iteration. Only the newest ring
  ## can reach a new cell, so each iteration looks at that ring alone, not at the whole buffer. A cell
  ## belongs to the one fire that claims it first, as in SpaDES.tools::spread2(), which picks at random
  ## among the fires next to a contested cell.
  owner <- integer(terra::ncell(r)) # 0L = unclaimed
  owner[loci] <- fire
  burned <- split(loci, factor(fire, levels = seq_len(nFires)))
  ## cells of each fire, one element per ring, burned cells first
  rings <- lapply(burned, list)
  size <- actualSize
  active <- rep(TRUE, nFires)
  frontier <- loci
  ## Buffer rows of fire `k` from the cells `cells`, flagging the burned cells
  bufferRows <- function(k, cells) {
    data.table(buffer = as.integer(cells %in% burned[[k]]), pixelID = cells, ids = fireIds[k])
  }
  out <- list()

  it <- 1L
  ## A buffer grows one cell per iteration, so it can take at most this many to cover the raster.
  maxIts <- max(terra::nrow(r), terra::ncol(r)) + 1L
  while (any(active) && it <= maxIts) {
    pairs <- terra::adjacent(r, frontier, directions = 8, pairs = TRUE)
    pairs <- pairs[owner[pairs[, "to"]] == 0L, , drop = FALSE]
    pairs <- pairs[sample.int(nrow(pairs)), , drop = FALSE]
    pairs <- pairs[!duplicated(pairs[, "to"]), , drop = FALSE]
    newCells <- pairs[, "to"]
    newFire <- owner[pairs[, "from"]]
    owner[newCells] <- newFire
    grew <- tabulate(newFire, nFires)
    size <- size + grew
    newRing <- split(newCells, factor(newFire, levels = seq_len(nFires)))

    bigger <- active & size > goalSize
    ## A buffer that did not grow has filled the landscape (raster edge) before reaching its target:
    ## keep it, at the largest size the landscape allows, rather than lose the fire.
    stalled <- active & !bigger & grew == 0L
    for (k in which(stalled)) {
      out[[length(out) + 1L]] <- bufferRows(k, unlist(rings[[k]], use.names = FALSE))
    }
    for (k in which(bigger)) {
      earlier <- unlist(rings[[k]], use.names = FALSE)
      last <- newRing[[k]]
      needMore <- goalSize[k] - length(earlier)
      cells <- if (needMore > 0) {
        c(earlier, last[sample.int(length(last), needMore)])
      } else {
        earlier[sample.int(length(earlier), goalSize[k])]
      }
      out[[length(out) + 1L]] <- bufferRows(k, cells)
    }
    done <- bigger | stalled
    active <- active & !done
    for (k in which(active & grew > 0L)) rings[[k]] <- c(rings[[k]], newRing[k])
    ## the next frontier is the cells just claimed by fires still growing ...
    frontier <- newCells[active[newFire]]
    if (any(done)) {
      ## ... and a finished fire's cells are free again for the others, so the cells of fires still
      ## growing that touch them are on the frontier too
      freed <- unlist(c(lapply(rings[done], unlist, use.names = FALSE), newRing[done]), use.names = FALSE)
      owner[freed] <- 0L
      if (any(active) && length(freed) > 0L) {
        touching <- unique(terra::adjacent(r, freed, directions = 8, pairs = TRUE)[, "to"])
        frontier <- unique(c(frontier, touching[owner[touching] != 0L]))
      }
      rings[done] <- list(NULL)
    }
    it <- it + 1L
  }
  ## fires still growing when the iterations ran out
  for (k in which(active)) {
    out[[length(out) + 1L]] <- bufferRows(k, unlist(rings[[k]], use.names = FALSE))
  }
  out3 <- if (length(out) > 0) {
    rbindlist(out)
  } else {
    emptyDT
  }
  setorderv(out3, "buffer", order = -1L)
  out4 <- out3[, list(buffer = buffer[1], ids = ids[1]), by = "pixelID"]
}

#' Buffer-size multiplier for fire footprints
#'
#' Returns a per-fire buffer target (in pixels) by scaling `size` so that
#' small fires get proportionally larger buffers and large fires plateau
#' at `baseMultiplier * size`. The result is then clamped from below by
#' `minSize`. Used when generating non-burned "control" pixels around
#' historical fire perimeters.
#'
#' @param size Numeric vector. Fire sizes in pixels.
#' @param minSize Numeric scalar. Absolute floor on the returned size; the
#'   buffer (plus burned pixels) is never smaller than this. Default `1000`.
#' @param baseMultiplier Numeric scalar. Asymptotic multiplier for large
#'   fires (as `size` grows, the effective multiplier tends to
#'   `baseMultiplier`). Default `5`.
#'
#' @return Integer vector the same length as `size`, giving the target
#'   buffered area in pixels for each fire.
#'
#' @export
multiplier <- function(size, minSize = 1000, baseMultiplier = 5) {
  pmax(minSize, round(pmax(baseMultiplier, 14 - log(size)) * size, 0))
}

#' Remove buffered fires in `fireBufferedListDT` that are outside `flammableRTM`
#'
#' @param fireBufferedDT data.table containing indices for buffered annual fires
#'
#' @template flammableRTM
#'
#' @return `fireBufferedDT` excluding fires with indices (burned or unburned) outside `flammableRTM`
#'
#' @export
removeBufferedFiresOutsideRTM <- function(fireBufferedDT, flammableRTM) {
  fireBufferedDT[, flammable := flammableRTM[fireBufferedDT$pixelID]]
  toExclude <- unique(fireBufferedDT[is.na(flammable)]$ids)

  # THIS CHANGE FROM ELIOT SEPT 3, 2025
  #  The previous seems egregious -- it would remove every fire in its entirety
  #   if there was even one pixel that was outside the RTM.
  #  This first came up as a problem on Newfoundland, there there are inlets that make
  #   small areas of "outside RTM"... i.e., the ocean.
  #  Maybe it is OK for situations where there is no ocean nearby -- i.e., lakes would be
  #   non-flammable, but inside RTM
  fireBufferedDT <- fireBufferedDT[!is.na(flammable)]
  # fireBufferedDT <- fireBufferedDT[!ids %in% toExclude]


  fireBufferedDT[, flammable := NULL]
  return(fireBufferedDT)
}
