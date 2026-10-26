## Hindcast validation: what happened in the ELF core in one year, and figures that set the simulated
## rep-years beside it. The simulated side comes from fireSense_summary (per-rep rows written by its
## `hindcast` event); the observed side is computed here with the fit's own definitions, by calling the
## helpers the ELF gate and the escape fit use (`ELFfireCounts()`, `ELFfireOverlay()`, `escapedIgnitions()`).

#' Observed fire in an ELF core in one year
#'
#' What the hindcast simulation is compared with: natural-cause ignitions, escapes, area burned and the fire
#' sizes of `year`, all inside the ELF core (the buffer is excluded). Ignitions and escapes are
#' counted on the core with the fit's definitions: NFDB points with `CAUSE` in `naturalCauses` of any size
#' are ignitions ([ELFfireCounts()]), and those that [escapedIgnitions()] marks (one per fire polygon of the
#' year of at least `escapeSizeHa` and larger than a pixel) are escapes. Area burned is the core-clipped area of the fire polygons
#' ([ELFfireOverlay()]), given for the natural-cause polygons and for polygons of every cause, because the
#' simulation ignites natural fires only but the spread model was fitted to polygons of every cause.
#'
#' @param coreRaster `SpatRaster` of the ELF core on the simulation grid: cells outside the core are `NA`
#'   or 0 (`rasterToMatchELF` in the fireSense simulation).
#' @param firePoints `SpatVector` of fire points with columns `YEAR`, `CAUSE` and `SIZE_HA`
#'   (`ignitionFirePoints`).
#' @param firePolys `SpatVector` of fire perimeters with column `YEAR`, or a list of them as
#'   `firePolys` is in the simulation (`NULL` elements are dropped).
#' @param year The year to describe.
#' @param escapeSizeHa Size (ha) a fire must reach to count as escaped.
#' @param pixelAreaHa Area of one pixel (ha). Default: the area of a `coreRaster` pixel.
#' @param naturalCauses Values of `CAUSE` in `firePoints` that are natural ignitions.
#' @param polyCauseCol Column of `firePolys` holding the cause. If it is missing, the natural-cause values
#'   of the result are `NA`, with a warning.
#' @param naturalPolyCauses Values of `polyCauseCol` that are natural.
#'
#' @return A list. `counts` is a one-row `data.table`: `year`, `ignitions`, `escapes`, `areaBurnedHa`
#'   (natural-cause polygons in the core), `areaBurnedHaAllCause`, `nFires50` and `nFires50AllCause`
#'   (clipped pieces of at least `escapeSizeHa`). `sizes` is a `data.table` of the clipped pieces, with
#'   `cause` (`"natural"` or `"allCause"`) and `ha`. `polys` is the `SpatVector` of the year's natural-cause
#'   polygons (not clipped; `NULL` without polygons or a cause column), to draw on a map. `points` is the
#'   `SpatVector` of the natural ignition points that `counts$ignitions` counts (in the core's CRS, with
#'   `escaped`), the ignition locations the simulation is compared with.
#'
#' @export
#' @importFrom data.table data.table rbindlist
#' @importFrom terra ifel
hindcastObserved <- function(coreRaster, firePoints, firePolys, year, escapeSizeHa = 50,
                             pixelAreaHa = prod(terra::res(coreRaster)) / 1e4,
                             naturalCauses = c("L", "N"), polyCauseCol = "FIRECAUS",
                             naturalPolyCauses = "Natural") {
  stopifnot(inherits(coreRaster, "SpatRaster"), inherits(firePoints, "SpatVector"),
            length(year) == 1L)
  if (is.list(firePolys)) {
    firePolys <- Filter(Negate(is.null), firePolys)
    firePolys <- if (length(firePolys)) do.call(rbind, unname(firePolys))
  }
  ## the core is zone "core" (value 2) of a one-layer ELF raster, as ELFfireOverlay() reads it
  core <- terra::ifel(!is.na(coreRaster) & coreRaster > 0, 2, 0)
  names(core) <- "core"

  coreOverlay <- function(polys) if (is.null(polys)) .emptyOverlay() else ELFfireOverlay(core, polys, year)
  allCause <- coreOverlay(firePolys)
  noCause <- !is.null(firePolys) && !polyCauseCol %in% names(firePolys)
  if (noCause) warning("firePolys has no column '", polyCauseCol, "': the natural-cause area is NA")
  naturalPolys <- if (is.null(firePolys) || noCause) NULL else
    firePolys[firePolys[[polyCauseCol]][, 1] %in% naturalPolyCauses & as.integer(firePolys$YEAR) == year, ]
  natural <- coreOverlay(naturalPolys)
  inCore <- function(o) o$ha[o$zone == "core"]
  big <- function(o) as.integer(sum(inCore(o) >= escapeSizeHa))

  ## escapes as the escape fit defines them: natural ignitions that became a polygon of the year
  naturalPts <- firePoints[firePoints$CAUSE %in% naturalCauses, ]
  yearPolys <- if (!is.null(firePolys)) firePolys[as.integer(firePolys$YEAR) == year, ]
  naturalPts$.escaped <- escapedIgnitions(naturalPts, yearPolys, escapeSizeHa, pixelAreaHa)
  ## the polygons are not counted here, so the all-cause overlay only satisfies the signature
  pts <- ELFfireCounts(core, naturalPts, firePolys = NULL, fireYears = year, pixelAreaHa = pixelAreaHa,
                       naturalCauses = naturalCauses, escapeSizeHa = escapeSizeHa, overlay = allCause,
                       returnPoints = TRUE)
  pts$points$escaped <- pts$points$.escaped
  pts$points$.escaped <- NULL
  sizesOf <- function(cause, o) data.table::data.table(cause = rep(cause, length(inCore(o))), ha = inCore(o))
  list(
    counts = data.table::data.table(
      year = as.integer(year), ignitions = pts$counts$naturalIgnitions, escapes = sum(pts$points$escaped),
      areaBurnedHa = if (noCause) NA_real_ else sum(inCore(natural)),
      areaBurnedHaAllCause = sum(inCore(allCause)),
      nFires50 = if (noCause) NA_integer_ else big(natural),
      nFires50AllCause = big(allCause)),
    sizes = data.table::rbindlist(list(sizesOf("natural", natural), sizesOf("allCause", allCause))),
    polys = naturalPolys,
    points = pts$points
  )
}

#' One simulated year in an ELF core
#'
#' The simulated counterpart of [hindcastObserved()], for one simulated year: the year's ignitions and
#' escapes (`ignitionsAndEscapes` from `fireSense_ignitionPredict`) and each fire of the year's burn-ID raster
#' (`rstAnnualBurnID` from `fireSense_burn`, one integer per fire), counted or clipped to the ELF core. A fire
#' that burned only outside the core is not listed. The simulation burns the buffer too, so these sizes are
#' smaller than `burnSummary$areaBurnedHa` for fires crossing the core edge.
#'
#' @param burnID `SpatRaster` of fire IDs (`NA` where nothing burned).
#' @param coreRaster `SpatRaster` of the ELF core; cells outside it are `NA` or 0. It need not share the grid
#'   of `burnID`: pixels are placed by their centres.
#' @param ignitions `data.table` with the ignition cell numbers of `burnID` in `pixelID` and, in `escaped`,
#'   whether each ignition escaped. `NULL` gives `NA` ignitions and escapes.
#' @param escapeSizeHa Size (ha) at which a fire counts in `nFires50`.
#'
#' @return A list. `counts` is a one-row `data.table`: `ignitions`, `escapes`, `areaBurnedHa` and `nFires50`,
#'   the columns of [hindcastObserved()]'s `counts` of the same names. `sizes` has one row per fire in the
#'   core: `fire_id` and `ha`.
#'
#' @export
#' @importFrom data.table data.table
#' @importFrom terra extract values xyFromCell
hindcastSimulated <- function(burnID, coreRaster, ignitions = NULL, escapeSizeHa = 50) {
  stopifnot(inherits(burnID, "SpatRaster"), inherits(coreRaster, "SpatRaster"))
  inCore <- function(cells) {
    v <- terra::extract(coreRaster, terra::xyFromCell(burnID, cells))[[1]]
    !is.na(v) & v > 0
  }
  ids <- terra::values(burnID, mat = FALSE)
  burned <- which(!is.na(ids))
  n <- table(ids[burned][inCore(burned)])
  sizes <- data.table::data.table(fire_id = as.integer(names(n)),
                                  ha = as.numeric(n) * prod(terra::res(burnID)) / 1e4)
  here <- if (!is.null(ignitions)) inCore(ignitions$pixelID)
  list(
    counts = data.table::data.table(
      ignitions = if (is.null(here)) NA_integer_ else sum(here),
      escapes = if (is.null(here)) NA_integer_ else sum(here & ignitions$escaped %in% TRUE),
      areaBurnedHa = sum(sizes$ha), nFires50 = sum(sizes$ha >= escapeSizeHa)),
    sizes = sizes
  )
}

#' Simulated rep-years beside the observed year
#'
#' One panel each for the number of ignitions, escapes, area burned in the core and fires of at least the
#' escape size: the histogram of the simulated rep-years with the observed value as a vertical line. Area
#' burned and fire count have two lines, the natural-cause and the all-cause observed value.
#'
#' @param reps `data.table` with one row per simulated rep-year and columns `ignitions`, `escapes`,
#'   `areaBurnedHa` and `nFires50` (the rows `fireSense_summary` writes to `hindcast_reps.csv`, pooled over
#'   parameter sets).
#' @param observed The result of [hindcastObserved()].
#' @param title Plot title, e.g. the ELF and year.
#'
#' @return A `ggplot`.
#'
#' @export
#' @importFrom data.table as.data.table melt rbindlist
plotHindcastDistributions <- function(reps, observed, title = NULL) {
  value <- panel <- line <- NULL # data.table NSE
  panels <- c(ignitions = "ignitions", escapes = "escapes", areaBurnedHa = "area burned in core (ha)",
              nFires50 = "fires >= escape size")
  long <- function(d) data.table::melt(data.table::as.data.table(d)[, lapply(.SD, as.numeric), .SDcols = names(panels)],
                                       measure.vars = names(panels), variable.name = "panel")
  sim <- long(reps)
  obs <- data.table::rbindlist(list(
    long(observed$counts)[, line := "natural cause"],
    data.table::data.table(panel = c("areaBurnedHa", "nFires50"), line = "all causes",
                           value = unlist(observed$counts[, c("areaBurnedHaAllCause", "nFires50AllCause"), with = FALSE]))
  ), use.names = TRUE)
  sim[, panel := factor(panel, names(panels), panels)]
  obs[, panel := factor(panel, names(panels), panels)]
  ggplot2::ggplot(sim, ggplot2::aes(x = value)) +
    ggplot2::geom_histogram(bins = 20, fill = "grey70", colour = "grey30") +
    ggplot2::geom_vline(data = obs[!is.na(value)], ggplot2::aes(xintercept = value, colour = line),
                        linewidth = 1) +
    ggplot2::facet_wrap(~panel, scales = "free") +
    ggplot2::labs(x = NULL, y = "simulated rep-years", colour = "observed", title = title)
}

#' Fire-size distribution, simulated and observed
#'
#' The empirical cumulative distribution of core-clipped fire sizes (at least `minHa`), pooled over the
#' simulated rep-years, with the observed natural-cause and all-cause fire sizes as step lines.
#'
#' @param simSizes Numeric vector of simulated fire sizes (ha) pooled over rep-years (the `ha` column of
#'   `hindcast_sizes.csv`).
#' @param observed The result of [hindcastObserved()].
#' @param minHa Sizes below this are left out of both sides (default: 50 ha, the escape size).
#' @param title Plot title.
#'
#' @return A `ggplot`.
#'
#' @export
plotHindcastSizes <- function(simSizes, observed, minHa = 50, title = NULL) {
  ha <- NULL # ggplot2 NSE
  d <- rbind(data.frame(source = rep("simulated", length(simSizes)), ha = simSizes),
             data.frame(source = ifelse(observed$sizes$cause == "natural", "observed, natural cause",
                                        "observed, all causes"), ha = observed$sizes$ha))
  d <- d[d$ha >= minHa, ]
  ggplot2::ggplot(d, ggplot2::aes(x = ha, colour = source)) +
    ggplot2::stat_ecdf(geom = "step") +
    ggplot2::scale_x_log10() +
    ggplot2::labs(x = "fire size in core (ha)", y = "cumulative proportion", colour = NULL, title = title)
}

#' Burn-frequency map with the observed fires on top
#'
#' The cumulative burn map of the simulated rep-years divided by their number (the fraction of rep-years
#' each pixel burned), with the outlines of the observed fire polygons drawn over it.
#'
#' @param burnMap `SpatRaster` count of rep-years in which each pixel burned.
#' @param nReps Number of rep-years the count is over.
#' @param firePolys `SpatVector` of the observed polygons to draw, or `NULL`.
#' @param coreRaster Optional `SpatRaster` of the ELF core: the map is masked to it.
#' @param title Plot title.
#'
#' @return A `ggplot`.
#'
#' @export
#' @importFrom terra as.data.frame crop crs geom mask project
plotHindcastMap <- function(burnMap, nReps, firePolys = NULL, coreRaster = NULL, title = NULL) {
  x <- y <- frequency <- geom <- part <- NULL # ggplot2 NSE
  freq <- burnMap / nReps
  if (!is.null(coreRaster)) freq <- terra::mask(terra::crop(freq, coreRaster), coreRaster)
  d <- terra::as.data.frame(freq, xy = TRUE, na.rm = TRUE)
  names(d)[3] <- "frequency"
  g <- ggplot2::ggplot() +
    ggplot2::geom_raster(data = d, ggplot2::aes(x = x, y = y, fill = frequency)) +
    ggplot2::scale_fill_viridis_c(name = "fraction of\nrep-years burned", limits = c(0, 1)) +
    ggplot2::coord_equal() + ggplot2::labs(x = NULL, y = NULL, title = title)
  if (!is.null(firePolys) && nrow(firePolys)) {
    outline <- as.data.frame(terra::geom(terra::project(firePolys, terra::crs(burnMap))))
    g <- g + ggplot2::geom_path(data = outline, ggplot2::aes(x = x, y = y, group = interaction(geom, part)),
                                colour = "red", linewidth = 0.4)
  }
  g
}
