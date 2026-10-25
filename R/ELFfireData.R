## ELF fire counts and the zero/few-fire gate.
##
## An ELF with no fire cannot be fitted: fireSense_dataPrepFit stops with "no ignitions present"
## (fireSense_dataPrepFit.R:1595-1597) or with no fire polygons (:1537-1541), and it does so only
## after a full cold run. Counting fires while the ELF maps are built removes those ELFs before they
## enter the queue, and names the thin ones while there is still a chance to merge them.
##
## Counts use the same definitions as the fit: natural-cause NFDB ignitions (L/N,
## fireSense_dataPrepFit.R:1595) and NBAC polygons of every cause (getFirePolygons applies no cause
## filter) whose area clipped to the ELF exceeds one pixel. Fire records are read by the caller with
## fireregimetools' loaders, the same ones fireSense_dataPrepFit uses (its R/fireRecords.R:21,34),
## so a release change moves both together.

#' Count fires in each ELF
#'
#' Counts natural-cause ignitions, escaped natural ignitions and fire polygons falling in
#' each ELF's fitting area (core plus buffer), per year. Here an escaped ignition is one whose
#' reported size [escapedFires()] keeps, a proxy that needs no polygons: the escape fit and the hindcast
#' use [escapedIgnitions()], the ignitions that became a fire polygon of at least `escapeSizeHa`.
#'
#' Points are counted on the grid of `rasWhole` with a single [terra::cellFromXY()]
#' lookup. Polygons are intersected with each ELF's area and their clipped area is
#' compared with `pixelAreaHa`, mirroring the size filter the spread fit applies
#' (`fireSense_dataPrepFit.R:966-982`).
#'
#' @param rasWhole `SpatRaster` of ELF layers on one common grid, one layer per ELF, as
#'   built by [makeELFs()] (`sim$ELFs$rasWhole`). Cells are 0 (outside), 1 (buffer) or
#'   2 (core).
#' @param firePoints `SpatVector` of fire points, as returned by
#'   `fireregimetools::load_nfdb_points()`. Needs columns `YEAR`, `CAUSE` and `SIZE_HA`.
#' @param firePolys `SpatVector` of fire perimeters, as returned by
#'   `fireregimetools::load_nbac_polys()`. Needs column `YEAR`. May be `NULL` when `overlay` is given.
#' @param fireYears integer vector of the years being fitted.
#' @param pixelAreaHa area of one fitting pixel, in hectares. A polygon counts when its
#'   area inside the ELF exceeds this.
#' @param naturalCauses causes treated as natural ignitions.
#' @param overlay the table from [ELFfireOverlay()] for `rasWhole`, `firePolys` and `fireYears`, or `NULL`
#'   (default) to compute it. Passing it lets the clipping be done once and cached.
#' @param escapeSizeHa size (ha) a fire must reach to count as escaped, as in [escapedFires()]
#'   (`pixelAreaHa` is its pixel area).
#' @param returnPoints if `TRUE`, also return the counted points (see Value).
#'
#' @return A `data.table` with one row per ELF and year: `ELF`, `year`,
#'   `naturalIgnitions`, `escapes`, `firePolygons`. With `returnPoints = TRUE`, a list instead: `counts` (that
#'   table) and `points`, the `SpatVector` of the natural points that were counted (in `fireYears` and inside at least one
#'   ELF), in the CRS of `rasWhole`, with a logical column `escaped`.
#'
#' @export
#' @importFrom data.table data.table rbindlist setkeyv
#' @importFrom terra cellFromXY crds crs project
ELFfireCounts <- function(rasWhole, firePoints, firePolys, fireYears,
                          pixelAreaHa, naturalCauses = c("L", "N"), escapeSizeHa = 50,
                          overlay = NULL, returnPoints = FALSE) {
  stopifnot(
    inherits(rasWhole, "SpatRaster"),
    inherits(firePoints, "SpatVector"),
    inherits(firePolys, "SpatVector") || !is.null(overlay),
    length(fireYears) > 0,
    is.numeric(pixelAreaHa), length(pixelAreaHa) == 1
  )

  ELFnames <- names(rasWhole)
  fireYears <- sort(unique(as.integer(fireYears)))

  ## Points: one cellFromXY on the common grid, then index each layer's values.
  natural <- firePoints[firePoints$CAUSE %in% naturalCauses, ]
  natural$.row <- seq_len(nrow(natural))
  isEscape <- natural$.row %in% escapedFires(natural, escapeSizeHa, pixelAreaHa)$.row
  natural <- terra::project(natural, terra::crs(rasWhole))
  ptCell <- terra::cellFromXY(rasWhole, terra::crds(natural))
  ptYear <- as.integer(natural$YEAR)
  keep <- !is.na(ptCell) & ptYear %in% fireYears
  ptCell <- ptCell[keep]
  ptYear <- ptYear[keep]
  ptEscape <- isEscape[keep]
  counted <- natural[keep, ]
  counted$escaped <- ptEscape
  inAny <- rep(FALSE, length(ptCell))

  ## Polygons: clipped area per ELF, so a fire straddling an edge is judged on the part inside.
  if (is.null(overlay)) overlay <- ELFfireOverlay(rasWhole, firePolys, fireYears)
  overlay <- overlay[overlay$zone == "fit" & overlay$ha > pixelAreaHa, ]

  counts <- lapply(ELFnames, function(elf) {
    inArea <- rasWhole[[elf]][ptCell][[1]] > 0
    inArea[is.na(inArea)] <- FALSE
    inAny <<- inAny | inArea
    ign <- table(factor(ptYear[inArea], levels = fireYears))
    esc <- table(factor(ptYear[inArea & ptEscape], levels = fireYears))

    data.table::data.table(
      ELF = elf,
      year = fireYears,
      naturalIgnitions = as.integer(ign),
      escapes = as.integer(esc),
      firePolygons = as.integer(table(factor(overlay$YEAR[overlay$ELF == elf], levels = fireYears)))
    )
  })

  counts <- data.table::rbindlist(counts)
  data.table::setkeyv(counts, c("ELF", "year"))
  if (returnPoints) list(counts = counts[], points = counted[inAny, ]) else counts[]
}

#' Clip fire polygons to each ELF, once
#'
#' Intersects the fire polygons with each ELF's fitting area (core plus buffer) and with its core, and
#' returns one row per clipped piece. [ELFfireCounts()] counts the pieces larger than a pixel in the
#' fitting area, and [ELFregionStats()] sums the core pieces into a burn rate, so the slow intersection
#' is done once for both. A polygon wholly inside a zone is not clipped; only those crossing its edge
#' are, which is most of the saving.
#'
#' @param rasWhole `SpatRaster` of ELF layers on one common grid, as in [ELFfireCounts()].
#' @param firePolys `SpatVector` of fire perimeters with column `YEAR`.
#' @param fireYears integer vector of the years to keep.
#'
#' @return A `data.table`: `ELF`, `zone` (`"fit"` for core plus buffer, `"core"`), `YEAR` and `ha`, the
#'   geodesic area of the piece. ELFs no polygon reaches have no rows.
#'
#' @export
#' @importFrom data.table data.table rbindlist
#' @importFrom terra as.polygons crs expanse intersect project relate
ELFfireOverlay <- function(rasWhole, firePolys, fireYears) {
  stopifnot(inherits(rasWhole, "SpatRaster"), inherits(firePolys, "SpatVector"))
  fireYears <- sort(unique(as.integer(fireYears)))
  polys <- terra::project(firePolys, terra::crs(rasWhole))
  polys <- polys[as.integer(polys$YEAR) %in% fireYears, ]
  empty <- .emptyOverlay()
  if (!nrow(polys)) return(empty)

  ## A polygon wholly inside a zone is not clipped: its own area is its piece. Only polygons that cross
  ## the zone's edge are intersected with it, which is the slow step (fire perimeters have thousands of
  ## vertices).
  polys$ha <- terra::expanse(polys, unit = "ha")
  ## Bounding boxes, once: a zone is tested only against the polygons whose box meets its own, because
  ## testing all of them against each of the ELFs' zones is slow.
  x <- y <- NULL # data.table NSE
  vertices <- data.table::as.data.table(terra::geom(polys))
  boxes <- vertices[, list(xmin = min(x), xmax = max(x), ymin = min(y), ymax = max(y)), by = "geom"]
  data.table::setorderv(boxes, "geom")
  clipTo <- function(area) {
    area <- area[as.logical(area[[1]][, 1]), ]
    if (!nrow(area)) return(NULL)
    e <- as.vector(terra::ext(area))
    near <- which(boxes$xmax >= e[["xmin"]] & boxes$xmin <= e[["xmax"]] &
                    boxes$ymax >= e[["ymin"]] & boxes$ymin <= e[["ymax"]])
    candidates <- polys[near, ]
    hit <- candidates[as.vector(terra::relate(candidates, area, "intersects")), ]
    if (!nrow(hit)) return(NULL)
    inside <- as.vector(terra::relate(hit, area, "within"))
    ## the intersect warns when nothing overlaps
    crossing <- if (any(!inside)) suppressWarnings(terra::intersect(hit[!inside, ], area))
    list(YEAR = c(as.integer(hit$YEAR[inside]), if (!is.null(crossing)) as.integer(crossing$YEAR)),
         ha = c(hit$ha[inside], if (!is.null(crossing)) terra::expanse(crossing, unit = "ha")))
  }
  rows <- lapply(names(rasWhole), function(elf) {
    elfRas <- rasWhole[[elf]]
    out <- list(fit = clipTo(terra::as.polygons(elfRas > 0, dissolve = TRUE)),
                core = clipTo(terra::as.polygons(elfRas == 2, dissolve = TRUE)))
    out <- out[!vapply(out, is.null, logical(1))]
    if (!length(out)) return(NULL)
    data.table::rbindlist(lapply(names(out), function(z)
      data.table::data.table(ELF = elf, zone = z, YEAR = out[[z]]$YEAR, ha = out[[z]]$ha)))
  })
  rows <- rows[!vapply(rows, is.null, logical(1))]
  if (!length(rows)) return(empty)
  data.table::rbindlist(rows)
}

#' Classify each ELF as zero, few or ok on its fire counts
#'
#' An ELF with no natural ignitions, or no fire polygons, over the whole fitting window
#' cannot be fitted and is `"zero"`. One below any threshold is `"few"`: it can be
#' fitted, but the thin end of the data drives the escape model's 5-fold CV
#' (`fireSense_ignitionFit.R:463,889`), which needs at least 10 ignited cell-years
#' before some fold is left holding a single row. The escape and spread fits use only escaped
#' fires (see [escapedFires()]), so an ELF with many small fires and few large ones is `"few"`
#' through `minEscapes`.
#'
#' @param counts `data.table` from [ELFfireCounts()].
#' @param minNaturalIgnitions,minFirePolygons,minEscapes thresholds below which an ELF is `"few"`.
#'
#' @return A `data.table` with one row per ELF: `ELF`, `naturalIgnitions`, `escapes`,
#'   `firePolygons`, `yearsWithFire`, `status` (`"zero"`, `"few"` or `"ok"`).
#'
#' @export
#' @importFrom data.table as.data.table fifelse setkeyv
ELFfitStatus <- function(counts, minNaturalIgnitions = 50, minFirePolygons = 50, minEscapes = 5) {
  stopifnot(
    is.data.frame(counts),
    all(c("ELF", "year", "naturalIgnitions", "escapes", "firePolygons") %in% names(counts))
  )

  naturalIgnitions <- escapes <- firePolygons <- status <- NULL # data.table NSE

  out <- data.table::as.data.table(counts)[
    , list(
      naturalIgnitions = sum(naturalIgnitions),
      escapes = sum(escapes),
      firePolygons = sum(firePolygons),
      yearsWithFire = sum(naturalIgnitions > 0 | firePolygons > 0)
    ),
    by = "ELF"
  ]

  out[, status := data.table::fifelse(
    naturalIgnitions == 0 | firePolygons == 0, "zero",
    data.table::fifelse(
      naturalIgnitions < minNaturalIgnitions | firePolygons < minFirePolygons |
        escapes < minEscapes,
      "few", "ok"
    )
  )]

  data.table::setkeyv(out, "ELF")
  out[]
}

#' ELFs that cannot be fitted
#'
#' @param status `data.table` from [ELFfitStatus()].
#'
#' @return `character` vector of ELF names with status `"zero"`.
#'
#' @export
ELFsExcluded <- function(status) {
  stopifnot(is.data.frame(status), all(c("ELF", "status") %in% names(status)))
  as.character(status$ELF[status$status == "zero"])
}

#' Measure flammable area in each ELF
#'
#' Counts the land-cover pixels an ELF has that the fit would treat as flammable, i.e.
#' those not in `nonflammableLCC`. Unclassified cover (VLCE2 code `0`) is nonflammable
#' by default, which is why ELFs outside VLCE2's mapped extent measure at or near zero
#' (8.2, 10.3.2 and 3.2.4 are 0; `no-tree-elfs.md` §2, `elf-8.2-ecoregion.md`).
#'
#' This does **not** gate anything. Per Eliot's 2026-09-11 decision there is no
#' flammable floor -- `fireSense_spreadFit` is allowed to try regardless -- so the result
#' is a report, and an input to deciding whether merging an ELF with a neighbour would
#' give it enough to fit.
#'
#' @param rasWhole `SpatRaster` of ELF layers, as in [ELFfireCounts()].
#' @param rstLCC `SpatRaster` of land cover, the same source the fit uses.
#' @param nonflammableLCC land-cover codes that cannot burn. The default,
#'   [fireSenseNonflammableLCC], matches `fireSense_dataPrepFit`'s `nonflammableLCC`
#'   parameter.
#'
#' @return A `data.table` with one row per ELF: `ELF`, `flammablePixels`,
#'   `totalPixels`, `flammableFraction`, `flammableAreaHa`.
#'
#' @export
#' @importFrom data.table data.table rbindlist setkeyv
#' @importFrom terra crs project res values
ELFflammableArea <- function(rasWhole, rstLCC, nonflammableLCC = fireSenseNonflammableLCC) {
  stopifnot(inherits(rasWhole, "SpatRaster"), inherits(rstLCC, "SpatRaster"))

  ## Land cover is normally finer than the ELF grid, so move the ELF footprint onto the
  ## land-cover grid rather than the reverse: resampling cover codes would invent classes.
  lccVals <- terra::values(rstLCC, mat = FALSE)
  pixelAreaHa <- prod(terra::res(rstLCC)) / 1e4

  out <- lapply(names(rasWhole), function(elf) {
    elfOnLCC <- terra::project(rasWhole[[elf]], rstLCC, method = "near")
    inELF <- terra::values(elfOnLCC, mat = FALSE) > 0
    inELF[is.na(inELF)] <- FALSE

    lccHere <- lccVals[inELF]
    lccHere <- lccHere[!is.na(lccHere)]
    flammable <- sum(!lccHere %in% nonflammableLCC)

    data.table::data.table(
      ELF = elf,
      flammablePixels = flammable,
      totalPixels = length(lccHere),
      flammableFraction = if (length(lccHere)) flammable / length(lccHere) else NA_real_,
      flammableAreaHa = flammable * pixelAreaHa
    )
  })

  out <- data.table::rbindlist(out)
  data.table::setkeyv(out, "ELF")
  out[]
}

## The table ELFfireOverlay() returns when nothing is clipped.
.emptyOverlay <- function() {
  data.table::data.table(ELF = character(0), zone = character(0), YEAR = integer(0), ha = numeric(0))
}
