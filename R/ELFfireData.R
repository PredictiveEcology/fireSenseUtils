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
#' Counts natural-cause ignitions and fire polygons falling in each ELF's fitting
#' area (core plus buffer), per year.
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
#'   `fireregimetools::load_nfdb_points()`. Needs columns `YEAR` and `CAUSE`.
#' @param firePolys `SpatVector` of fire perimeters, as returned by
#'   `fireregimetools::load_nbac_polys()`. Needs column `YEAR`.
#' @param fireYears integer vector of the years being fitted.
#' @param pixelAreaHa area of one fitting pixel, in hectares. A polygon counts when its
#'   area inside the ELF exceeds this.
#' @param naturalCauses causes treated as natural ignitions.
#'
#' @return A `data.table` with one row per ELF and year: `ELF`, `year`,
#'   `naturalIgnitions`, `firePolygons`.
#'
#' @export
#' @importFrom data.table data.table rbindlist setkeyv
#' @importFrom terra cellFromXY crds crs project
ELFfireCounts <- function(rasWhole, firePoints, firePolys, fireYears,
                          pixelAreaHa, naturalCauses = c("L", "N")) {
  stopifnot(
    inherits(rasWhole, "SpatRaster"),
    inherits(firePoints, "SpatVector"),
    inherits(firePolys, "SpatVector"),
    length(fireYears) > 0,
    is.numeric(pixelAreaHa), length(pixelAreaHa) == 1
  )

  ELFnames <- names(rasWhole)
  fireYears <- sort(unique(as.integer(fireYears)))

  ## Points: one cellFromXY on the common grid, then index each layer's values.
  natural <- firePoints[firePoints$CAUSE %in% naturalCauses, ]
  natural <- terra::project(natural, terra::crs(rasWhole))
  ptCell <- terra::cellFromXY(rasWhole, terra::crds(natural))
  ptYear <- as.integer(natural$YEAR)
  keep <- !is.na(ptCell) & ptYear %in% fireYears
  ptCell <- ptCell[keep]
  ptYear <- ptYear[keep]

  ## Polygons: clipped area per ELF, so a fire straddling an edge is judged on the part inside.
  polys <- terra::project(firePolys, terra::crs(rasWhole))
  polys <- polys[as.integer(polys$YEAR) %in% fireYears, ]

  counts <- lapply(ELFnames, function(elf) {
    inArea <- rasWhole[[elf]][ptCell][[1]] > 0
    inArea[is.na(inArea)] <- FALSE
    ign <- table(factor(ptYear[inArea], levels = fireYears))

    data.table::data.table(
      ELF = elf,
      year = fireYears,
      naturalIgnitions = as.integer(ign),
      firePolygons = .ELFpolygonCounts(rasWhole[[elf]], polys, fireYears, pixelAreaHa)
    )
  })

  counts <- data.table::rbindlist(counts)
  data.table::setkeyv(counts, c("ELF", "year"))
  counts[]
}

## Fire polygons per year whose area inside one ELF exceeds one pixel. The ELF's area is the
## footprint of its non-zero cells, so this needs no polygon version of the ELF map.
#' @importFrom terra as.polygons expanse intersect
.ELFpolygonCounts <- function(elfRas, polys, fireYears, pixelAreaHa) {
  zeroes <- rep(0L, length(fireYears))
  area <- terra::as.polygons(elfRas > 0, dissolve = TRUE)
  area <- area[as.logical(area[[1]][, 1]), ]
  if (nrow(area) == 0 || nrow(polys) == 0) {
    return(zeroes)
  }
  clipped <- terra::intersect(polys, area)
  if (nrow(clipped) == 0) {
    return(zeroes)
  }
  big <- terra::expanse(clipped, unit = "ha") > pixelAreaHa
  as.integer(table(factor(as.integer(clipped$YEAR)[big], levels = fireYears)))
}

#' Classify each ELF as zero, few or ok on its fire counts
#'
#' An ELF with no natural ignitions, or no fire polygons, over the whole fitting window
#' cannot be fitted and is `"zero"`. One below either threshold is `"few"`: it can be
#' fitted, but the thin end of the data drives the escape model's 5-fold CV
#' (`fireSense_IgnitionFit.R:463,889`), which needs at least 10 ignited cell-years
#' before some fold is left holding a single row.
#'
#' @param counts `data.table` from [ELFfireCounts()].
#' @param minNaturalIgnitions,minFirePolygons thresholds below which an ELF is `"few"`.
#'
#' @return A `data.table` with one row per ELF: `ELF`, `naturalIgnitions`,
#'   `firePolygons`, `yearsWithFire`, `status` (`"zero"`, `"few"` or `"ok"`).
#'
#' @export
#' @importFrom data.table as.data.table fifelse setkeyv
ELFfitStatus <- function(counts, minNaturalIgnitions = 50, minFirePolygons = 50) {
  stopifnot(
    is.data.frame(counts),
    all(c("ELF", "year", "naturalIgnitions", "firePolygons") %in% names(counts))
  )

  naturalIgnitions <- firePolygons <- status <- NULL # data.table NSE

  out <- data.table::as.data.table(counts)[
    , list(
      naturalIgnitions = sum(naturalIgnitions),
      firePolygons = sum(firePolygons),
      yearsWithFire = sum(naturalIgnitions > 0 | firePolygons > 0)
    ),
    by = "ELF"
  ]

  out[, status := data.table::fifelse(
    naturalIgnitions == 0 | firePolygons == 0, "zero",
    data.table::fifelse(
      naturalIgnitions < minNaturalIgnitions | firePolygons < minFirePolygons,
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
