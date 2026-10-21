## Which ELFs have too few fires, which are too small, and which merge.
##
## Fire records are fetched and read as fireSense_dataPrepFit reads them for the fit (its R/fireRecords.R):
## archives through reproducible, shapefiles through fireregimetools, fires of every size. So the gate and
## the fit count the same fires.

#' Path to the shapefile in a fire-record archive
#'
#' Downloads the archive if needed. The shapefile's name carries the release, so callers key their
#' Cache on the name.
#'
#' @param url URL of the archive.
#' @param destinationPath Directory to download and extract into.
#' @param ... Passed to `reproducible::preProcess()`.
#'
#' @return Character path(s) of the `.shp` file(s) in the archive.
#' @export
fireRecordShapefile <- function(url, destinationPath, ...) {
  files <- reproducible::preProcess(url = url, destinationPath = destinationPath, fun = NA, ...)$targetFilePath
  grep("\\.shp$", files, value = TRUE)
}

#' Path to the land-cover map the size rule uses
#'
#' SCANFI v2 land cover for 2020 in Canada LCC class codes, the same map and source as
#' [LandR::prepInputs_SCANFI_LCC_FAO()] starts from. The file's name carries its release, so a copy in
#' `destinationPath` or `getOption("reproducible.inputPaths")` is used as it is; otherwise it is downloaded
#' from LandR's SCANFI mirror, which needs no Google login.
#'
#' @param destinationPath Directory to download into.
#'
#' @return Character path of the `.tif`.
#' @export
sizeRuleLandCover <- function(destinationPath) {
  source <- LandR:::.scanfiLCCFAOSource(2020, "V2")
  local <- file.path(c(destinationPath, getOption("reproducible.inputPaths")), source$targetFile)
  local <- local[file.exists(local)]
  ## normalized, as preProcess() returns it: callers Cache on this path, so a local copy and a download must match
  if (length(local)) return(normalizePath(local[1]))
  ## LandR's wrapper fetches the Drive id from its public SCANFI mirror, so no Google login is needed
  LandR:::.withSCANFIAccess(
    reproducible::preProcess(url = source$url, targetFile = source$targetFile,
                             destinationPath = destinationPath, fun = NA)$targetFilePath,
    what = "the SCANFI land cover map", dataYear = 2020, dataVersion = "V2", urlArg = "url")
}

## A fire-record archive read with a fireregimetools loader through prepInputs(), which caches the loaded
## object itself: the key holds the url, the release (the shapefile's name), and the loader call with the
## years, the size limit and the study area. `loader` is the quoted fireregimetools function.
.loadFireRecords <- function(url, destinationPath, loader, fireYears, studyArea) {
  shp <- fireRecordShapefile(url, destinationPath)
  e <- as.vector(terra::ext(studyArea))
  call <- bquote(.(loader)(
    targetFile, study_area = terra::as.polygons(terra::ext(.(e[["xmin"]]), .(e[["xmax"]]), .(e[["ymin"]]), .(e[["ymax"]])),
                                                crs = .(terra::crs(studyArea))),
    fire_years = .(as.integer(fireYears)), min_size_ha = 0))
  reproducible::prepInputs(url = url, targetFile = basename(shp), destinationPath = destinationPath,
                           fun = call, useCache = TRUE)
}

## The clipped NBAC pieces of every ELF, from the archive: loads the polygons only when the table is not
## cached. Cache keys on the ELF maps, the release and the years (and this code), not on the polygons.
.nbacOverlay <- function(rasWhole, nbacUrl, destinationPath, fireYears) {
  polys <- .loadFireRecords(nbacUrl, destinationPath, quote(fireregimetools::load_nbac_polys), fireYears,
                            terra::as.polygons(terra::ext(rasWhole), crs = terra::crs(rasWhole)))
  ELFfireOverlay(rasWhole, polys, fireYears)
}

#' Count fires per ELF and merge ELFs with too few
#'
#' Merges each ELF with too few fires with a neighbour that shares its base ([ELFmergePlan()]), then each
#' region still smaller than `minRegionAreaKm2` with a similar neighbour ([ELFsizePlan()]). Arctic and Prairie grassland ELFs
#' ([ELFsNotFitted()]) are neither counted nor merged.
#'
#' The fire records are read with `reproducible::prepInputs()`, which caches them, and the clipped fire
#' polygons of every ELF ([ELFfireOverlay()], the slow step) are cached on the ELF maps, the NBAC release
#' and `fireYears`, so a change to the plan rules does not repeat either.
#'
#' @param ELFs List from [makeELFs()].
#' @param fireYears Integer vector of years to count fires over.
#' @param pixelAreaHa Area of one pixel of the ELF rasters, in ha.
#' @param destinationPath Directory the fire records are downloaded to.
#' @param nfdbUrl,nbacUrl URLs of the NFDB point and NBAC polygon archives.
#' @param minNaturalIgnitions Fewer natural-cause ignitions than this is too few.
#' @param minFirePolygons Fewer fire polygons than this is too few.
#' @param escapeSizeHa Size (ha) a fire must reach to count as escaped.
#' @param minEscapes Fewer escaped natural fires than this is too few.
#' @param landCoverFile Path to a land-cover map in class codes (SCANFI's), or `NULL`. With
#'   `minRegionAreaKm2` it switches on the size rule.
#' @param minRegionAreaKm2 After the fire-count plan, regions with a core smaller than this (km2) are merged
#'   with a similar neighbour. `NULL` or `NA` turns the rule off.
#' @param maxLandCoverDist,maxBurnRatio How similar a neighbour must be for the size rule: largest land-cover
#'   (Bray-Curtis) distance and largest burn-rate ratio.
#'
#' @return List: `ELFs` (the merged maps), `status` ([ELFfitStatus()]), `plan` (the fire-count and size
#'   merges and the skips, in [ELFmergePlan()]'s format) and `excluded` (names of ELFs not fitted).
#'
#' @export
fewFireELFs <- function(ELFs, fireYears, pixelAreaHa, destinationPath,
                        nfdbUrl = nfdbPointUrl(), nbacUrl = latestNBACUrl(),
                        minNaturalIgnitions = 50, minFirePolygons = 50,
                        escapeSizeHa = 50, minEscapes = 5, landCoverFile = NULL,
                        minRegionAreaKm2 = 35000, maxLandCoverDist = 0.35, maxBurnRatio = 6) {
  ## Arctic and Prairie grassland ELFs are out permanently: not counted, not merged (ELFsNotFitted())
  ids <- setdiff(names(ELFs$rasWhole), ELFsNotFitted(names(ELFs$rasWhole)))
  rasWhole <- terra::rast(unname(ELFs$rasWhole[ids]))
  names(rasWhole) <- ids
  studyArea <- terra::as.polygons(terra::ext(rasWhole), crs = terra::crs(rasWhole))

  points <- .loadFireRecords(nfdbUrl, destinationPath, quote(fireregimetools::load_nfdb_points), fireYears, studyArea)
  nbacShp <- fireRecordShapefile(nbacUrl, destinationPath)
  overlay <- reproducible::Cache(
    .nbacOverlay(rasWhole, nbacUrl, destinationPath, fireYears),
    omitArgs = "destinationPath", .functionName = ".nbacOverlay", useCache = TRUE,
    .cacheExtra = list(basename(nbacShp), ELFfireOverlay))
  counts <- ELFfireCounts(rasWhole, points, NULL, fireYears, pixelAreaHa,
                          escapeSizeHa = escapeSizeHa, overlay = overlay)
  status <- ELFfitStatus(counts, minNaturalIgnitions = minNaturalIgnitions,
                         minFirePolygons = minFirePolygons, minEscapes = minEscapes)
  plan <- ELFmergePlan(status, ELFneighbours(rasWhole),
                       minNaturalIgnitions = minNaturalIgnitions,
                       minFirePolygons = minFirePolygons, minEscapes = minEscapes)
  if (!is.null(landCoverFile) && length(minRegionAreaKm2) == 1L && !is.na(minRegionAreaKm2)) {
    ## the size rule sees the regions that remain: fire-count merges applied, skipped ELFs left out
    remaining <- mergeELFs(ELFs, plan)$rasWhole
    keep <- setdiff(names(remaining), c(ELFsSkipped(plan), ELFsNotFitted(names(remaining))))
    remaining <- terra::rast(unname(remaining[keep]))
    names(remaining) <- keep
    merges <- plan[plan$action == "merge", ]
    members <- stats::setNames(lapply(keep, function(id) {
      hit <- match(id, merges$ELF)
      if (is.na(hit)) id else merges$members[[hit]]
    }), keep)
    regionStats <- ELFregionStats(remaining, terra::rast(landCoverFile), fireYears = fireYears,
                                  overlay = overlay, members = members)
    regionStats$members <- unname(members)
    sizePlan <- ELFsizePlan(regionStats, ELFneighbours(remaining), minAreaKm2 = minRegionAreaKm2,
                            maxLandCoverDist = maxLandCoverDist, maxBurnRatio = maxBurnRatio)
    plan <- combinePlans(plan, sizePlan, status)
  }
  for (i in seq_len(nrow(plan))) {
    message("fireSense_ELFs: ", plan$action[i], " ", paste(plan$members[[i]], collapse = " + "),
            if (plan$action[i] == "merge") paste0(" -> ", plan$ELF[i]), " (", plan$reason[i], ")")
  }
  list(ELFs = mergeELFs(ELFs, plan), status = status, plan = plan, excluded = ELFsSkipped(plan))
}

#' Add the size merges to the fire-count plan
#'
#' A size merge covers every ELF in its `members`, so a fire-count merge whose members are all in one is
#' replaced by it. Skips stay. The size merges get the fire counts of their members from `status`.
#'
#' @param firePlan Plan from [ELFmergePlan()].
#' @param sizePlan Plan from [ELFsizePlan()].
#' @param status [ELFfitStatus()] of the original ELFs.
#'
#' @return A `data.table` in `ELFmergePlan()`'s format.
#' @export
combinePlans <- function(firePlan, sizePlan, status) {
  if (!nrow(sizePlan)) return(firePlan)
  covered <- vapply(seq_len(nrow(firePlan)), function(i) {
    firePlan$action[i] == "merge" &&
      any(vapply(sizePlan$members, function(m) all(firePlan$members[[i]] %in% m), logical(1)))
  }, logical(1))
  total <- function(col) vapply(sizePlan$members, function(m) sum(status[[col]][match(m, status$ELF)]), status[[col]][1])
  sizeRows <- data.table::data.table(
    action = sizePlan$action, ELF = sizePlan$ELF, members = sizePlan$members,
    naturalIgnitions = total("naturalIgnitions"), escapes = total("escapes"), firePolygons = total("firePolygons"),
    reason = sizePlan$reason)
  data.table::rbindlist(list(firePlan[!covered, ], sizeRows))
}
