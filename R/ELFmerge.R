## Merging ELFs that have too few fires.
##
## An ELF with too few fires over the whole fire record cannot be fitted reliably (ELFfitStatus()).
## Eliot, 2026-09-14: merge it with a neighbour that shares its base -- a piece of a split ecoprovince
## (3.1.2) with another piece of the same province (3.1.1), a whole ecoprovince (12.1) with another of
## the same ecozone (12.2) -- taking the one that shares the longest border. If the pair still has too
## few fires, neither is fitted. A merged ELF is named by the shared base and the last parts of its
## members: 3.2.1 with 3.2.4 is "3.2.1_4", 12.1 with 12.2 is "12.1_2".

.ELFparent <- function(ids) sub("\\.[^.]*$", "", ids)
.ELFdepth <- function(ids) lengths(strsplit(ids, ".", fixed = TRUE))
.ELFlast <- function(ids) sub("^.*\\.", "", ids)

## Ecozones 1 and 2 are out of fireSense permanently (Eliot, 2026-09-14: "Elf 1.xx and 2.xx are out.
## Omit permanently"). runELFs() never queues them, and they are neither counted nor merged.
.arcticELFsPattern <- "^1\\.|^2\\."

#' ELFs permanently left out of fireSense
#'
#' ELFs of ecozones 1 and 2 (ids `1.*` and `2.*`) are never fitted: [runELFs()] leaves them out of
#' every list and map, and [ELFmergePlan()] neither merges them nor uses them as partners.
#'
#' @param ids ELF ids.
#'
#' @return The ids in `ids` that are left out.
#'
#' @export
ELFsArctic <- function(ids) {
  ids[grepl(.arcticELFsPattern, ids)]
}

#' Shared core borders between ELFs
#'
#' Measures how much core border each pair of ELFs shares on the common grid of `rasWhole`: the
#' number of horizontally or vertically adjacent cell pairs where one cell is core (`2`) in one
#' ELF and the other is core in the other.
#'
#' @param rasWhole ELF layers on one common grid, as a multi-layer `SpatRaster` or a named list of
#'   `SpatRaster`s (`sim$ELFs$rasWhole`). Cells are 0 (outside), 1 (buffer) or 2 (core).
#'
#' @return A `data.table` with one row per pair of ELFs whose cores touch: `ELF1`, `ELF2`,
#'   `sharedEdges` (cell edges) and `sharedLength` (in the grid's map units).
#'
#' @export
#' @importFrom data.table data.table rbindlist
#' @importFrom terra ncell rast res values
ELFneighbours <- function(rasWhole) {
  if (is.list(rasWhole)) rasWhole <- terra::rast(unname(rasWhole))
  stopifnot(inherits(rasWhole, "SpatRaster"), !is.null(names(rasWhole)))
  ids <- names(rasWhole)

  ## The ELF whose core each cell is in; 0 for none. Cores do not overlap; if they did, the first wins.
  owner <- integer(terra::ncell(rasWhole))
  for (k in seq_along(ids)) {
    v <- terra::values(rasWhole[[k]], mat = FALSE)
    owner[owner == 0L & v %in% 2] <- k
  }
  m <- matrix(owner, nrow = terra::nrow(rasWhole), ncol = terra::ncol(rasWhole), byrow = TRUE)

  pairsOf <- function(a, b) {
    keep <- a > 0L & b > 0L & a != b
    data.table::data.table(i = pmin(a[keep], b[keep]), j = pmax(a[keep], b[keep]))
  }
  edges <- data.table::rbindlist(list(
    pairsOf(m[, -ncol(m), drop = FALSE], m[, -1L, drop = FALSE]),
    pairsOf(m[-nrow(m), , drop = FALSE], m[-1L, , drop = FALSE])
  ))
  counts <- edges[, list(sharedEdges = .N), by = c("i", "j")]

  data.table::data.table(
    ELF1 = ids[counts$i],
    ELF2 = ids[counts$j],
    sharedEdges = counts$sharedEdges,
    sharedLength = counts$sharedEdges * terra::res(rasWhole)[1]
  )
}

#' Name of a merged ELF
#'
#' @param members ids of the ELFs being merged; they must share one base.
#'
#' @return The shared base followed by the members' last parts in numeric order, e.g.
#'   `c("3.2.4", "3.2.1")` gives `"3.2.1_4"` and `c("12.2", "12.1")` gives `"12.1_2"`.
#'
#' @export
ELFmergedName <- function(members) {
  stopifnot(length(members) >= 2L, length(unique(.ELFparent(members))) == 1L)
  members <- members[order(as.numeric(.ELFlast(members)))]
  paste0(.ELFparent(members[1]), ".", paste(.ELFlast(members), collapse = "_"))
}

#' Decide which ELFs with too few fires to merge, and which not to fit
#'
#' Every ELF that [ELFfitStatus()] calls `"zero"` or `"few"` is considered once, in ELF order. Its
#' partner is the neighbouring ELF with the same base and the same depth -- another piece of the
#' same split ecoprovince, or another whole ecoprovince of the same ecozone -- that shares the
#' longest core border (ties: more fires, then ELF order) and is not already part of a merge. If
#' the two together reach both thresholds they are merged; otherwise neither is fitted. An ELF with
#' no such neighbour is not fitted.
#'
#' @param status `data.table` from [ELFfitStatus()].
#' @param neighbours `data.table` from [ELFneighbours()].
#' @param minNaturalIgnitions,minFirePolygons thresholds a merged ELF must reach, as in
#'   [ELFfitStatus()].
#'
#' @return A `data.table` with one row per decision: `action` (`"merge"` or `"skip"`), `ELF` (the
#'   merged ELF's name, `NA` for a skip), `members` (list of the ELF ids involved),
#'   `naturalIgnitions` and `firePolygons` (their totals) and `reason`.
#'
#' @export
#' @importFrom data.table data.table rbindlist
ELFmergePlan <- function(status, neighbours, minNaturalIgnitions = 50, minFirePolygons = 50) {
  stopifnot(
    is.data.frame(status),
    all(c("ELF", "naturalIgnitions", "firePolygons", "status") %in% names(status)),
    is.data.frame(neighbours),
    all(c("ELF1", "ELF2", "sharedLength") %in% names(neighbours))
  )
  ids <- as.character(status$ELF)
  ign <- stats::setNames(status$naturalIgnitions, ids)
  polys <- stats::setNames(status$firePolygons, ids)
  thin <- ids[status$status %in% c("zero", "few") & !ids %in% ELFsArctic(ids)]
  thin <- thin[order(numeric_version(thin))]

  decision <- function(action, members, reason) {
    data.table::data.table(
      action = action,
      ELF = if (action == "merge") ELFmergedName(members) else NA_character_,
      members = list(members[order(numeric_version(members))]),
      naturalIgnitions = sum(ign[members]),
      firePolygons = sum(polys[members]),
      reason = reason
    )
  }

  used <- character(0)
  out <- list()
  for (elf in thin) {
    if (elf %in% used) next
    touching <- neighbours[neighbours$ELF1 == elf | neighbours$ELF2 == elf, ]
    other <- ifelse(touching$ELF1 == elf, touching$ELF2, touching$ELF1)
    ok <- other %in% ids & !other %in% used &
      .ELFparent(other) == .ELFparent(elf) & .ELFdepth(other) == .ELFdepth(elf)
    if (!any(ok)) {
      out[[length(out) + 1L]] <- decision("skip", elf, "too few fires; no neighbour shares its base")
      used <- c(used, elf)
      next
    }
    other <- other[ok]
    shared <- touching$sharedLength[ok]
    partner <- other[order(-shared, -(ign[other] + polys[other]), numeric_version(other))][1]
    members <- c(elf, partner)
    enough <- sum(ign[members]) >= minNaturalIgnitions && sum(polys[members]) >= minFirePolygons
    out[[length(out) + 1L]] <- if (enough) {
      decision("merge", members, paste0("too few fires; merged with ", partner, ", the longest shared border"))
    } else {
      decision("skip", members, paste0("too few fires even together with ", partner))
    }
    used <- c(used, members)
  }

  if (!length(out)) {
    return(data.table::data.table(action = character(0), ELF = character(0), members = list(),
                                  naturalIgnitions = integer(0), firePolygons = integer(0),
                                  reason = character(0)))
  }
  data.table::rbindlist(out)
}

#' ELFs not fitted after merging
#'
#' @param plan `data.table` from [ELFmergePlan()].
#'
#' @return `character` vector of the ELF ids in the plan's skip decisions.
#'
#' @export
ELFsSkipped <- function(plan) {
  as.character(unlist(plan$members[plan$action == "skip"]))
}

#' The ELF to run for an ELF id
#'
#' @param ELF an ELF id.
#' @param plan `data.table` from [ELFmergePlan()].
#'
#' @return The merged ELF's name if `ELF` was merged, otherwise `ELF`.
#'
#' @export
ELFrunName <- function(ELF, plan) {
  hit <- which(plan$action == "merge" & vapply(plan$members, function(m) ELF %in% m, logical(1)))
  if (length(hit)) plan$ELF[hit[1]] else ELF
}

#' Merge ELF maps as planned
#'
#' Replaces the members of each merge in [ELFmergePlan()] by one ELF. Its `rasWhole` layer is the
#' cell-by-cell maximum of the members' layers, so a cell that is core in any member is core. A
#' buffer drawn around the union of cores equals the union of the members' buffers, so the buffer
#' does not need redrawing. Its `rasCentered` layer is that map projected onto a Lambert conformal
#' conic centred on the merged core, as [makeELFs()] centres each ELF. Its `poly` rows are rebuilt
#' from the merged layer the way [makeELFs()] builds them.
#'
#' @param ELFs list with `rasWhole` and `rasCentered` (named lists of `SpatRaster`) and optionally
#'   `poly`, as returned by [makeELFs()].
#' @param plan `data.table` from [ELFmergePlan()].
#'
#' @return `ELFs`, with each merge's members replaced by the merged ELF.
#'
#' @export
#' @importFrom terra as.polygons rast
mergeELFs <- function(ELFs, plan) {
  merges <- plan[plan$action == "merge", ]
  for (k in seq_len(nrow(merges))) {
    members <- merges$members[[k]]
    nam <- merges$ELF[k]
    whole <- max(terra::rast(unname(ELFs$rasWhole[members])))
    names(whole) <- nam

    ELFs$rasWhole[members] <- NULL
    ELFs$rasWhole[[nam]] <- whole
    ELFs$rasCentered[members] <- NULL
    ELFs$rasCentered[[nam]] <- .centredELF(whole)

    if (!is.null(ELFs$poly)) {
      a <- whole
      a[a[] == 0] <- NA
      vec <- terra::as.polygons(a)
      vec[, "ID"] <- nam
      vec[, "buffer"] <- vec[, nam]
      vec[, nam] <- NULL
      ELFs$poly <- rbind(ELFs$poly[!ELFs$poly$ID %in% members, ], vec)
    }
  }
  ELFs
}

## A merged ELF on a Lambert conformal conic centred on its core, the projection bufferOut() gives
## each ELF, at the same 5 km resolution.
#' @importFrom terra as.polygons classify ext ifel project trim
.centredELF <- function(whole, resolution = 5000) {
  core <- terra::as.polygons(terra::ifel(whole == 2, 1L, NA), dissolve = TRUE)
  exts <- round(terra::ext(terra::project(core, "epsg:4943")))[]
  middle <- mean(c(exts[["xmin"]], exts[["xmax"]]))
  prj <- paste0("+proj=lcc +lat_0=", exts[["ymin"]], " +lon_0=", middle, " +lat_1=", exts[["ymax"]],
                " +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs")
  out <- terra::project(terra::classify(whole, cbind(0, NA)), prj, res = resolution, method = "near")
  out <- terra::trim(out)
  names(out) <- names(whole)
  out
}
