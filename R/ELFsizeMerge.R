## Merging regions that are too small, with a similar neighbour.
##
## Eliot, 2026-10-10: a region whose core is smaller than `minAreaKm2` merges with a neighbour of the same
## group whose land cover and burn rate are similar. Candidates are judged by what the merge leaves behind:
## the one that leaves the fewest small regions without any partner wins, so the most similar pair is
## passed over when merging it would strand a third region (11.2 + 11.3 would strand 11.1). The rule knows
## nothing about ELFs except the default grouping (`.ELFparent()`): it takes a table of regions and a table of
## neighbours.

#' Area, land cover and burn rate of each region
#'
#' Statistics for [ELFsizePlan()]. A region is the core (cells equal to 2) of one layer of `rasWhole`.
#' Land-cover shares come from a regular grid of up to `maxSample` points in the core, so the result is
#' deterministic. The burn rate is the area of `firePolys` inside the core, as a percentage of the core per year.
#'
#' @param rasWhole `SpatRaster` of region layers on one grid, as in [ELFfireCounts()]. Cells are 0, 1 or 2
#'   (core). Its map units must be metres.
#' @param landCover `SpatRaster` of land-cover class codes.
#' @param firePolys `SpatVector` of fire perimeters with column `YEAR`, or `NULL` (burn rate `NA`).
#' @param fireYears years over which `firePolys` are counted; their number is the divisor of the burn rate.
#' @param classes land-cover codes whose shares are returned. The default is SCANFI's: water, rock, bryoids,
#'   shrubs, herbs, coniferous, broadleaf, mixedwood.
#' @param maxSample largest number of sample points per region.
#'
#' @return A `data.table` with one row per region: `ELF`, `areaKm2`, `burnRate` (%/yr) and `landCover`
#'   (list column of named numeric vectors of shares, names `classes`).
#'
#' @export
#' @importFrom data.table data.table rbindlist
#' @importFrom terra as.polygons classify crs expanse extract intersect project spatSample
ELFregionStats <- function(rasWhole, landCover, firePolys = NULL, fireYears = NULL,
                           classes = c(20, 30, 40, 50, 100, 210, 220, 230), maxSample = 20000L) {
  stopifnot(inherits(rasWhole, "SpatRaster"), inherits(landCover, "SpatRaster"))
  if (!is.null(firePolys)) {
    stopifnot(!is.null(fireYears))
    fireYears <- sort(unique(as.integer(fireYears)))
    firePolys <- terra::project(firePolys, terra::crs(rasWhole))
    firePolys <- firePolys[as.integer(firePolys$YEAR) %in% fireYears, ]
  }
  rows <- lapply(names(rasWhole), function(elf) {
    core <- terra::as.polygons(terra::classify(rasWhole[[elf]], cbind(c(0, 1), NA)), dissolve = TRUE)
    areaKm2 <- sum(terra::expanse(core, unit = "km"))
    pts <- terra::spatSample(core, maxSample, method = "regular")
    lcc <- terra::extract(landCover, terra::project(pts, terra::crs(landCover)), ID = FALSE)[[1]]
    shares <- tabulate(match(lcc, classes), length(classes))
    shares <- stats::setNames(shares / max(sum(shares), 1), classes)
    burn <- NA_real_
    if (!is.null(firePolys)) {
      burnedKm2 <- 0
      if (nrow(firePolys)) {
        clipped <- suppressWarnings(terra::intersect(firePolys, core)) # warns when nothing overlaps
        if (nrow(clipped)) burnedKm2 <- sum(terra::expanse(clipped, unit = "km"))
      }
      burn <- 100 * burnedKm2 / areaKm2 / length(fireYears)
    }
    data.table::data.table(ELF = elf, areaKm2 = areaKm2, burnRate = burn, landCover = list(shares))
  })
  data.table::rbindlist(rows)
}

## Bray-Curtis distance between two share vectors
.landCoverDist <- function(p, q) sum(abs(p - q)) / sum(p + q)

## Larger over smaller burn rate; a rate of 0 is taken as 1e-6, a missing rate as no difference
.burnRatio <- function(p, q) {
  if (is.na(p) || is.na(q)) return(1)
  max(p, q) / max(min(p, q), 1e-6)
}

## Candidate merges: pairs of neighbours of one group, at least one of them small, similar enough.
.sizeCandidates <- function(R, N, minAreaKm2, maxLandCoverDist, maxBurnRatio, group) {
  small <- R$id[R$area < minAreaKm2]
  empty <- data.table::data.table(a = character(0), b = character(0), lcDist = numeric(0), burnRatio = numeric(0))
  if (!length(small) || !nrow(N)) return(empty)
  N <- N[(N$a %in% small | N$b %in% small) & group(N$a) == group(N$b), ]
  if (!nrow(N)) return(empty)
  at <- function(id) match(id, R$id)
  N$lcDist <- mapply(function(i, j) .landCoverDist(R$landCover[[at(i)]], R$landCover[[at(j)]]), N$a, N$b)
  N$burnRatio <- mapply(function(i, j) .burnRatio(R$burn[at(i)], R$burn[at(j)]), N$a, N$b)
  N[N$lcDist <= maxLandCoverDist & N$burnRatio <= maxBurnRatio, ]
}

## The regions and neighbour table after merging regions i and j: area summed, land cover and burn rate
## weighted by area, edges combined. The new region is named from all the original members.
.sizeMerge <- function(R, N, i, j) {
  ri <- match(i, R$id); rj <- match(j, R$id)
  w <- R$area[c(ri, rj)] / sum(R$area[c(ri, rj)])
  members <- c(R$members[[ri]], R$members[[rj]])
  members <- members[order(numeric_version(members))]
  id <- ELFmergedName(members)
  burn <- if (anyNA(R$burn[c(ri, rj)])) NA_real_ else sum(w * R$burn[c(ri, rj)])
  new <- data.table::data.table(id = id, area = sum(R$area[c(ri, rj)]), burn = burn,
                                landCover = list(w[1] * R$landCover[[ri]] + w[2] * R$landCover[[rj]]),
                                members = list(members), why = list(character(0)))
  R2 <- rbind(R[-c(ri, rj), ], new)
  N$a[N$a %in% c(i, j)] <- id
  N$b[N$b %in% c(i, j)] <- id
  N <- N[N$a != N$b, ]
  list(R = R2, N = unique(data.table::data.table(a = pmin(N$a, N$b), b = pmax(N$a, N$b))), id = id)
}

#' Merge regions that are too small with a similar neighbour
#'
#' A region with `areaKm2 < minAreaKm2` is small. Its candidate partners are the neighbours of the same
#' `group` whose land-cover (Bray-Curtis) distance is at most `maxLandCoverDist` and whose burn-rate ratio
#' (larger over smaller) is at most `maxBurnRatio`. For each candidate the merge is simulated (areas summed,
#' land cover and burn rate weighted by area, neighbours combined) and the small regions left with no
#' candidate partner, the orphans, are counted. The candidate with the fewest orphans is accepted, ties
#' going to the smaller land-cover distance and then the smaller burn ratio. This repeats until no candidate
#' is left, so a merged region that is still small can merge again and a group can grow beyond two regions.
#'
#' @param stats `data.table` from [ELFregionStats()]: `ELF`, `areaKm2`, `burnRate`, `landCover`. An optional
#'   list column `members` names the original ids a region stands for (default: its own id), so a region
#'   that is already a merge keeps a name [ELFmergedName()] can build.
#' @param neighbours `data.table` from [ELFneighbours()], with `ELF1` and `ELF2`.
#' @param minAreaKm2 regions with a smaller core area are small. `NULL` or `NA` turns the rule off.
#' @param maxLandCoverDist,maxBurnRatio the similarity limits above.
#' @param group function of the region ids giving the group within which regions may merge. The default
#'   is the ELF base ([ELFmergePlan()]'s rule).
#'
#' @return A `data.table` with one row per merged region, as [ELFmergePlan()] returns for merges: `action`
#'   (`"merge"`), `ELF` (the merged name), `members` (list of the original ids) and `reason`, plus `areaKm2`.
#'   [mergeELFs()] applies it. Zero rows if nothing merges.
#'
#' @export
#' @importFrom data.table data.table rbindlist
ELFsizePlan <- function(stats, neighbours, minAreaKm2 = 35000, maxLandCoverDist = 0.35,
                        maxBurnRatio = 6, group = .ELFparent) {
  stopifnot(
    is.data.frame(stats), all(c("ELF", "areaKm2", "burnRate", "landCover") %in% names(stats)),
    is.data.frame(neighbours), all(c("ELF1", "ELF2") %in% names(neighbours))
  )
  noPlan <- data.table::data.table(action = character(0), ELF = character(0), members = list(),
                                   reason = character(0), areaKm2 = numeric(0))
  if (length(minAreaKm2) != 1L || is.na(minAreaKm2) || !nrow(stats)) return(noPlan)

  ids <- as.character(stats$ELF)
  R <- data.table::data.table(
    id = ids, area = stats$areaKm2, burn = stats$burnRate, landCover = stats$landCover,
    members = if ("members" %in% names(stats)) stats$members else as.list(ids),
    why = rep(list(character(0)), length(ids)))
  keep <- as.character(neighbours$ELF1) %in% ids & as.character(neighbours$ELF2) %in% ids
  N <- data.table::data.table(a = pmin(as.character(neighbours$ELF1), as.character(neighbours$ELF2))[keep],
                              b = pmax(as.character(neighbours$ELF1), as.character(neighbours$ELF2))[keep])
  N <- unique(N)

  repeat {
    C <- .sizeCandidates(R, N, minAreaKm2, maxLandCoverDist, maxBurnRatio, group)
    if (!nrow(C)) break
    C$orphans <- mapply(function(i, j) {
      m <- .sizeMerge(R, N, i, j)
      small <- m$R$id[m$R$area < minAreaKm2]
      cc <- .sizeCandidates(m$R, m$N, minAreaKm2, maxLandCoverDist, maxBurnRatio, group)
      sum(!small %in% c(cc$a, cc$b))
    }, C$a, C$b)
    o <- order(C$orphans, C$lcDist, C$burnRatio, numeric_version(sub("_.*$", "", C$a)),
               numeric_version(sub("_.*$", "", C$b)))
    C <- C[o, ]
    i <- C$a[1]; j <- C$b[1]
    ## the reason is told from the smaller region's side: it merged with `other`
    ai <- R$area[match(i, R$id)]; aj <- R$area[match(j, R$id)]
    other <- if (ai <= aj) j else i
    why <- c(R$why[[match(i, R$id)]], R$why[[match(j, R$id)]],
             sprintf("smaller than %s km2; merged with %s (land cover distance %.2f, burn ratio %.1f)",
                     format(minAreaKm2, scientific = FALSE), other, C$lcDist[1], C$burnRatio[1]))
    m <- .sizeMerge(R, N, i, j)
    R <- m$R; N <- m$N
    R$why[[match(m$id, R$id)]] <- why
  }

  merged <- R[lengths(R$members) > 1L & vapply(R$why, length, 1L) > 0L, ]
  if (!nrow(merged)) return(noPlan)
  o <- order(numeric_version(vapply(merged$members, `[`, "", 1L)))
  merged <- merged[o, ]
  data.table::data.table(action = "merge", ELF = merged$id, members = merged$members,
                         reason = vapply(merged$why, paste, "", collapse = "; "), areaKm2 = merged$area)
}
