utils::globalVariables(c(
  ".N", "ids", "Year"
))

#' Outer wrapper on spread fire polygon data munging that does several things:
#' 1. ensure buffered fires are entirely in studyArea
#' 2. ensure every fire has a corresponding ignition point, and vice versa
#' 3. ensure these points are flammable
#'
#' @param firePolys the semi-processed fire polys, with field matching pointsIDcolumn
#' @template flammableRTM
#' @param spreadFirePoints the ignition points corresponding to `firePolys`
#' @param areaMultiplier Either a scalar that will buffer `areaMultiplier * fireSize`
#' or a function of `fireSize`. See `?fireSenseUtils::bufferToArea`.
#' @param minSize an alternative to areaMultiplier, typically used when fires are small
#' @param pointsIDcolumn the name of the column denoting fire ids in both
#' spreadFirePoints and firePolys
#'
#' @export
#' @importFrom purrr transpose
#' @importFrom data.table rbindlist
harmonizeFireData <- function(firePolys, flammableRTM, spreadFirePoints,
                              areaMultiplier, minSize, pointsIDcolumn = "FIRE_ID") {
  ## safety to ensure missing years actually removed
  spreadFirePoints[is.na(names(spreadFirePoints))] <- NULL

  fireYears <- names(firePolys)
  fireBufferedListDT <- bufferToArea(
    poly = firePolys,
    polyName = fireYears,
    rasterToMatch = flammableRTM,
    verb = TRUE,
    areaMultiplier = areaMultiplier,
    field = pointsIDcolumn,
    minSize = minSize
    )

  #TODO: is this old cache behavior? or a problem with bad data?
  if (!any(sapply(fireBufferedListDT, is.data.table))) {
    fireBufferedListDT <- lapply(fireBufferedListDT, as.data.table)
  }

  # remove polygons where the buffered fire area exceeds rasterToMatch
  # these are fires that are too close to the border, and therefore we lack the
  # vegetation information to characterize the burn properly, even if the fire is
  # inside the border. This is why studyArea must be buffered. 17-7-2024

  #removal happens outside of bufferToArea. It could be done inside,
  #but this way the removal of fires is transparent and can be quantified.
  #temporarily add Year to count fires and ensure we have correct year if entire year excluded
  fireBufferedListDT <- lapply(names(fireBufferedListDT), function(year){
    annual <- fireBufferedListDT[[year]]
    annual[, Year := year]
  })
  totalFires <- nrow(rbindlist(fireBufferedListDT)[, .N, .(Year, ids)])

  fireBufferedListDT <- lapply(fireBufferedListDT, FUN = removeBufferedFiresOutsideRTM,
                               flammableRTM = flammableRTM)
  newYears <- rbindlist(fireBufferedListDT)[, .N, .(Year)]

  firesRemaining = nrow(rbindlist(fireBufferedListDT)[, .N, .(Year, ids)])
  #messaging
  firePeriod <- paste0(gsub(x = fireYears, "year", "")[1], "-",
                       gsub(x = rev(fireYears), "year", "")[1])
  message(paste0("between ", firePeriod, " there were ", totalFires, " escaped fires",
                 " and ", totalFires-firesRemaining, " were replaced for bordering the studyArea border"))
  fireBufferedListDT <- lapply(fireBufferedListDT, function(x) {x[, Year := NULL]})
  #ensure year derived from name, not data.table itself
  names(fireBufferedListDT) <- newYears$Year

  harmonized <- harmonizeBufferAndPoints(
    cent = spreadFirePoints,
    buff = fireBufferedListDT,
    ras = flammableRTM,
    idCol = pointsIDcolumn
  )

  ## ensure mismatched (e.g. points w/ no polys) and now missing years actually removed.
  emptyYearsPoints <- which(vapply(harmonized, is.null, logical(1))) |> names()
  emptyYearsPolys <- which(vapply(fireBufferedListDT, function(x) nrow(x) == 0, logical(1))) |> names()
  stopifnot(emptyYearsPoints == emptyYearsPolys)
  emptyYears <- unique(c(emptyYearsPoints, emptyYearsPolys))
  if (length(emptyYears > 0)) {
    harmonized[[emptyYears]] <- NULL
    fireBufferedListDT[[emptyYears]] <- NULL
  }

  #make sure that cleanUpSpreadFirePoints removes points with no polygons
  harmonized <- Map(f = cleanUpSpreadFirePoints,
                    firePoints = harmonized,
                    bufferDT = fireBufferedListDT,
                    MoreArgs = list(flammableRTM = flammableRTM)) |>
    purrr::transpose()

  #one last safety check - ensuring points have corresponding fires
  nfires_poly <- sapply(harmonized$FireBuffered, FUN = function(x){nrow(x[, .N, .(ids)])})
  #polygons can have multiple ignition points, but the reverse is not true due to table structure
  nfires_point <- sapply(harmonized$SpatialPoints, FUN = function(x){length(unique(x$FIRE_ID))})
  if (!identical(nfires_poly, nfires_point)) {
    stop('spread fire point and poly harmonization error in dataPrepFit. Please debug harmonizeFireData')
  }

  #fire polys should be cleaned too - though at this point it stops being used.

  return(list(fireBufferedListDT = harmonized$FireBuffered,
              firePolys = firePolys, #should be returned because some years may have been converted to NULL
              spreadFirePoints = harmonized$SpatialPoints))
}
