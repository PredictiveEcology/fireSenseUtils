utils::globalVariables(c(
  "plotID"
))

#' Read lightning data
#'
#' Data is not publicly available; it is currently stored in a User Access Controlled location.
#' Data comes from William Burrows. Please contact him to obtain data.
#' This function will get the data (which is in a `.txt` format) and do all the necessary
#' sequence of calculations to create a `SpatRaster`.
#'
#' @param filename The local filename to read with `data.table::fread`
#'
#' @param ... Passed to `postProcess`, e.g., `to`, `maskTo`
#'
#' @returns This will return a `SpatRaster` object with resolution the same as that
#'   supplied by `to` or `projectTo`. If no `SpatRaster` is supplied, it will
#'   default to resolution of `10000` m (10 km x 10 km).
#'
#' @export
#' @importFrom data.table fread dcast
#' @importFrom terra res rast rasterize aggregate focal
#' @references Kochtubajda, B., & Burrows, W. R. (2020). Cloud-to-Ground Lightning in Canada:
#' 20 Years of CLDN Data. Atmosphere-Ocean, 58(5), 316–332.
#' <https://doi.org/10.1080/07055900.2020.1845117>
#'
#' @examples
#' \dontrun{
#' # The example fetches from a Google Drive folder that requires
#' # `drive_auth()` credentials. R CMD check runs examples
#' # non-interactively where no Drive token is available, so the body
#' # is gated under \dontrun{} to keep docs complete without breaking
#' # checks.
#' if (require(googledrive, quietly = TRUE)) {
#'   crsToUse <- paste("+proj=lcc +lat_0=0 +lon_0=-95 +lat_1=49 +lat_2=77",
#'                     "+x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs")
#'   ras <- terra::rast(
#'     terra::ext(-1115000, -376750, 7267000, 7874000),
#'     resolution = 250,
#'     vals = 1,
#'     crs = crsToUse
#'   )
#'   ## Google Drive folder:
#'   ## <https://drive.google.com/drive/folders/1ftpKKUO8BZcJ9ba4LDWcmF5uQkpRGmku>
#'   data <- list(
#'     lightningDays = "1jeKJquhVJsesoNk2EPP1QZkttX3Zwp5c",
#'     lightningDensity = "12fnhfKtER-JXkl06M4_yZ3GvpZWtQlIr",
#'     positiveCG = "1bn6cQ23tvPicFLHn1tz4Z3AqDzJI4r60",
#'     positiveCGdensity = "1GNixhXj1Ex1jT0tWXfhmxef-dX3ze1a4"
#'   )
#'   lightning <- Map(
#'     url = data,
#'     f = function(url) {
#'       reproducible::prepInputs(
#'         url = url,
#'         fun = readLightningData(targetFile, to = ras),
#'         destinationPath = file.path(tempdir(), "lightning_data")
#'       )
#'     }
#'   )
#' }
#' }
readLightningData <- function(filename, ...) {
  a <- data.table::fread(filename)
  ldName <- "LightningDensity"
  shortColnames <- c("Lat", "Long", ldName)
  re <- rep(seq(NCOL(a)/length(shortColnames)), each = length(shortColnames))
  cn <- paste0(shortColnames, re)
  colnames(a) <- cn
  b <- melt(a, measure.vars = c(grep("Lat", cn, value = TRUE), grep("Long", cn, value = TRUE),
                                grep(ldName, cn, value = TRUE)))
  re2 <- rep(seq(prod(dim(a))/length(shortColnames)), length(shortColnames))
  b[, variable := gsub("[0-9]", "", variable)]
  b[, plotID := re2]
  d <- data.table::dcast(b, plotID ~ variable, value.var = "value")
  d <- na.omit(d)
  ld = sf::st_as_sf(d[, -"plotID"], coords = c("Long", "Lat"),
                    crs = 4326, agr = "constant")
  dots <- list(...)
  res10k <- 1e4

  if (length(dots)) {
    ld <- postProcess(ld, ...)
    # if (length(to) > 1) {
    gridded <- mapply(x = dots, function(x) reproducible::.isGridded(x), SIMPLIFY = TRUE)
    # area <- sapply(dots, SpaDES.project:::areas, USE.NAMES = TRUE)
    # largest <- which.max(area)
    whToUse <- intersect(names(dots), c("projectTo", "to"))[1]
    if (all(!gridded)) {
      rtm <- Map(x = dots[whToUse], function(x) terra::rast(x, res = res10k))
    } else {
      rtm <- dots[whToUse]
    }
    rtm <- rtm[[1]]
    rtm10k <- terra::aggregate(rtm, fact = res10k/terra::res(rtm)[1])
    b <- terra::rasterize(ld, rtm10k, field = ldName)
    namesOrig <- names(ld)[1]
    ld <- terra::focal(b, 3, "mean", na.policy="only")
    names(ld)[1] <- namesOrig

    ld <- postProcess(ld, ...) # bring it back to resolution of to if supplied

  } else {
    dims <- dim(a)
    ldsv <- terra::vect(ld)
    crsToUse <- "+proj=lcc +lat_0=0 +lon_0=-95 +lat_1=49 +lat_2=77 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs"
    ld2 <- postProcess(ldsv, projectTo = crsToUse, field = ldName)
    ext <- round(terra::ext(ld2), digits = -4)
    rtm10k <- rast(ext, res = as.integer(res10k), vals = 1, crs = crsToUse)
    b <- terra::rasterize(ld2, rtm10k, field = ldName)
    ld <- terra::focal(b, 3, "mean", na.policy="only")
  }
  ld
}

