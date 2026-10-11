# Extracted from test-ELFs.R:234

# prequel ----------------------------------------------------------------------
sq <- function(x0, x1, y0 = 0, y1 = 3e5)
  sprintf("POLYGON ((%1$s %3$s, %2$s %3$s, %2$s %4$s, %1$s %4$s, %1$s %3$s))", x0, x1, y0, y1)
makeTestEcos <- function() {
  provs <- c("4.1", "4.2", "4.3")
  withIsland <- paste0("MULTIPOLYGON (((0 0, 300000 0, 300000 300000, 0 300000, 0 0)), ",
                       "((400000 310000, 420000 310000, 420000 330000, 400000 330000, 400000 310000)))")
  ecos <- terra::vect(c(withIsland, sq(3e5, 6e5), sq(6e5, 12e5)), crs = "EPSG:3978")
  ecos$ECOZONE <- "4"
  ecos$ECOREGION <- ecos$ECOPROVINC <- ecos$ECODISTRIC <- provs
  ecos
}
makeTestFRU <- function(ecos) {
  x <- terra::rasterize(ecos, terra::rast(ecos, resolution = 5000), field = "ECOPROVINC")
  x <- as.numeric(x) + 1
  names(x) <- "FRU"
  x
}
localELFsDir <- function(env = parent.frame()) {
  dp <- withr::local_tempdir(tmpdir = dirname(tempdir()), .local_envir = env)
  withr::local_options(reproducible.useCache = FALSE,
                       reproducible.cachePath = file.path(dp, "cache"), .local_envir = env)
  dp
}
cellAt <- function(r, x, y) terra::extract(r, cbind(x, y))[[1]]
makeTestELFRasters <- function() {
  r <- terra::rast(nrows = 4, ncols = 10, xmin = 0, xmax = 10000, ymin = 0, ymax = 4000,
                   crs = "EPSG:3978")
  a <- terra::setValues(r, rep(c(2, 2, 2, 2, 1, 1, 0, 0, 0, 0), 4))
  b <- terra::setValues(r, rep(c(0, 0, 0, 0, 1, 1, 2, 2, 2, 2), 4))
  list(rasWhole = list(A = a, B = b))
}
overlapELFs <- function() {
  r <- terra::rast(nrows = 1, ncols = 12, xmin = 0, xmax = 12000, ymin = 0, ymax = 1000,
                   crs = "EPSG:3978")
  a <- terra::setValues(r, c(2, 2, 2, 1, 1, 1, 1, 1, 0, 0, 0, 0))
  b <- terra::setValues(r, c(0, 0, 0, 0, 1, 1, 1, 1, 1, 2, 2, 2))
  poly <- terra::vect(c(sq(0, 8000, 0, 1000), sq(4000, 12000, 0, 1000)), crs = "EPSG:3978")
  poly$ID <- c("A", "B")
  list(ras = list(rasWhole = list(A = a, B = b)), poly = poly)
}
labelsOf <- function(r) as.character(terra::as.data.frame(r, na.rm = FALSE)[[1]])
localRunELFs <- function(userName, excluded = NULL, env = parent.frame()) {
  sim <- SpaDES.core::simInit()
  rasList <- setNames(as.list(1:4), c("1.1", "4.1", "4.2", "6.1"))
  sim$ELFs <- list(rasCentered = rasList, rasWhole = rasList)
  sim$spreadFitPreRun <- setNames(data.frame(c("2.1", "4.1", "4.2")), polygonIDTxt)
  sim$ELFsExcluded <- excluded
  SpaDES.core::outputs(sim) <- data.frame(objectName = c("ELFs", "ELFsPostMerge", "studyAreas"),
                                          file = c("ELF_polygons.rds", "ELFs_postMerge.rds", "studyAreas4.rds"))
  sim@outputs$file <- sub("\\.rds$", ".png", sim@outputs$file)  # the pngs Plots() saves
  uploaded <- new.env()
  local_mocked_bindings(Cache = function(x, ...) x, cacheId = function(x) "id", .env = env)
  local_mocked_bindings(simInitAndSpades2 = function(l) sim, .package = "SpaDES.core", .env = env)
  local_mocked_bindings(user = function() userName, .package = "SpaDES.project", .env = env)
  local_mocked_bindings(drive_update = function(file, media) uploaded$media <- media,
                        .package = "googledrive", .env = env)
  local_mocked_bindings(getRemoteMetadata = function(...) list(remoteHash = "h"),
                        .package = "reproducible", .env = env)
  uploaded
}

# test -------------------------------------------------------------------------
prj <- list(modules = "fireSense_ELFs", paths = list(modulePath = withr::local_tempdir()),
              params = list(fireSense_ELFs = list()))
uploaded <- localRunELFs("emcintir")
expect_identical(runELFs(prj), c("4.1", "4.2"))
