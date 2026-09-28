## harmonizeFireData() drops a fire year whose fires all lie outside the study area. The fire
## polygons of the years after it were then filtered with another year's fire IDs, and lost.

circleFire <- function(x, y, r, id) {
  sf::st_buffer(sf::st_sf(FIRE_ID = id, geometry = sf::st_sfc(sf::st_point(c(x, y)), crs = 3978)), r)
}
firePoint <- function(x, y, id) {
  sf::st_sf(FIRE_ID = id, geometry = sf::st_sfc(sf::st_point(c(x, y)), crs = 3978))
}

test_that("harmonizeFireData keeps each year's fires when an earlier year is dropped", {
  ## 40 x 40 km at 400 m; east of x = 36 km is outside the study area (NA)
  rtm <- terra::rast(nrows = 100, ncols = 100, extent = c(0, 40000, 0, 40000),
                     crs = "EPSG:3978", vals = 1L)
  rtm[terra::xFromCell(rtm, seq_len(terra::ncell(rtm))) > 36000] <- NA
  polys <- list(year2001 = circleFire(8000, 30000, 1500, 1),
                year2002 = circleFire(38500, 20000, 800, 2),   # all outside
                year2003 = circleFire(10000, 10000, 1500, 3))
  pts <- list(year2001 = firePoint(8000, 30000, 1),
              year2002 = firePoint(38500, 20000, 2),
              year2003 = firePoint(10000, 10000, 3))
  set.seed(1)

  ## the misaligned years also made Map() recycle, with a warning
  expect_no_warning(capture.output(
    res <- suppressMessages(harmonizeFireData(polys, rtm, pts, areaMultiplier = 5, minSize = 100))))

  expect_named(res$firePolys, c("year2001", "year2003"))
  expect_identical(res$firePolys$year2003$FIRE_ID, 3)
  expect_identical(res$firePolys$year2001$FIRE_ID, 1)
  expect_identical(res$spreadFirePoints$year2003$FIRE_ID, 3)
})

## Every package function harmonizeFireData() reaches, directly or through the functions it calls
## (including S3 methods and functions passed as arguments), must be in harmonizeFireDataDeps();
## otherwise a cached call misses its change.
test_that("harmonizeFireDataDeps() lists every package function harmonizeFireData() reaches", {
  ns <- asNamespace("fireSenseUtils")
  own <- Filter(function(o) is.function(get(o, ns)), ls(ns, all.names = TRUE))
  reached <- character()
  todo <- "harmonizeFireData"
  while (length(todo)) {
    f <- todo[1]
    todo <- todo[-1]
    if (f %in% reached) next
    reached <- c(reached, f)
    fn <- get(f, ns)
    calls <- intersect(codetools::findGlobals(fn), own)
    if (any(grepl("UseMethod", deparse(body(fn)))))
      calls <- c(calls, grep(paste0("^", f, "\\."), own, value = TRUE))
    todo <- c(todo, setdiff(calls, reached))
  }
  deps <- harmonizeFireDataDeps()
  inDeps <- vapply(setdiff(reached, "harmonizeFireData"), function(f)
    any(vapply(deps, identical, logical(1), get(f, ns))), logical(1))
  expect_true(all(inDeps), info = paste("missing:", paste(names(inDeps)[!inDeps], collapse = ", ")))
})
