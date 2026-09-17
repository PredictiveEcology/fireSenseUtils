## bufferToArea() and rasterFireBufferDT() pick buffer pixels at random, in forked workers when cores > 1.
## Each forked child used to seed itself independently, so the same call under the same seed gave
## different buffers on every run, and different buffers again with one core. fireSense_dataPrepFit's
## spread-fit data then changed between runs of the same study area (61 of 5,000 buffer pixels differed in
## one year of ELF 11.2). Every element now runs under a seed drawn from the session's stream, forked or not.

skipUnlessFork <- function() {
  skip_on_cran()
  skip_on_os("windows") # no fork
  skip_if_not_installed("Require")
  skip_if(Require:::isRstudio(), "these functions do not fork in RStudio")
  skip_if(parallelly::availableCores(constraints = "connections", omit = 1) < 2, "needs 2 cores to fork")
}

## data.tables carry an external pointer that identical() compares by address
asFrames <- function(x) lapply(x, function(d) if (is.null(d)) NULL else as.data.frame(d))

test_that("bufferToArea gives the same buffers for the same seed, forked or not", {
  ## R CMD check caps availableCores() at 2; with omit = 1 nothing would fork
  withr::local_envvar(c("_R_CHECK_LIMIT_CORES_" = NA))
  skipUnlessFork()
  rtm <- terra::rast(nrows = 300, ncols = 300, extent = c(0, 120000, 0, 120000),
                     crs = "EPSG:3978", vals = 1L)
  withr::with_seed(123, {
    polys <- lapply(1:3, function(y) {
      pts <- sf::st_as_sf(data.frame(x = runif(6, 20000, 100000), y = runif(6, 20000, 100000),
                                     FIRE_ID = as.numeric(seq_len(6) + y * 10)),
                          coords = c("x", "y"), crs = 3978)
      sf::st_buffer(pts, dist = runif(6, 500, 2500))
    })
  })
  names(polys) <- paste0("year", 2001:2003)
  run <- function(seed, cores) {
    withr::with_seed(seed, asFrames(
      bufferToArea(poly = polys, rasterToMatch = rtm, areaMultiplier = 10, field = "FIRE_ID",
                   minSize = 500, cores = cores, verb = FALSE)))
  }

  serial <- run(1, cores = 1)
  expect_false(identical(serial, run(2, cores = 1))) # the seed matters, so the checks below can fail
  forked <- run(1, cores = 2)
  expect_identical(forked, run(1, cores = 2))
  expect_identical(forked, serial)
})

## A small study area with three fire years, shared by the rasterFireBufferDT tests
fireRasterFixture <- function() {
  flammable <- terra::rast(nrows = 300, ncols = 300, extent = c(0, 120000, 0, 120000),
                           crs = "EPSG:3978", vals = 1L)
  fireRaster <- terra::rast(flammable)
  terra::values(fireRaster) <- NA_integer_
  withr::with_seed(123, {
    for (yr in 2001:2003) for (k in 1:4) {
      r0 <- sample(20:260, 1)
      c0 <- sample(20:260, 1)
      fireRaster[r0:(r0 + sample(5:15, 1)), c0:(c0 + sample(5:15, 1))] <- yr
    }
  })
  list(fireRaster = fireRaster, flammable = flammable)
}

test_that("rasterFireBufferDT gives the same buffers for the same seed, forked or not", {
  withr::local_envvar(c("_R_CHECK_LIMIT_CORES_" = NA))
  skipUnlessFork()
  fx <- fireRasterFixture()
  run <- function(seed, cores) {
    asFrames(
      rasterFireBufferDT(years = 2001:2003, fireRaster = fx$fireRaster, flammableRTM = fx$flammable,
                         bufferForFireRaster = 1000, areaMultiplier = 10, minSize = 500,
                         verb = 0, cores = cores, seed = seed))
  }

  serial <- run(1, cores = 1)
  expect_false(identical(serial, run(2, cores = 1))) # the seed matters, so the checks below can fail
  forked <- run(1, cores = 2)
  expect_identical(forked, run(1, cores = 2))
  expect_identical(forked, serial)
})

## fireSense_dataPrepFit caches this call. Two runs with identical inputs used to draw their seeds from
## whatever state the session's RNG was in, so they built different buffers -- and different spread-fit
## covariates -- under cache keys that differed only by accident (SpaDES.core #456 would have merged them).
## By default the seed now comes from the inputs, so the same inputs give the same buffers in any session.
test_that("rasterFireBufferDT's default buffers depend on its inputs, not on the session's RNG", {
  fx <- fireRasterFixture()
  run <- function(ambientSeed, areaMultiplier = 10) {
    withr::with_seed(ambientSeed, asFrames(
      rasterFireBufferDT(years = 2001:2003, fireRaster = fx$fireRaster, flammableRTM = fx$flammable,
                         bufferForFireRaster = 1000, areaMultiplier = areaMultiplier, minSize = 500,
                         verb = 0, cores = 1)))
  }

  base <- run(1)
  expect_identical(run(2), base)
  expect_identical(run(3), base)
  ## ... while the inputs still matter
  expect_false(identical(run(1, areaMultiplier = 5), base))

  ## and the call leaves the session's random stream where it was
  set.seed(42)
  before <- .Random.seed
  run2 <- asFrames(rasterFireBufferDT(years = 2001:2003, fireRaster = fx$fireRaster,
                                      flammableRTM = fx$flammable, bufferForFireRaster = 1000,
                                      areaMultiplier = 10, minSize = 500, verb = 0, cores = 1))
  expect_true(identical(.Random.seed, before))
  expect_identical(run2, base)
})
