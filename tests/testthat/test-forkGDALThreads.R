## GDAL keeps one process-wide worker-thread pool, created at the first multi-threaded raster write.
## A child forked after that inherited the pool but none of its threads, so bufferToArea()'s forked
## workers waited for them forever with 0 CPU (fits jobs stalled at prepSpreadFitData). The workers
## must finish.
test_that("bufferToArea forked workers do not hang after a multi-threaded write in the parent", {
  skip_on_os("windows") # no fork
  skip_if_not_installed("Require")
  skip_if(Require:::isRstudio(), "bufferToArea does not fork in RStudio")
  skip_if(parallelly::availableCores(constraints = "connections", omit = 1) < 2,
          "needs 2 cores to fork")

  old <- terra::terraOptions(print = FALSE)
  withr::defer(terra::terraOptions(todisk = old$todisk, threads = old$threads, tempdir = old$tempdir))
  terra::terraOptions(todisk = TRUE, threads = 4, tempdir = withr::local_tempdir())

  rtm <- terra::rast(nrows = 500, ncols = 500, extent = c(0, 120000, 0, 120000),
                     crs = "EPSG:3978", vals = 1L)
  terra::writeRaster(rtm, file.path(terra::terraOptions(print = FALSE)$tempdir, "rtm.tif"))

  set.seed(1)
  polys <- lapply(1:2, function(y) {
    pts <- sf::st_as_sf(data.frame(x = runif(5, 20000, 100000), y = runif(5, 20000, 100000),
                                   FIRE_ID = as.numeric(seq_len(5) + y * 10)),
                        coords = c("x", "y"), crs = 3978)
    sf::st_buffer(pts, dist = runif(5, 500, 2000))
  })
  names(polys) <- paste0("year", 2001:2002)

  ## in a child, so a hang fails the test instead of stalling the suite
  job <- parallel::mcparallel(
    bufferToArea(poly = polys, rasterToMatch = rtm, areaMultiplier = 10, field = "FIRE_ID",
                 minSize = 500, cores = 2, verb = FALSE))
  res <- NULL
  deadline <- Sys.time() + 120
  while (is.null(res) && Sys.time() < deadline)
    res <- parallel::mccollect(job, wait = FALSE, timeout = 1)
  if (is.null(res)) {
    system2("pkill", c("-KILL", "-P", job$pid))
    tools::pskill(job$pid, tools::SIGKILL)
    parallel::mccollect(job)
    fail("bufferToArea's forked workers did not finish within 120 s")
  } else {
    out <- res[[1]]
    expect_false(inherits(out, "try-error"), label = paste(out))
    expect_named(out, names(polys))
  }
})
