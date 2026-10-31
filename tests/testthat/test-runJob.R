## runJob() with SpaDES.core and SpaDES.project mocked: no simulation, no queue and no Drive is touched.
fakeSim <- function(phase = NULL, mode = "fit", heldOutFold = NA) {
  list(runName = "6.6.1/run1",
       params = list(.globals = list(heldOutFold = heldOutFold),
                     fireSense_spreadFit = c(list(mode = mode), if (!is.null(phase)) list(phase = phase))))
}

## `skip`: what queueSkipForFits() returns; `fail`: whether outSaveTarUpload() errors. `calls` records every call.
mockJob <- function(skip = NULL, fail = FALSE, env = parent.frame()) {
  calls <- new.env()
  calls$spades <- list()
  calls$skip <- list()
  calls$wait <- list()
  calls$upload <- list()
  testthat::local_mocked_bindings(
    simInitAndSpades2 = function(sim, ...) {
      calls$spades <- c(calls$spades, list(sim))
      "simOut"
    },
    .package = "SpaDES.core", .env = env)
  testthat::local_mocked_bindings(
    queueSkipForFits = function(fitQueue, id) {
      calls$skip <- c(calls$skip, list(list(fitQueue, id)))
      skip
    },
    queueWaitWhile = function(queueFile, id, ...) {
      calls$wait <- c(calls$wait, list(list(queueFile, id)))
      invisible(NULL)
    },
    outSaveTarUpload = function(runName, sim, gFolder) {
      calls$upload <- c(calls$upload, list(list(runName = runName, sim = sim, gFolder = gFolder)))
      if (fail) stop("boom")
      invisible(NULL)
    },
    .package = "SpaDES.project", .env = env)
  calls
}

test_that("runJob() passes the phase's events to SpaDES.core and returns the simList", {
  calls <- mockJob()
  sim <- fakeSim()
  expect_identical(runJob(sim, "6.6.1", phase = "inputs"), "simOut")
  expect_identical(calls$spades[[1]]$events, phaseEvents("inputs"))
  runJob(sim, "6.6.1", phase = "fit")
  expect_identical(calls$spades[[2]]$events, phaseEvents("fit"))
  runJob(sim, "6.6.1", phase = "all")
  expect_null(calls$spades[[3]]$events)
  expect_null(sim$events) # the sim is copied, not changed
})

test_that("runJob() takes the phase from the spreadFit parameter, and 'all' without one", {
  calls <- mockJob()
  runJob(fakeSim(phase = "fit"), "6.6.1")
  expect_identical(calls$spades[[1]]$events, phaseEvents("fit"))
  runJob(fakeSim(), "6.6.1")
  expect_null(calls$spades[[2]]$events)
  expect_error(runJob(fakeSim(), "6.6.1", phase = "predict"), "should be one of")
})

test_that("a fit that validates stops after cross-validation", {
  calls <- mockJob()
  runJob(fakeSim(mode = c("fit", "validate")), "6.6.1", phase = "fit")
  runJob(fakeSim(heldOutFold = 2L), "6.6.1", phase = "fit")
  runJob(fakeSim(mode = "fit", heldOutFold = NA), "6.6.1", phase = "fit")
  expect_identical(calls$spades[[1]]$events, phaseEvents("fit", validates = TRUE))
  expect_identical(calls$spades[[2]]$events, phaseEvents("fit", validates = TRUE))
  expect_identical(calls$spades[[3]]$events, phaseEvents("fit", validates = FALSE))
})

test_that("an inputs job skips an ELF the fit queue has, and says why", {
  calls <- mockJob(skip = "the fit queue has claimed it")
  expect_message(out <- runJob(fakeSim(), "6.6.1", phase = "inputs", fitQueue = "fits.rds"),
                 "the fit queue has claimed it")
  expect_null(out)
  expect_identical(calls$skip, list(list("fits.rds", "6.6.1")))
  expect_length(calls$spades, 0L) # simInit was never called
})

test_that("an inputs job runs when the fit queue has no claim, and only inputs jobs ask", {
  calls <- mockJob(skip = NULL)
  expect_identical(runJob(fakeSim(), "6.6.1", phase = "inputs", fitQueue = "fits.rds"), "simOut")
  runJob(fakeSim(), "6.6.1", phase = "inputs") # no queue
  runJob(fakeSim(), "6.6.1", phase = "fit", fitQueue = "fits.rds")
  runJob(fakeSim(), "6.6.1", phase = "all", fitQueue = "fits.rds")
  expect_length(calls$skip, 1L)
  expect_length(calls$spades, 4L)
})

test_that("a fit job waits on the inputs queue, and only fit jobs wait", {
  calls <- mockJob()
  runJob(fakeSim(), "6.6.1", phase = "fit", afterQueue = "inputs.rds")
  runJob(fakeSim(), "6.6.1", phase = "fit")
  runJob(fakeSim(), "6.6.1", phase = "inputs", afterQueue = "inputs.rds")
  runJob(fakeSim(), "6.6.1", phase = "all", afterQueue = "inputs.rds")
  expect_identical(calls$wait, list(list("inputs.rds", "6.6.1")))
  expect_length(calls$spades, 4L)
})

test_that("only a phase 'all' job that is not a hindcast uploads, and only with an uploadDir", {
  calls <- mockJob()
  expect_identical(runJob(fakeSim(), "6.6.1", phase = "all", uploadDir = "folderID"), "simOut")
  expect_identical(calls$upload,
                   list(list(runName = "6.6.1/run1", sim = "simOut", gFolder = "folderID")))
  runJob(fakeSim(), "6.6.1", phase = "all", uploadDir = "folderID", hindcast = TRUE)
  runJob(fakeSim(), "6.6.1", phase = "inputs", uploadDir = "folderID")
  runJob(fakeSim(), "6.6.1", phase = "fit", uploadDir = "folderID")
  runJob(fakeSim(), "6.6.1", phase = "all")
  expect_length(calls$upload, 1L)
})

test_that("an upload failure only warns", {
  calls <- mockJob(fail = TRUE)
  expect_warning(out <- runJob(fakeSim(), "6.6.1", phase = "all", uploadDir = "folderID"),
                 "outSaveTarUpload failed: boom")
  expect_identical(out, "simOut")
})

test_that("runJob() sources diagFile and profiles to rprofPrefix when they are given", {
  calls <- mockJob()
  diagFile <- withr::local_tempfile(fileext = ".R")
  writeLines("options(fireSenseUtils.test.diag = TRUE)", diagFile)
  withr::local_options(fireSenseUtils.test.diag = NULL)
  prefix <- withr::local_tempfile()
  runJob(fakeSim(), "6.6.1", phase = "inputs", diagFile = diagFile, rprofPrefix = prefix)
  expect_true(getOption("fireSenseUtils.test.diag"))
  expect_true(file.exists(paste0(prefix, "_6.6.1_", Sys.getpid(), ".out")))
})
