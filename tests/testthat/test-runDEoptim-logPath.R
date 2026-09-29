## The default log file must not carry a module's name: the modules are being renamed
## (fireSense_SpreadFit -> fireSense_spreadFit) and the package cannot know which one called it.

test_that("runDEoptim()'s default logPath names no module", {
  logPath <- eval(formals(runDEoptim)$logPath)
  expect_false(grepl("fireSense_", basename(logPath), ignore.case = TRUE))
  expect_match(basename(logPath), "^runDEoptim_.*\\.log$")
})

test_that("a caller's logPath is used as given", {
  form <- ~ 0 + CMD
  lower <- stats::setNames(c(0.25, 0.1, 0, 0), c("maxAsymptote", "inflectionPoint1", "CMD", "yearSpreadSD"))
  got <- NULL
  testthat::local_mocked_bindings(
    clusterSetup = function(..., logPath) {
      got <<- logPath
      list(itermax = 5, trace = FALSE, strategy = 2L, NP = 40L)
    },
    DEoptimIterative = function(fn, lower, upper, control, ...) list(),
    .package = "clusters")
  withr::local_options(reproducible.useCache = FALSE)
  mine <- file.path(withr::local_tempdir(), "fireSense_spreadFit_x.log")
  suppressMessages(suppressWarnings(runDEoptim(
    landscape = NULL, annualDTx1000 = NULL, nonAnnualDTx1000 = NULL,
    fireBufferedListDT = NULL, historicalFires = NULL,
    itermax = 5, trace = FALSE, strategy = 2L, cores = c("hostA", "hostB"),
    paths = list(cachePath = withr::local_tempdir()),
    lower = lower, upper = lower + 1, mutuallyExclusive = NULL, formulaToFit = form,
    objFunCoresInternal = 1L, covMinMax = NULL, maxFireSpread = 0.3, Nreps = 1L, thresh = 512,
    logPath = mine, .verbose = FALSE)))
  expect_identical(got, mine)
})
