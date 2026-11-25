# Extracted from test-runDEoptim-paramNames.R:29

# setup ------------------------------------------------------------------------
library(testthat)
test_env <- simulate_test_env(package = "fireSenseUtils", path = "..")
attach(test_env, warn.conflicts = FALSE)

# test -------------------------------------------------------------------------
form <- ~ 0 + CMD + youngAge + dom_agb_Pice_mar
covs <- c("CMD", "youngAge", "dom_agb_Pice_mar")
lower <- stats::setNames(c(0.25, 0.1, 0, 0, 0, 0),
                           c("maxAsymptote", "inflectionPoint1", covs, "yearSpreadSD"))
testthat::local_mocked_bindings(
    clusterSetup = function(...) list(itermax = 5, trace = FALSE, strategy = 2L, NP = 40L),
    DEoptimIterative = function(fn, lower, upper, control, ...) list(),
    .package = "clusters")
withr::local_options(reproducible.useCache = FALSE)
msgs <- testthat::capture_messages(suppressWarnings(runDEoptim(
    landscape = NULL, annualDTx1000 = NULL, nonAnnualDTx1000 = NULL,
    fireBufferedListDT = NULL, historicalFires = NULL,
    itermax = 5, trace = FALSE, strategy = 2L, cores = c("hostA", "hostB"),
    paths = list(cachePath = withr::local_tempdir()),
    lower = lower, upper = lower + 1, mutuallyExclusive = NULL, formulaToFit = form,
    objFunCoresInternal = 1L, covMinMax = NULL, maxFireSpread = 0.3, Nreps = 1L, thresh = 512,
    logPath = file.path(withr::local_tempdir(), "fit.log"), .verbose = FALSE)))
msgs <- paste(msgs, collapse = "")
expect_match(msgs, paste0("Fitting 6 parameters: logistic: maxAsymptote, inflectionPoint1; ",
                            "covariates: CMD, youngAge, dom_agb_Pice_mar; year effect: yearSpreadSD"),
               fixed = TRUE)
expect_match(msgs, "objectiveFunction threshold SNLL to run all years after first 2 years: 512",
               fixed = TRUE)
