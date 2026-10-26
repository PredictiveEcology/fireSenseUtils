# Extracted from test-bufferIndependentObjective.R:84

# setup ------------------------------------------------------------------------
library(testthat)
test_env <- simulate_test_env(package = "fireSenseUtils", path = "..")
attach(test_env, warn.conflicts = FALSE)

# prequel ----------------------------------------------------------------------
fixture <- readRDS(test_path("fixtures", "bufferIndependentObjective.rds"))
h <- 0.25
shift <- log(dnorm(0) / h)

# test -------------------------------------------------------------------------
sc <- fixture$scores
total <- 0
for (y in unique(sc$obs$yr)) {
    o <- sc$obs[sc$obs$yr == y, ]
    s <- sc$sims[sc$sims$yr == y, ]
    len <- vapply(o$ids, function(i) max(s$sim[s$ids == i]) + 1L, numeric(1))
    first <- cumsum(c(1, head(len, -1)))
    cellsOf <- function(k) first[k] + seq_len(len[k]) - 1
    buf <- data.table::rbindlist(lapply(seq_along(o$ids), function(k)
      data.table::data.table(ids = o$ids[k], pixelID = cellsOf(k), buffer = 0L))) # pixel ids as numeric
    buf[, `:=`(edge = pixelID %in% first, open = pixelID %in% first)]
    fires <- data.table::data.table(cells = first + 1, size = o$size, ids = o$ids)
    r <- terra::rast(nrows = ceiling(sum(len) / 1000), ncols = 1000)
    nCall <- 0L
    local_mocked_bindings(
      cropToCells = function(r, cells, ...) list(r = r, ncell = terra::ncell(r), toCrop = identity, toFull = identity),
      paramsSeparate = function(...) list(logisticPars = c(0.27, 1, 1), covPars = 1),
      spreadProbFromIntegerCovs = function(...) data.table::data.table(pixelID = buf$pixelID, cov = 0),
      logisticAll = function(...) seq(0.06, 0.08, length.out = nrow(buf)),
      .package = "fireSenseUtils")
    local_mocked_bindings(
      spreadCpp = function(landscape, loci, ...) {
        nCall <<- nCall + 1L
        data.table::rbindlist(lapply(seq_along(o$ids), function(k) {
          z <- s[s$ids == o$ids[k] & s$rep == nCall, ]
          cells <- if (z$cut) cellsOf(k)[seq_len(z$sim)] else cellsOf(k)[1L + seq_len(z$sim)]
          data.table::data.table(initialLocus = loci[k], indices = cells)
        }))
      },
      .package = "SpaDES.tools")
    out <- fireSenseUtils:::objFunInner(
      yr = sub("year", "", y), annDTx1000 = NULL, par = 1, parsModel = NULL,
      annualFires = fires, nonAnnualDTx1000 = NULL, annualFireBufferedDT = buf,
      indexNonAnnual = NULL, colsToUse = "cov", covMinMax = NULL, mutuallyExclusive = NULL, doAssertions = FALSE,
      maxFireSpread = 1, lowerSpreadProb = 0.05, cells = numeric(terra::ncell(r)), lanscape1stQuantileThresh = 1,
      weighted = FALSE, r = r, Nreps = sc$Nreps, doSNLL_FSTest = TRUE, doMADTest = FALSE, plot.it = FALSE,
      verbose = 0, sizeLik = "kernel")
    expect_equal(out$runaways, sum(s$cut))
    total <- total + out$SNLL_FS
  }
