## A toy fit for the spread-fit validation tests: two fire years, a fuel and youngAge, one climate covariate.
toyValidationInputs <- function() {
  dt <- data.table::data.table
  landscape <- terra::rast(nrows = 40, ncols = 40, xmin = 0, xmax = 40, ymin = 0, ymax = 40, vals = 1)
  cellOf <- function(row, col) (row - 1L) * 40L + col
  block <- function(r0, c0, n) as.integer(outer(seq(r0, r0 + n - 1L), seq(c0, c0 + n - 1L), cellOf))
  ## each fire: a 15 x 15 buffer, of which the central 5 x 5 burned; the ignition at the centre
  fire <- function(r0, c0, id) {
    px <- block(r0, c0, 15L)
    dt(pixelID = px, buffer = as.integer(px %in% block(r0 + 5L, c0 + 5L, 5L)), ids = id)
  }
  fb <- list(year2001 = rbind(fire(1L, 1L, 1L), fire(20L, 20L, 2L)),
             year2002 = fire(10L, 5L, 3L))
  hist <- list(year2001 = data.frame(size = c(25L, 25L), cells = c(cellOf(8L, 8L), cellOf(27L, 27L)),
                                     ids = 1:2, date = "year2001"),
               year2002 = data.frame(size = 25L, cells = cellOf(17L, 12L), ids = 3L, date = "year2002"))
  set.seed(10)
  allPx <- sort(unique(unlist(lapply(fb, `[[`, "pixelID"))))
  nonAnnual <- list(`2000` = dt(   # named by data year; fire years 2001 and 2002 use the one strictly before
    pixelID = allPx,
    agb = as.integer(round(ifelse(stats::runif(length(allPx)) < 0.3, 0, stats::runif(length(allPx), 0, 5000)) * 1000)),
    youngAge = as.integer(stats::runif(length(allPx)) < 0.1) * 1000L))
  annual <- lapply(fb, function(b) dt(pixelID = b$pixelID, clim = as.integer(stats::runif(nrow(b), 0, 100) * 1000)))
  list(landscape = landscape, annualDTx1000 = annual, nonAnnualDTx1000 = nonAnnual,
       historicalFires = hist, fireBufferedListDT = fb,
       formulaToFit = "~ 0 + clim + youngAge + agb",
       covMinMax = data.table::data.table(clim = c(0, 100), youngAge = c(0, 1), agb = c(0, 1e4)),
       mutuallyExclusive = list(youngAge = "agb"),
       par = c(maxAsymptote = 0.26, clim = 2, youngAge = -2, agb = 3))
}
