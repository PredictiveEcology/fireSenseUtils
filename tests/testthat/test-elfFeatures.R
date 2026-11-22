## elfSpreadCoefs(), elfFeatures() and elfOutlines(): the per-ELF feature table and outlines for an interactive map.
## The ledger rows are small synthetic ones, as latestSpreadFits() returns them; Drive and the fire records are not read.

## A spread ledger row: `params` is a data.table with one row per parameter set (here two, so the mean is visible)
spreadRow <- function(id, climate = "CMD", link = "powerLawSmooth", scale = "linear", nf = c(nfLCC_100_40_60 = "100,40,60"),
                      extra = list(), spreadClimate = "ceiling", writtenAt = "2026-11-19 12:00:00") {
  one <- c(list(maxAsymptote = 0.25, hillSlope1 = 1, inflectionPoint1 = 1, `(Intercept)` = 4, youngAge = -2,
                dom_agb_Pice_mar = 3, sec_agb_Pice_gla = 50, yearSpreadSD = 1.5), stats::setNames(list(30), climate),
           stats::setNames(as.list(c(-1, 40)[seq_along(nf)]), names(nf)), extra)
  params <- data.table::as.data.table(rbind(one, lapply(one, `+`, 1)))
  params <- data.table::as.data.table(lapply(params, unlist))
  row <- data.frame(polygonID = id, spreadLink = link, spreadClimate = spreadClimate, fuelScale = scale,
                    writtenAt = as.POSIXct(writtenAt, tz = "UTC"))
  row$params <- list(params)
  row$nonForestedLCCGroups <- list(lapply(strsplit(sub("^nfLCC_", "", names(nf)), "_"), as.numeric) |>
                                     stats::setNames(names(nf)))
  row
}
bindRows <- function(...) do.call(rbind, list(...))

test_that("elfSpreadCoefs() gives the mean of each coefficient over a row's parameter sets, one row per ELF and coefficient", {
  rows <- bindRows(spreadRow("4.1"), spreadRow("4.3"))
  out <- elfSpreadCoefs(rows)
  expect_s3_class(out, "data.frame")
  expect_true(all(c("ELF", "coefficient", "value", "excluded") %in% names(out)))
  v <- function(elf, co) out$value[out$ELF == elf & out$coefficient == co]
  expect_equal(v("4.1", "maxAsymptote"), 0.75)       # mean of 0.25 and 1.25: the second set is the first plus 1
  expect_equal(v("4.1", "(Intercept)"), 4.5)
  expect_equal(v("4.3", "yearSpreadSD"), 2)
  expect_equal(v("4.1", "CMD"), 30.5)
  expect_true(all(is.na(out$excluded)))
})

test_that("elfSpreadCoefs() adds a row per land-cover class of a non-forest group, each with its group's value", {
  rows <- bindRows(spreadRow("4.1", nf = c(nfLCC_100_40_60 = "x", nfLCC_50 = "y")))
  out <- elfSpreadCoefs(rows)
  v <- function(co) out$value[out$coefficient == co]
  expect_equal(v("nfLCC_100_40_60"), -0.5)             # the group's own row stays
  for (cls in c("100", "40", "60")) expect_equal(v(paste0("nfLCC_", cls)), -0.5)
  expect_equal(v("nfLCC_50"), 40.5)                    # one class: the group is the class; not listed twice
  expect_identical(sum(out$coefficient == "nfLCC_50"), 1L)
})

test_that("elfSpreadCoefs() gives NA, excluded as 'different model', to ELFs whose model differs from the most common one", {
  rows <- bindRows(spreadRow("a"), spreadRow("b"), spreadRow("c"),
                   spreadRow("d", climate = "cumMDC"),             # another climate term
                   spreadRow("e", link = "logistic"),              # another link
                   spreadRow("f", scale = "log"),                  # another fuel scale
                   spreadRow("g", extra = list(newTerm = 1)),      # another non-fuel coefficient
                   spreadRow("h", spreadClimate = "activation"))   # another climate mode
  out <- elfSpreadCoefs(rows)
  expect_setequal(unique(out$ELF), letters[1:8])
  bad <- out[out$ELF %in% c("d", "e", "f", "g", "h"), ]
  expect_true(all(bad$excluded == "different model"))
  expect_true(all(is.na(bad$value)))
  expect_true(all(is.na(out$excluded[out$ELF %in% c("a", "b", "c")])))
  expect_false(anyNA(out$value[out$ELF %in% c("a", "b", "c")]))
  ## a different nonfuel set only: fuel terms may differ freely between ELFs
  rows <- bindRows(spreadRow("a"), spreadRow("b", extra = list(dom_agb_Pinu_con = 2, treedWetland_agb = 1)))
  expect_true(all(is.na(elfSpreadCoefs(rows)$excluded)))
})

test_that("elfSpreadCoefs() copes with rows that record no spreadClimate, and with no rows", {
  rows <- bindRows(spreadRow("a"), spreadRow("b"))
  rows$spreadClimate <- NULL
  expect_false(anyNA(elfSpreadCoefs(rows)$value))
  expect_identical(NROW(elfSpreadCoefs(NULL)), 0L)
})

test_that("elfCoefLabel() names a coefficient in words, with LandR species names", {
  expect_identical(elfCoefLabel(c("dom_agb_Pice_mar", "sec_agb_Pinu_con", "nfLCC_40", "youngAge", "CMD", "yearSpreadSD",
                                  "(Intercept)", "treedWetland_agb")),
                   c("Dominant AGB: Black spruce (Pice_mar)", "Secondary AGB: Lodgepole pine (Pinu_con)",
                     "Non-forest LCC 40", "Young stands", "CMD", "Year SD", "Intercept", "Treed wetland AGB"))
  ## a species group of codes LandR does not know keeps its codes
  expect_identical(elfCoefLabel("dom_agb_Pc_gl.Ab_la"), "Dominant AGB: Pc_gl + Ab_la")
})

fitDatesTable <- function(ELFs) {
  data.frame(ELF = ELFs, spreadFitDate = as.POSIXct(c("2026-11-19 10:00:00", NA, "2026-11-21 08:00:00")[seq_along(ELFs)], tz = "UTC"),
             spreadLedger = "s.rds", ignitionFitDate = as.POSIXct(c("2026-10-05 12:00:00", "2026-10-06 12:00:00", NA)[seq_along(ELFs)], tz = "UTC"),
             ignitionLedger = "i.rds")
}

test_that("elfFeatures() is one long table: Observed rates, Fit dates and Spread coefficients, with no row for a missing value", {
  local_mocked_bindings(
    ELFfireRates = function(ELFs, plan, fireYears, pixelAreaHa, destinationPath, ...)
      data.table::data.table(ELF = c("3.1.1_2", "5.1"), naturalFiresPerYear = c(5, 1), naturalFiresPer100km2PerYear = c(0.5, 0.25),
                             haBurnedPerYear = c(1000, 50), proportionBurnedPerYear = c(0.01, 0.002)),
    elfFitDates = function(spreadFolder, ignitionFolder, destinationPath, ELFs = NULL) fitDatesTable(c("3.1.1_2", "5.1", "9.9")),
    latestSpreadFits = function(...) bindRows(spreadRow("3.1.1_2"), spreadRow("5.1"), spreadRow("9.9")))
  plan <- data.table::data.table(action = "merge", ELF = "3.1.1_2", members = list(c("3.1.1", "3.1.2")), reason = "too few fires")
  out <- elfFeatures(toyELFs(), plan, fireYears = 2001:2003, pixelAreaHa = 2500, destinationPath = tempdir())
  expect_identical(names(out), c("ELF", "group", "feature", "label", "units", "type", "value"))
  expect_setequal(unique(out$group), c("Observed", "Fit", "Spread coefficients"))
  expect_false(anyNA(out$value))
  expect_false("9.9" %in% out$ELF)                                  # not on the map (after the merges)
  v <- function(elf, f) out$value[out$ELF == elf & out$feature == f]
  expect_equal(v("3.1.1_2", "percentBurnedPerYear"), 1)             # a percent, not a proportion
  expect_equal(v("5.1", "naturalFiresPer100km2PerYear"), 0.25)
  ## dates: days since 1970-01-01, type "date"
  expect_equal(v("3.1.1_2", "spreadFitDate"), as.numeric(as.Date("2026-11-19")))
  expect_identical(unique(out$type[out$feature == "spreadFitDate"]), "date")
  expect_identical(unique(out$type[out$feature == "naturalFiresPerYear"]), "numeric")
  expect_length(v("5.1", "spreadFitDate"), 0L)                      # no spread date: no row
  expect_identical(unique(out$label[out$feature == "ignitionFitDate"]), "Ignition and escape fit date")
  ## coefficients: a per-class row for each non-forest class, none for the group
  expect_true("nfLCC_40" %in% sub("^spread_", "", out$feature))
  expect_false(any(grepl("nfLCC_100_40_60", out$feature)))
  expect_equal(v("5.1", "spread_maxAsymptote"), 0.75)
  expect_identical(unique(out$label[out$feature == "spread_youngAge"]), "Young stands")
})

test_that("elfOutlines() gives one dissolved core outline per ELF, in longitude-latitude, with its name", {
  out <- elfOutlines(toyELFs(), tolerance = 100)
  expect_s4_class(out, "SpatVector")
  expect_identical(sort(out$ELF), sort(names(toyELFs()$rasWhole)))
  expect_equal(nrow(out), 4)
  expect_true(terra::is.lonlat(out))
  expect_true("name" %in% names(out))
  ## the core of 3.1.1 is 3 of 12 columns of the 5 km grid: 15 km x 20 km, not its buffer's 20 km x 20 km
  one <- out[out$ELF == "3.1.1", ]
  expect_equal(terra::expanse(one, unit = "km"), 15 * 20, tolerance = 0.05)
  ## a ring-shaped core is one polygon per ELF, not one per piece
  expect_equal(length(unique(out$ELF)), nrow(out))
})
