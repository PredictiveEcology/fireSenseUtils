## elfSpreadCoefs(), elfFeatures() and elfOutlines(): the per-ELF feature table and outlines for an interactive map.
## The ledger rows are small synthetic ones, as latestSpreadFits() returns them; Drive and the fire records are not read.

## A spread ledger row: `params` is a data.table with one row per parameter set (here two, so the mean is visible)
spreadRow <- function(id, climate = "CMD", link = "powerLawSmooth", scale = "linear", nf = c(nfLCC_100_40_60 = "100,40,60"),
                      extra = list(), spreadClimate = "ceiling", writtenAt = "2026-11-19 12:00:00", sppEquiv = NULL) {
  one <- c(list(maxAsymptote = 0.25, hillSlope1 = 1, inflectionPoint1 = 1, `(Intercept)` = 4, youngAge = -2,
                dom_agb_Pice_mar = 3, sec_agb_Pice_gla = 50, yearSpreadSD = 1.5), stats::setNames(list(30), climate),
           stats::setNames(as.list(c(-1, 40)[seq_along(nf)]), names(nf)), extra)
  params <- data.table::as.data.table(rbind(one, lapply(one, `+`, 1)))
  params <- data.table::as.data.table(lapply(params, unlist))
  row <- data.frame(polygonID = id, spreadLink = link, spreadClimate = spreadClimate, fuelScale = scale,
                    writtenAt = as.POSIXct(writtenAt, tz = "UTC"))
  row$params <- list(params)
  if (!is.null(sppEquiv)) row$sppEquiv <- list(sppEquiv)
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

test_that("elfSpreadCoefs() gives NA, excluded as 'different model', to ELFs whose model form differs from the most common one", {
  rows <- bindRows(spreadRow("a"), spreadRow("b"), spreadRow("c"),
                   spreadRow("e", link = "logistic"),              # another link
                   spreadRow("f", scale = "log"),                  # another fuel scale
                   spreadRow("h", spreadClimate = "activation"),   # another climate mode
                   spreadRow("i", extra = list(CMD_yr = 1, CMD_dev = 2)))   # the climate split
  out <- elfSpreadCoefs(rows)
  expect_setequal(unique(out$ELF), c("a", "b", "c", "e", "f", "h", "i"))
  bad <- out[out$ELF %in% c("e", "f", "h", "i"), ]
  expect_true(all(bad$excluded == "different model"))
  expect_true(all(is.na(bad$value)))
  expect_true(all(is.na(out$excluded[out$ELF %in% c("a", "b", "c")])))
  expect_false(anyNA(out$value[out$ELF %in% c("a", "b", "c")]))
})

test_that("elfSpreadCoefs() keeps an ELF that fitted another climate variable or other fuel terms; each stays its own coefficient", {
  rows <- bindRows(spreadRow("a"), spreadRow("b"),
                   spreadRow("d", climate = "cumMDC"),                                       # another climate variable
                   spreadRow("g", extra = list(newTerm = 1, dom_agb_Pinu_con = 2, treedWetland_agb = 1)))
  out <- elfSpreadCoefs(rows)
  expect_true(all(is.na(out$excluded)))
  has <- function(elf, co) co %in% out$coefficient[out$ELF == elf]
  expect_true(has("d", "cumMDC") && !has("d", "CMD"))
  expect_true(has("a", "CMD") && !has("a", "cumMDC"))
  expect_true(has("g", "newTerm") && has("g", "dom_agb_Pinu_con"))
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
  ## a class of codes LandR does not know keeps its code
  expect_identical(elfCoefLabel("dom_agb_Pc_gl.Ab_la"), "Dominant AGB: Pc_gl.Ab_la")
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
  expect_identical(names(out), c("ELF", "group", "feature", "label", "units", "type", "value", "note"))
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

test_that("a grouped fuel term gives one coefficient per member species, noted with the species it is shared with", {
  ## the ledger row's sppEquiv maps each species to its fuel class; the term is named for the class
  sppEquiv <- data.frame(LandR = c("Popu_tre", "Betu_pap", "Popu_bal", "Pice_mar", "Pinu_con", "Pinu_con"),
                         FuelClass = c(rep("Pp_tr.Bt_pa.Pp_ba", 3), "Pice_mar", "Pinu_con", "Pinu_con"))
  rows <- bindRows(spreadRow("5.1", extra = list(`dom_agb_Pp_tr.Bt_pa.Pp_ba` = 7, sec_agb_Pinu_con = 5), sppEquiv = sppEquiv))
  out <- elfSpreadCoefs(rows)
  expect_true("note" %in% names(out))
  expect_false(any(grepl("\\.", out$coefficient[grepl("_agb_", out$coefficient)])))    # no group-level coefficient
  dom <- out[grepl("^dom_agb_", out$coefficient), ]
  expect_setequal(dom$coefficient, paste0("dom_agb_", c("Popu_tre", "Betu_pap", "Popu_bal", "Pice_mar")))
  expect_equal(dom$value[dom$coefficient == "dom_agb_Popu_tre"], 7.5)      # the term's mean
  expect_equal(dom$value[dom$coefficient == "dom_agb_Pice_mar"], 3.5)
  note <- function(co) out$note[out$coefficient == co]
  expect_identical(note("dom_agb_Popu_tre"), "shared with Betu_pap, Popu_bal")
  expect_identical(note("dom_agb_Betu_pap"), "shared with Popu_tre, Popu_bal")
  expect_identical(note("dom_agb_Pice_mar"), "")
  expect_identical(note("sec_agb_Pinu_con"), "")                             # a species in two rows of sppEquiv is still pure
  expect_identical(note("youngAge"), "")

  local_mocked_bindings(
    ELFfireRates = function(...) data.table::data.table(ELF = "5.1", naturalFiresPerYear = 1, naturalFiresPer100km2PerYear = 1,
                                                        haBurnedPerYear = 1, proportionBurnedPerYear = 0.01),
    latestSpreadFits = function(...) rows)
  feat <- elfFeatures(toyELFs(), NULL, fireYears = 2001:2003, pixelAreaHa = 2500, destinationPath = tempdir(),
                      dates = fitDatesTable(c("5.1", "9.9", "3.1.1")))
  expect_identical(names(feat), c("ELF", "group", "feature", "label", "units", "type", "value", "note"))
  one <- feat[feat$feature == "spread_dom_agb_Popu_tre", ]
  expect_identical(one$label, "Dominant AGB: Trembling aspen (Popu_tre)")
  expect_identical(one$note, "shared with Betu_pap, Popu_bal")
  expect_false(any(grepl("spread_dom_agb_Pp_tr", feat$feature, fixed = TRUE)))
  expect_identical(sum(grepl("^spread_dom_agb_", feat$feature)), 4L)
  expect_identical(unique(feat$note[feat$group == "Observed"]), "")
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

test_that("elfFeatures(dates =) uses the fit dates it is given, and does not read them again", {
  local_mocked_bindings(
    ELFfireRates = function(...) data.table::data.table(ELF = "5.1", naturalFiresPerYear = 1, naturalFiresPer100km2PerYear = 1,
                                                        haBurnedPerYear = 1, proportionBurnedPerYear = 0.01),
    elfFitDates = function(...) stop("the ledgers were read again"),
    latestSpreadFits = function(...) NULL)
  out <- elfFeatures(toyELFs(), NULL, fireYears = 2001:2003, pixelAreaHa = 2500, destinationPath = tempdir(),
                     dates = fitDatesTable(c("5.1", "9.9", "3.1.1")))
  v <- function(elf, f) out$value[out$ELF == elf & out$feature == f]
  expect_equal(v("5.1", "spreadFitDate"), as.numeric(as.Date("2026-11-19")))
  expect_equal(v("5.1", "ignitionFitDate"), as.numeric(as.Date("2026-10-05")))
  expect_false("9.9" %in% out$ELF)
})

test_that("pixelAreaHa() is the area of a cell in hectares, and noELFPlan() a plan with no rows", {
  expect_equal(pixelAreaHa(terra::rast(resolution = 240, crs = "EPSG:3978", extent = terra::ext(0, 2400, 0, 2400))), 5.76)
  expect_equal(pixelAreaHa(toyGrid()), 2500)                       # 5 km cells
  p <- noELFPlan()
  expect_identical(nrow(p), 0L)
  expect_true(all(c("action", "ELF", "members") %in% names(p)))
  expect_identical(names(mergeELFs(toyELFs(), p)$rasWhole), names(toyELFs()$rasWhole))   # merges nothing
})
