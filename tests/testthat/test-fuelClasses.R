test_that("combine_fuel_classes behaves reasonably", {
  withr::local_package("terra")

  ## testing this logic is easier than testing assessFuelClasses due to GLM
  tempDF <- data.table(
    species = c(
      "Pice_mar", "Pinu_con", "Popu_tre",
      "Betu_pap", "Pice_eng", "Pice_gla", "Abie_las"
    ),
    coef = c(-0.0132, 0.00154, -0.021, -0.0017, -0.0191, -0.0092, -0.0063),
    sign = c(
      "negative", "positive", "negative", "negative",
      "negative", "negative", "negative"
    ),
    FuelClass = c(
      "BlkSprc", "LdJkPine", "PopBrch", "PopBrch",
      "SprcFrLrch", "SprcFrLrch", "SprcFrLrch"
    ),
    above10PctRelB = c(0.066, 0.625, 0.223, 0.051, 0.095, 0.352, 0.560)
  )

  ## birch and aspen will be combined, fir will not be combined with Picea
  out <- combine_fuel_classes(tempDF)
  expect_true("Bt_pa.Pp_tr" %in% out$assignedFuelClass)
  expect_true("Abie_las" %in% out$assignedFuelClass)

  ## with birch positive instead of negative
  tempDF2 <- data.table(
    species = c(
      "Pice_mar", "Pinu_con", "Popu_tre",
      "Betu_pap", "Pice_eng", "Pice_gla", "Abie_las"
    ),
    coef = c(-0.0132, 0.00154, 0.021, -0.0017, -0.0191, -0.0092, -0.0063),
    sign = c(
      "negative", "positive", "negative", "positive",
      "negative", "negative", "negative"
    ),
    FuelClass = c(
      "BlkSprc", "LdJkPine", "PopBrch", "PopBrch",
      "SprcFrLrch", "SprcFrLrch", "SprcFrLrch"
    ),
    above10PctRelB = c(0.066, 0.625, 0.223, 0.051, 0.095, 0.352, 0.560)
  )

  ## birch should not be combined as sign is now different, fir combined instead
  out <- combine_fuel_classes(tempDF2)
  expect_true("Betu_pap" %in% out$assignedFuelClass)
  expect_false("Abie_las" %in% out$assignedFuelClass)

  tempDF3 <- data.table(
    species = c(
      "Pice_mar", "Pinu_con", "Popu_tre", "Pinu_ban",
      "Betu_pap", "Pice_eng", "Pice_gla", "Abie_las"
    ),
    coef = c(-0.0132, 0.00154, 0.021, 0.013, -0.0017, -0.0191, -0.0092, -0.0063),
    sign = c(
      "negative", "positive", "negative", "positive",
      "negative", "negative", "negative", "negative"
    ),
    FuelClass = c(
      "BlkSprc", "LdJkPine", "PopBrch", "LdJkPine",
      "PopBrch", "SprcFrLrch", "SprcFrLrch", "SprcFrLrch"
    ),
    above10PctRelB = c(0.066, 0.625, 0.223, 0.10, 0.051, 0.095, 0.352, 0.560)
  )

  ## jack and lodgepole should be combined along with engelmann and white spruce, birch and aspen
  ## Abie_las should NOT be combined
  out <- combine_fuel_classes(tempDF3)
  expect_true("Abie_las" %in% out$assignedFuelClass)
  expect_true("Pn_ba.Pn_co" %in% out$assignedFuelClass)
  expect_true("Bt_pa.Pp_tr" %in% out$assignedFuelClass)

  ## test that rare species are grouped first, regardless of sign
  tempDF4 <- data.table(
    species = c(
      "Pice_mar", "Pinu_con", "Popu_tre",
      "Betu_pap", "Pice_eng", "Pice_gla", "Abie_las"
    ),
    coef = c(-0.0132, 0.00154, 0.021, -0.0017, -0.0191, -0.0092, -0.0063),
    sign = c(
      "negative", "positive", "negative", "negative",
      "negative", "negative", "negative"
    ),
    FuelClass = c(
      "BlkSprc", "LdJkPine", "PopBrch", "PopBrch",
      "SprcFrLrch", "SprcFrLrch", "SprcFrLrch"
    ),
    above10PctRelB = c(0.066, 0.625, 0.223, 0.041, 0.095, 0.352, 0.560)
  )
  out <- combine_fuel_classes(tempDF4)

  expect_true("Bt_pa.Pp_tr" %in% out$assignedFuelClass)
  expect_false("Popu_tre" %in% out$assignedFuelClass)
  expect_true("Pc_en.Pc_gl" %in% out$assignedFuelClass)

  tempDF5 <- data.table(
    species = c(
      "Popu_bal", "Pinu_con", "Popu_tre",
      "Betu_pap", "Pice_eng", "Pice_gla", "Abie_las"
    ),
    coef = c(0.0132, 0.00154, 0.021, 0.0017, -0.0191, -0.0092, -0.0063),
    sign = c(
      "positive", "positive", "positive", "positive",
      "negative", "negative", "negative"
    ),
    FuelClass = c(
      "PopBrch", "LdJkPine", "PopBrch", "PopBrch",
      "SprcFrLrch", "SprcFrLrch", "SprcFrLrch"
    ),
    above10PctRelB = c(0.066, 0.625, 0.223, 0.041, 0.095, 0.352, 0.560)
  )
  out <- combine_fuel_classes(tempDF5)

  ## aspen is not merged because it has more biomass than birch and balsam poplar
  expect_true("Popu_tre" %in% out$assignedFuelClass)

  ## spruce is merged because engelmann is the second least abundant
  expect_true("Pc_en.Pc_gl" %in% out$assignedFuelClass)

  ## add yet another populus family - make sure sign is respected
  tempDF6 <- data.table(
    species = c(
      "Popu_bal", "Betu_all", "Pinu_con", "Popu_tre",
      "Betu_pap", "Pice_eng", "Pice_gla", "Abie_las"
    ),
    coef = c(0.0132, -0.023, 0.00154, 0.021, 0.0017, -0.0191, -0.0092, -0.0063),
    sign = c(
      "positive", "negative", "positive", "positive",
      "positive", "negative", "negative", "negative"
    ),
    FuelClass = c(
      "PopBrch", "PopBrch", "LdJkPine", "PopBrch",
      "PopBrch", "SprcFrLrch", "SprcFrLrch", "SprcFrLrch"
    ),
    above10PctRelB = c(0.066, 0.051, 0.625, 0.223, 0.041, 0.095, 0.352, 0.560)
  )
  out <- combine_fuel_classes(tempDF6)

  ## aspen is merged because it has the correct sign unlike birch
  expect_false("Popu_tre" %in% out$assignedFuelClass)

  ## spruce is merged because engelmann is the third least abundant
  expect_true("Pc_en.Pc_gl" %in% out$assignedFuelClass)

  ## one final test with different b threshold
  tempDF7 <- data.table(
    species = c(
      "Popu_bal", "Betu_all", "Pinu_con", "Popu_tre",
      "Betu_pap", "Pice_eng", "Pice_gla", "Abie_las"
    ),
    coef = c(0.0132, -0.023, 0.00154, 0.021, 0.0017, -0.0191, -0.0092, -0.0063),
    sign = c(
      "positive", "negative", "positive", "positive",
      "positive", "negative", "negative", "negative"
    ),
    FuelClass = c(
      "PopBrch", "PopBrch", "LdJkPine", "PopBrch",
      "PopBrch", "SprcFrLrch", "SprcFrLrch", "SprcFrLrch"
    ),
    above10PctRelB = c(0.066, 0.051, 0.625, 0.223, 0.041, 0.095, 0.352, 0.560)
  )
  out <- combine_fuel_classes(tempDF7, lowThreshold = 0.1)

  ## swamp birch is merged regardless of sign because it is under the b threshold
  expect_false("Betu_all" %in% out$assignedFuelClass)

  ## there are fewer than 5 because so many classes are rare
  expect_true(length(unique(out$assignedFuelClass)) == 4)

  ## cleanup
  withr::deferred_run()
})

test_that("cohortsToFuelClasses names a species mapped to two fuel classes", {
  withr::local_package("data.table")

  ## `unique(sppEquiv[, .(FuelClass, <sppEquivCol>)])` is over both columns, but the join
  ## keys on the species column alone -- so a species with two fuel classes survives as two
  ## rows and every one of its cohorts is multiplied. data.table then stops with
  ## "Join results in N rows; more than nrow(x)+nrow(i)", which names neither the species
  ## nor this function, and a 15-worker run lost five study areas to it before anyone could
  ## tell what it meant.
  ##
  ## Real case: LandR::sppEquivalencies_CA maps Pseu_men (Douglas-fir) to both "DgFrPoPine"
  ## and "CedrMplOther", so every area containing Douglas-fir failed and every area without
  ## it was fine.
  sppEquiv <- data.table(
    LandR     = c("Pice_mar", "Pinu_con", "Pseu_men",   "Pseu_men"),
    FuelClass = c("BlkSprc",  "LdJkPine", "DgFrPoPine", "CedrMplOther")
  )
  cohortData <- data.table(
    pixelGroup  = rep(1:2, each = 2),
    speciesCode = c("Pice_mar", "Pseu_men", "Pinu_con", "Pseu_men"),
    age         = c(50L, 60L, 70L, 80L),
    B           = c(100L, 200L, 300L, 400L)
  )

  err <- tryCatch(
    cohortsToFuelClasses(cohortData = cohortData, pixelGroupMap = NULL, flammableRTM = NULL,
                         sppEquiv = sppEquiv, sppEquivCol = "LandR", cutoffForYoungAge = 15L,
                         requiredFuelClasses = unique(sppEquiv$FuelClass)),
    error = function(e) conditionMessage(e))

  ## it stops, and the message is actionable: the species, both of its classes, and what to do
  expect_type(err, "character")
  expect_match(err, "Pseu_men", fixed = TRUE)
  expect_match(err, "DgFrPoPine", fixed = TRUE)
  expect_match(err, "CedrMplOther", fixed = TRUE)
  expect_match(err, "single FuelClass", fixed = TRUE)
  ## and it is NOT the opaque data.table message the user used to get
  expect_false(grepl("nrow(x)+nrow(i)", err, fixed = TRUE))

  ## a clean table gets past this guard (it fails later on the NULL rasters, not here)
  clean <- sppEquiv[!(LandR == "Pseu_men" & FuelClass == "CedrMplOther")]
  err2 <- tryCatch(
    cohortsToFuelClasses(cohortData = cohortData, pixelGroupMap = NULL, flammableRTM = NULL,
                         sppEquiv = clean, sppEquivCol = "LandR", cutoffForYoungAge = 15L,
                         requiredFuelClasses = unique(clean$FuelClass)),
    error = function(e) conditionMessage(e))
  expect_false(grepl("more than one FuelClass", paste(err2, collapse = " "), fixed = TRUE))
})

test_that("cohortsToFuelClasses survives 0-row cohortData", {
  withr::local_package("data.table")
  withr::local_package("terra")

  ## ELFs 3.2.1, 3.2.4, 3.2.5 and 3.3.2 have no tree species at all, so `cohortData` arrives
  ## with zero rows. FireSense can still fit on nonForest fuel only, so this must not stop.
  pixelGroupMap <- rast(nrows = 4, ncols = 4, vals = rep(1:2, 8))
  flammableRTM <- rast(pixelGroupMap, vals = 1)
  sppEquiv <- data.table(
    LandR     = c("Pice_mar", "Pinu_con"),
    FuelClass = c("BlkSprc",  "LdJkPine")
  )
  required <- c("BlkSprc", "LdJkPine")

  noCohorts <- data.table(pixelGroup = integer(), speciesCode = character(),
                          age = integer(), B = integer())
  out <- cohortsToFuelClasses(cohortData = noCohorts, pixelGroupMap = pixelGroupMap,
                              flammableRTM = flammableRTM, sppEquiv = sppEquiv,
                              sppEquivCol = "LandR", cutoffForYoungAge = 15L,
                              requiredFuelClasses = required)

  ## one zero-filled layer per required class (plus the youngAge layer the function always adds)
  expect_s4_class(out, "SpatRaster")
  expect_true(all(required %in% names(out)))
  expect_equal(sort(names(out)), sort(c(required, youngAgeTxt)))
  expect_true(all(values(out[[required]]) == 0))
  ## geometry comes from pixelGroupMap
  expect_true(compareGeom(out, pixelGroupMap, stopOnError = FALSE))

  ## the with-cohorts path is unchanged: biomass lands in the right class and pixel
  cohortData <- data.table(pixelGroup = c(1L, 2L), speciesCode = c("Pice_mar", "Pinu_con"),
                           age = c(50L, 60L), B = c(100L, 200L))
  out2 <- cohortsToFuelClasses(cohortData = cohortData, pixelGroupMap = pixelGroupMap,
                               flammableRTM = flammableRTM, sppEquiv = sppEquiv,
                               sppEquivCol = "LandR", cutoffForYoungAge = 15L,
                               requiredFuelClasses = required)
  expect_equal(sort(names(out2)), sort(c(required, youngAgeTxt)))
  expect_equal(unname(values(out2[[required]])[1, ]), c(100, 0))
  expect_equal(unname(values(out2[[required]])[2, ]), c(0, 200))
})

test_that("cohortsToFuelClasses with no tree species and no required classes gives youngAge only", {
  withr::local_package("data.table")
  withr::local_package("terra")

  ## The fitting path (fireSense_dataPrepFit -> fireSenseCovariatesCreate) passes NO
  ## requiredFuelClasses. With an ELF that has no tree species, sppEquiv has zero rows too, so
  ## there is not a single tree fuel-class layer to stack. The result must still carry the
  ## youngAge layer, on pixelGroupMap's geometry and NA mask, so the nonForest covariates
  ## can be built on top of it.
  pixelGroupMap <- rast(nrows = 4, ncols = 4, vals = 0L)
  pixelGroupMap[1:3] <- NA
  flammableRTM <- rast(pixelGroupMap, vals = 1)
  flammableRTM[1:3] <- NA
  sppEquiv <- data.table(LandR = character(), FuelClass = character())
  noCohorts <- data.table(pixelGroup = integer(), speciesCode = character(),
                          age = integer(), B = integer())

  out <- suppressWarnings( ## max(age) over zero rows
    cohortsToFuelClasses(cohortData = noCohorts, pixelGroupMap = pixelGroupMap,
                         flammableRTM = flammableRTM, sppEquiv = sppEquiv,
                         sppEquivCol = "LandR", cutoffForYoungAge = 15L,
                         requiredFuelClasses = NULL)
  )

  expect_s4_class(out, "SpatRaster")
  expect_identical(names(out), youngAgeTxt)
  expect_true(compareGeom(out, pixelGroupMap, stopOnError = FALSE))
  expect_identical(is.na(values(out, mat = FALSE)), is.na(values(pixelGroupMap, mat = FALSE)))
  expect_true(all(values(out, mat = FALSE) == 0, na.rm = TRUE))
})

test_that("assessFuelClasses with no tree species returns the non-forest groups only", {
  withr::local_package("data.table")
  set.seed(1)
  ## a landscape of non-forest pixels only (B is NA everywhere), three land covers that burn
  ## at different rates so the k-means on the glm coefficients has something to cluster
  n <- 300L
  landscape <- data.table(
    cell = seq_len(n), speciesCode = NA_character_,
    lcc = rep(c(40L, 50L, 80L), each = n / 3L),
    B = NA_integer_, totalBiomass = NA_integer_, year = 2020L
  )
  landscape[, burned := rbinom(.N, 1, c(`40` = 0.05, `50` = 0.3, `80` = 0.6)[as.character(lcc)])]
  noSpp <- data.table(LandR = character(0), FuelClass = character(0))

  out <- assessFuelClasses(landscape = landscape, fuelCol = "FuelClass", sppEquiv = noSpp,
                           sppEquivCol = "LandR", nonforestLCC = c(40L, 50L, 80L))

  expect_named(out, c("modSppEquiv", "nonForestedLCCGroups", "missingLCCgroup"))
  expect_identical(nrow(out$modSppEquiv), 0L)
  expect_true(all(c("species", "assignedFuelClass", "FuelClass") %in% names(out$modSppEquiv)))
  expect_length(out$nonForestedLCCGroups, 2L)
  expect_setequal(unlist(out$nonForestedLCCGroups), c(40, 50, 80))
  expect_true(out$missingLCCgroup %in% names(out$nonForestedLCCGroups))

  ## with species declared but no forested pixel left, the land-cover-code mismatch stop stays
  withSpp <- data.table(LandR = "Pice_mar", FuelClass = "BlkSprc")
  expect_error(
    assessFuelClasses(landscape = landscape, fuelCol = "FuelClass", sppEquiv = withSpp,
                      sppEquivCol = "LandR", nonforestLCC = c(40L, 50L, 80L)),
    "no forested pixels with a species remain"
  )
})

test_that("lccFlammableShare is each class's share of the flammable pixels", {
  skip_if_not_installed("terra")
  r <- terra::rast(nrows = 2, ncols = 5, vals = c(40, 40, 50, 81, 81, 81, 0, NA, 50, 50))
  sh <- lccFlammableShare(r)
  expect_equal(sh, c(`40` = 2, `50` = 3, `81` = 3) / 8)
  ## a flammable map restricts it further
  fl <- terra::rast(r, vals = c(1, 1, 1, 1, 1, 1, 1, 1, 0, 0))
  expect_equal(lccFlammableShare(r, fl), c(`40` = 2, `50` = 1, `81` = 3) / 6)
})

nfLandscape <- function(rates, n = 600L) {
  set.seed(2)
  landscape <- data.table::data.table(cell = seq_len(n), speciesCode = NA_character_,
                                      lcc = rep(as.integer(names(rates)), each = n / length(rates)),
                                      B = NA_integer_, totalBiomass = NA_integer_, year = 2020L)
  landscape[, burned := rbinom(.N, 1, rates[as.character(lcc)])]
  landscape[]
}
noSpp <- function() data.table::data.table(LandR = character(0), FuelClass = character(0))

test_that("a non-forest class under minCovariateProp of the ELF joins the group nearest its burn rate", {
  withr::local_package("data.table")
  ## 40 burns like 80; it is 1% of the ELF's flammable pixels, so it is not clustered on its own
  rates <- c(`40` = 0.6, `50` = 0.05, `80` = 0.6, `100` = 0.3)
  share <- c(`40` = 0.01, `50` = 0.3, `80` = 0.4, `100` = 0.29)
  out <- assessFuelClasses(landscape = nfLandscape(rates), fuelCol = "FuelClass", sppEquiv = noSpp(),
                           sppEquivCol = "LandR", nonforestLCC = c(40L, 50L, 80L, 100L),
                           lccShare = share, minCovariateProp = 0.05)
  grp <- out$nonForestedLCCGroups
  expect_setequal(unlist(grp), c(40, 50, 80, 100))
  in40 <- names(grp)[vapply(grp, function(g) 40 %in% g, logical(1))]
  expect_length(in40, 1L)
  expect_true(80 %in% grp[[in40]])
  expect_false(identical(grp[[in40]], 40))
})

test_that("when every non-forest class is rare, they form one group", {
  withr::local_package("data.table")
  rates <- c(`40` = 0.6, `50` = 0.05)
  out <- assessFuelClasses(landscape = nfLandscape(rates), fuelCol = "FuelClass", sppEquiv = noSpp(),
                           sppEquivCol = "LandR", nonforestLCC = c(40L, 50L),
                           lccShare = c(`40` = 0.02, `50` = 0.01), minCovariateProp = 0.05)
  expect_length(out$nonForestedLCCGroups, 1L)
  expect_setequal(unlist(out$nonForestedLCCGroups), c(40, 50))
  expect_identical(out$missingLCCgroup, names(out$nonForestedLCCGroups))
})

test_that("missingForest joins the group with its nearest class even when one code is a prefix of another", {
  withr::local_package("data.table")
  ## codes 10 and 100: the pattern "10" is a substring of "100". 10 burns alone; 50 and 100 burn alike.
  ## The forest-less pixels (lcc 99, not in nonforestLCC) burn like 10.
  n <- 800L
  landscape <- data.table(cell = seq_len(n), speciesCode = NA_character_,
                          lcc = rep(c(10L, 50L, 100L, 99L), each = n / 4L),
                          B = NA_integer_, totalBiomass = NA_integer_, year = 2020L)
  set.seed(3)
  landscape[, burned := rbinom(.N, 1, c(`10` = 0.02, `50` = 0.5, `100` = 0.6, `99` = 0.02)[as.character(lcc)])]
  out <- assessFuelClasses(landscape = landscape, fuelCol = "FuelClass", sppEquiv = noSpp(),
                           sppEquivCol = "LandR", nonforestLCC = c(10L, 50L, 100L))
  grp <- out$nonForestedLCCGroups
  expect_length(out$missingLCCgroup, 1L)
  expect_true(10 %in% grp[[out$missingLCCgroup]])
  expect_false(100 %in% grp[[out$missingLCCgroup]])
})

test_that("assessFuelClasses gives the same groups for the same data whatever the global seed, and leaves the RNG alone", {
  withr::local_package("data.table")
  ## 12 classes with irregular burn rates: k-means from random starts gives several partitions
  burnRates <- c(0.07, 0.07, 0.21, 0.28, 0.30, 0.39, 0.56, 0.64, 0.71, 0.71, 0.83, 0.88)
  codes <- seq(20L, by = 10L, length.out = length(burnRates))
  rates <- stats::setNames(burnRates, codes)
  landscape <- nfLandscape(rates, n = 600L * length(codes)) # sets its own seed: build it before varying the global one
  run <- function(seed) {
    set.seed(seed)
    assessFuelClasses(landscape = landscape, fuelCol = "FuelClass", sppEquiv = noSpp(),
                      sppEquivCol = "LandR", nonforestLCC = codes,
                      targetNonForestClasses = 3)$nonForestedLCCGroups
  }
  ref <- run(1)
  for (s in 2:25) expect_identical(run(s), ref)

  set.seed(99)
  before <- .Random.seed
  set.seed(99)
  assessFuelClasses(landscape = landscape, fuelCol = "FuelClass", sppEquiv = noSpp(),
                    sppEquivCol = "LandR", nonforestLCC = codes, targetNonForestClasses = 3)
  expect_identical(.Random.seed, before)
})
