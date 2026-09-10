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
