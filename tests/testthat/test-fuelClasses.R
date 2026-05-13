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
