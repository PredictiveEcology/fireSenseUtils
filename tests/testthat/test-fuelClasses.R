require("terra")
test_that("combine_fuel_classes behaves reasonably", {
  #testing this logic is easier than testing assessFuelClasses due to GLM
  tempDF <- data.table(species =c ("Pice_mar", "Pinu_con", "Popu_tre",
                                    "Betu_pap", "Pice_eng", "Pice_gla", "Abie_las"),
                        coef = c(-0.0132, 0.00154, -0.021, -0.0017, -0.0191, -0.0092, -0.0063),
                        sign = c("negative", "positive", "negative", "negative",
                                 "negative", "negative", "negative"),
                        FuelClass = c("BlkSprc", "LdJkPine", "PopBrch", "PopBrch",
                                      "SprcFrLrch", "SprcFrLrch", "SprcFrLrch"),
                        above10PctRelB = c(0.066, 0.625, 0.223, 0.051, 0.095, 0.352, 0.560))

   #birch and aspen will be combined, fir will not be combined with Picea
   out <- combine_fuel_classes(tempDF)
   expect_true("PopBrch" %in% out$assignedFuelClass)
   expect_true("Abie_las" %in% out$assignedFuelClass)

   #with birch positive instead of negative
   tempDF2 <- data.table(species =c ("Pice_mar", "Pinu_con", "Popu_tre",
                                     "Betu_pap", "Pice_eng", "Pice_gla", "Abie_las"),
                         coef = c(-0.0132, 0.00154, 0.021, -0.0017, -0.0191, -0.0092, -0.0063),
                         sign = c("negative", "positive", "negative", "positive",
                                  "negative", "negative", "negative"),
                         FuelClass = c("BlkSprc", "LdJkPine", "PopBrch", "PopBrch",
                                       "SprcFrLrch", "SprcFrLrch", "SprcFrLrch"),
                         above10PctRelB = c(0.066, 0.625, 0.223, 0.051, 0.095, 0.352, 0.560))
   #birch should not be combined as sign is now different, fir combined instead
   out <- combine_fuel_classes(tempDF2)
   expect_true("Betu_pap" %in% out$assignedFuelClass)
   expect_false("Abie_las" %in% out$assignedFuelClass)

   tempDF3 <- data.table(species =c ("Pice_mar", "Pinu_con", "Popu_tre", "Pinu_ban",
                                     "Betu_pap", "Pice_eng", "Pice_gla", "Abie_las"),
                         coef = c(-0.0132, 0.00154, 0.021, 0.013, -0.0017, -0.0191, -0.0092, -0.0063),
                         sign = c("negative", "positive", "negative", "positive",
                                  "negative", "negative", "negative", "negative"),
                         FuelClass = c("BlkSprc", "LdJkPine", "PopBrch", "LdJkPine",
                                       "PopBrch", "SprcFrLrch", "SprcFrLrch", "SprcFrLrch"),
                         above10PctRelB = c(0.066, 0.625, 0.223, 0.10, 0.051, 0.095, 0.352, 0.560))
  #jack and lodgepole should be combined along with engelmann and white spruce, birch and aspen
  #Abie_las should NOT be combined
  out <- combine_fuel_classes(tempDF3)
  expect_true("Abie_las" %in% out$assignedFuelClass)
  expect_true("LdJkPine" %in% out$assignedFuelClass)
  expect_true("PopBrch" %in% out$assignedFuelClass)

})
