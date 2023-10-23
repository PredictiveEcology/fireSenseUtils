require("terra")
test_that("makeLandcoverDT works", {
  
  flammableRTM <- terra::rast(xmin = 0, xmax = 15, ymin = 0, ymax = 15, res = c(1,1))
  flammableRTM[] <- sample(0:1, ncell(flammableRTM), replace = TRUE)
  rstLCC <- terra::rast(xmin = 0, xmax = 20, ymin = 0, ymax = 20, res = c(1,1))
  rstLCC[] <- sample(1:8, size = ncell(rstLCC), replace = TRUE)
  nonForestLCCgroups <- list("foo1" = 1:2, "foo2" = 3:4)
  forestedLCC <- c(5, 6)
  out <- makeLandcoverDT(rstLCC = rstLCC, flammableRTM = flammableRTM,
                         forestedLCC = forestedLCC, 
                         nonForestedLCCGroups = nonForestLCCgroups)  
  expect_true(sum(flammableRTM[]) == nrow(out))
      
})
