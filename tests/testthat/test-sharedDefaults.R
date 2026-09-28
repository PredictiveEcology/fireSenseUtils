## The shared defaults must keep the values the modules used before they took their defaults
## from these constants; changing one changes every fit and prediction.
test_that("shared module defaults keep their values", {
  expect_identical(fireSenseForestedLCC, c(81, 210, 220, 230, 240))
  expect_identical(fireSenseYoungAgeCutoff, 15)
  expect_identical(fireSenseNonForestCanBeYoungAge, TRUE)
  expect_identical(fireSenseFlammabilityThreshold, 0.1)
  expect_identical(fireSenseFuelClassCol, "FuelClass")
  expect_identical(fireSenseIgAggFactor, 4)
  expect_identical(fireSenseSCANFIVersion, "V3")
})

test_that("functions sharing an argument default to the shared constant", {
  expect_identical(formals(makeTSD)$cutoffForYoungAge, quote(fireSenseYoungAgeCutoff))
  expect_identical(formals(castCohortData)$cutoffForYoungAge, quote(fireSenseYoungAgeCutoff))
  expect_identical(formals(cohortsToFuelClasses)$fuelClassCol, quote(fireSenseFuelClassCol))
  expect_identical(formals(makeFireSenseLCC)$flammabilityThreshold, quote(fireSenseFlammabilityThreshold))
  expect_identical(formals(makeFireSenseLCC)$scanfiVersion, quote(fireSenseSCANFIVersion))
})

test_that("forested and non-flammable codes do not overlap", {
  expect_length(intersect(fireSenseForestedLCC, fireSenseNonflammableLCC), 0)
})
