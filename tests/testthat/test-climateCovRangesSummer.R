## Summer temperature and precipitation can be the spread climate covariate (fireSense_dataPrepFit
## `spreadClimateCandidates`); a covariate with no entry in climateCovRanges stops the spread fit.
test_that("climateCovRanges has Tmaxsm and PPTsm, ", {
  expect_true(all(c("Tmaxsm", "PPTsm") %in% names(climateCovRanges)))
  expect_identical(climateCovRanges$Tmaxsm, c(0, 10))
  expect_identical(climateCovRanges$PPTsm, c(0, 100))
})

test_that("the names are those the spread tables use: ClimateNA Tmax_sm and PPT_sm without the underscore", {
  skip_if_not_installed("climateData")
  nms <- climateNAnames(c("Tmaxsm", "PPTsm"))
  expect_identical(nms, c("Tmax_sm", "PPT_sm"))
  expect_identical(gsub("_", "", nms), c("Tmaxsm", "PPTsm"))
})
