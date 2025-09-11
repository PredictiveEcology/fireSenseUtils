test_that("makeFireSenseLCC works", {
  ## uses the NTEMS LCC which can be a lengthy download time
  testthat::skip_if_not(interactive())
  testthat::skip_on_cran()
  testthat::skip_on_ci()

  withr::local_package("terra")
  withr::local_package("LandR")

  dPath <- tempdir()

  SA  <- LandR::randomStudyArea(size = 1e8, seed = 1000)
  RTM <- terra::mask(terra::rast(SA, vals = 1, resolution = c(250, 250)), SA)
  flamLCC <- makeFireSenseLCC(neededYear = 2011, rasterToMatch = RTM,
                              studyArea = SA, destinationPath = dPath)

  flamLCC2 <- makeFireSenseLCC(neededYear = 2011, rasterToMatch = RTM,
                               studyArea = SA, destinationPath = dPath,
                               flammabilityThreshold = 0.4)

  #expect that the difference of the two is not 0
  flamDiff <- flamLCC - flamLCC2
  expect_true(sum(as.vector(flamLCC), na.rm = TRUE) > sum(as.vector(flamLCC2), na.rm = TRUE))
  # rm(flamLCC, flamLCC2, flamDiff)
  flamOtherForm <- makeFireSenseLCC(neededYear = 2011,
                                    rasterToMatch = RTM,
                                    ## don't pass SA - should still run
                                    destinationPath = dPath,
                                    flammabilityThreshold = 0.25)
  ## if SA is not provided, mask is with RTM, which leads to different NA cell
  expect_true(sum(is.na(flamOtherForm[])) != sum(is.na(flamLCC[])))

  ## cleanup
  withr::deferred_run()
})
