## Tests for covsX1000AndSetDF
##
## toX1000 is exercised directly in test-helpers-utils.R; here we focus on
## the integration wrapper: column ordering, list shape, paramOrder filtering,
## and that fireBufferedList / fireLociList pass through as data.frames.

library(data.table)

test_that("covsX1000AndSetDF: returns the four expected list elements", {
  annualList    <- list(year2001 = data.table(pixelID = 1:2, cov1 = c(0.1, 0.2)))
  nonAnnualList <- list(`2001_2002` = data.table(pixelID = 1:2, vegPC1 = c(0.3, 0.4)))
  fireBufferedList <- list(year2001 = data.table(pixelID = 1:2, buffer = c(1L, 0L), ids = 1L))
  fireLociList     <- list(year2001 = data.table(pixelID = 1L, size = 5L))
  paramOrder       <- c(cov1 = 1, vegPC1 = 2)

  out <- covsX1000AndSetDF(annualList, nonAnnualList, fireBufferedList,
                           fireLociList, paramOrder, toX1000Integer = TRUE)

  expect_named(out, c("annualDTx1000", "nonAnnualDTx1000",
                      "fireBufferedListDT", "historicalFires"))
})

test_that("covsX1000AndSetDF: annual covariates are multiplied by 1000 and integerised", {
  annualList    <- list(year2001 = data.table(pixelID = 1:2, cov1 = c(0.1, 0.2)))
  nonAnnualList <- list(`2001` = data.table(pixelID = 1:2, vegPC1 = c(0.5, 0.5)))
  fireBufferedList <- list(year2001 = data.table(pixelID = 1L, buffer = 1L, ids = 1L))
  fireLociList     <- list(year2001 = data.table(pixelID = 1L, size = 5L))
  paramOrder       <- c(cov1 = 1, vegPC1 = 2)

  out <- covsX1000AndSetDF(annualList, nonAnnualList, fireBufferedList,
                           fireLociList, paramOrder, toX1000Integer = TRUE)
  expect_equal(out$annualDTx1000[[1]]$cov1, c(100L, 200L))
  expect_equal(out$nonAnnualDTx1000[[1]]$vegPC1, c(500L, 500L))
})

test_that("covsX1000AndSetDF: errors cleanly when toX1000Integer = FALSE (inverse not implemented)", {
  annualList    <- list(year2001 = data.table(pixelID = 1L, cov1 = 0.1))
  nonAnnualList <- list(`2001` = data.table(pixelID = 1L, vegPC1 = 0.5))
  fireBufferedList <- list(year2001 = data.table(pixelID = 1L, buffer = 1L, ids = 1L))
  fireLociList     <- list(year2001 = data.table(pixelID = 1L, size = 5L))
  paramOrder       <- c(cov1 = 1, vegPC1 = 2)

  expect_error(
    covsX1000AndSetDF(annualList, nonAnnualList, fireBufferedList,
                      fireLociList, paramOrder, toX1000Integer = FALSE),
    regexp = "not yet implemented"
  )
})

test_that("covsX1000AndSetDF: paramOrder filters which columns are retained for reordering", {
  ## 'extraCov' is not in paramOrder, so it stays at the end of the table.
  annualList    <- list(year2001 = data.table(pixelID = 1:2, cov1 = c(0.1, 0.2),
                                              extraCov = c(0.5, 0.5)))
  nonAnnualList <- list(`2001` = data.table(pixelID = 1:2, vegPC1 = c(0.5, 0.5)))
  fireBufferedList <- list(year2001 = data.table(pixelID = 1L, buffer = 1L, ids = 1L))
  fireLociList     <- list(year2001 = data.table(pixelID = 1L, size = 5L))
  paramOrder       <- c(cov1 = 1, vegPC1 = 2)

  out <- covsX1000AndSetDF(annualList, nonAnnualList, fireBufferedList,
                           fireLociList, paramOrder, toX1000Integer = TRUE)
  ## First two columns should be pixelID then cov1
  expect_equal(colnames(out$annualDTx1000[[1]])[1:2], c("pixelID", "cov1"))
})
