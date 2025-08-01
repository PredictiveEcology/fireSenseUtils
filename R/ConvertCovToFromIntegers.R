#' Convert covariates to/from integers for smaller footprint; make `data.frame`
#'
#' @param annualList Named list of `data.table`s of covariates that vary annually.
#'   Must have a `pixelID` column. List name is year.
#' @param nonAnnualList Named list of `data.table`s of covariates that do not vary each year.
#'   Must have a `pixelID` column. List name is a concatenated set of years separated by underscore.
#' @param fireBufferedList Named list of `data.table`s of buffered perimeters with `pixelID`.
#'   Must have a `pixelID` column. List name is `yearXXXX` where `XXXX` is the 4 digit year.
#'   Column names must be `pixelID`, `buffer` with every row a `1` or `0` as to whether that pixel
#'   is in or outside the buffer (`1` is in the buffer!), and `ids` with a unique identifier
#'   for each fire.
#' @param fireLociList Named list of `data.table` of `pixelID` columns where fires occurred.
#'   Columns are `size` which is number of pixels in the fire
#' @param paramOrder Character vector of the order of covariates that will be used.
#' @param toX1000Integer Logical. If `TRUE`, then it multiples by 1000 and coerces to integer.
#'   If `FALSE`, then it does the inverse.
#'
#' @return a data.table with relevant columns made mutually exclusive
#'
#' @export
#' @importFrom data.table setDF
#'
covsX1000AndSetDF <- function(annualList, nonAnnualList, fireBufferedList, fireLociList,
                              paramOrder, toX1000Integer = TRUE) {

  annualCols <- colnames(annualList[[1]])
  nonAnnualCols <- colnames(nonAnnualList[[1]])
  annualCols <- annualCols[annualCols %in% names(paramOrder)]
  nonAnnualCols <- nonAnnualCols[nonAnnualCols %in% names(paramOrder)]

  annualList <- lapply(annualList, setcolorder, neworder = c("pixelID", annualCols))
  nonAnnualCols <- lapply(nonAnnualList, setcolorder, neworder = c("pixelID", nonAnnualCols))

  annualDT <- lapply(annualList, setDF)
  nonAnnualDT <- lapply(nonAnnualList, setDF)
  if (isTRUE(toX1000Integer)) {
    annualDTx1000 <- toX1000(annualDT)
    nonAnnualDTx1000 <- toX1000(nonAnnualDT)
  } else {
    browser()
  }
  fireBufferedListDT <- lapply(fireBufferedList, setDF)
  historicalFires <- lapply(fireLociList, setDF)
  list(annualDTx1000 = annualDTx1000,
       nonAnnualDTx1000 = nonAnnualDTx1000,
       fireBufferedListDT = fireBufferedListDT,
       historicalFires = historicalFires)
}

#' Convert numeric values to integers x 1000
#'
#' This simply converts to integers times 1000 so that the values can be saved
#' more quickly to disk, for example.
#'
#' @param lst a data.frame of (at least some) numeric columns
#'
#' @param omitCols character A vector of column names to not convert
#'
#' @return The same data.frame, but with columns converted to integer and times 1000
#'
#' @export
#' @importFrom data.table setDF setDT
#' @importFrom LandR asInteger
#'
toX1000 <- function(lst, omitCols = "pixelID") {
  annualDTx1000 <- lapply(lst, function(dt) {
    setDT(dt)
    cns <- setdiff(colnames(dt), omitCols)
    for (colnam in cns)
      set(dt, NULL, colnam, asInteger(dt[[colnam]] * 1000))
    setDF(dt)
  })
}
