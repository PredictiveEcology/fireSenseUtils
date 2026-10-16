#' The data year a year's vegetation comes from
#'
#' The one rule for matching a year to a vegetation or land-cover data year: the data year closest
#' to, and not after, the year. `fireSense_dataPrepFit` assigns each fire year its data year with it
#' and decides a fit's model structure on it; `fireSense_dataPrepPredict` takes the vegetation of
#' `time(sim)` with it; so a fit and every prediction, hindcast or forecast, read the same data for
#' the same year.
#'
#' @param dataYears integer years with data, in any order.
#' @param years integer years.
#'
#' @return `dataYearFor()`: integer, one data year per element of `years`. Stops if a year precedes
#'   the first data year: it has no vegetation data.
#' @export
#' @rdname dataYearFor
#' @examples
#' dataYearFor(c(1985, 1990, 2000, 2010, 2020), c(1985, 1999, 2015, 2024))
dataYearFor <- function(dataYears, years) {
  dataYears <- sort(as.integer(dataYears))
  if (any(years < min(dataYears))) {
    stop("years before the first dataYear (", min(dataYears), ") have no vegetation data: ",
         paste(years[years < min(dataYears)], collapse = ", "))
  }
  dataYears[findInterval(years, dataYears)]
}

#' @param dataYears integer years with data, in any order.
#' @param years integer years.
#'
#' @return `dataYearBefore()`: integer, the latest data year strictly before each element of `years`.
#'   Stops if a year is not after the first data year.
#' @details `dataYearBefore()` is the rule of a fit: SCANFI vegetation of year Y already shows year
#'   Y's burns, so `fireSense_dataPrepFit` builds a fire year's vegetation table from the data year
#'   strictly before it, and the spread objective must read the same table. Prediction and forecast
#'   keep `dataYearFor()` (not after).
#' @export
#' @rdname dataYearFor
#' @examples
#' dataYearBefore(c(1985, 1990, 2000, 2010, 2020), c(1986, 1990, 1991, 2020))
dataYearBefore <- function(dataYears, years) {
  dataYears <- sort(as.integer(dataYears))
  if (anyNA(dataYears) || any(years <= min(dataYears))) {
    stop("years not after the first dataYear (", min(dataYears), ") have no vegetation data: ",
         paste(years[years <= min(dataYears)], collapse = ", "))
  }
  dataYears[findInterval(years, dataYears, left.open = TRUE)]
}

#' @param x a list (anything with names) with one element per data year, named `year<dataYear>`
#'   (`yearTxt`) as `fireSense_dataPrepFit` names `cohortDatas`, `rstLCCs`, ...; any name whose
#'   digits are the year will do.
#' @param year one integer year.
#'
#' @return `dataYearElement()`: the name of the element of `x` that `year` uses, the one of
#'   `dataYearFor(<the years in names(x)>, year)`.
#' @export
#' @rdname dataYearFor
#' @examples
#' dataYearElement(list(year2010 = "a", year2020 = "b"), 2015)
dataYearElement <- function(x, year) {
  dataYears <- as.integer(gsub("[^0-9]", "", names(x)))
  names(x)[match(dataYearFor(dataYears, year), dataYears)]
}
