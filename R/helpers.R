#' Four- and five-parameter logistic functions
#'
#' @param x DESCRIPTION NEEDED
#' @param par DESCRIPTION NEEDED
#'
#' @return DESCRIPTION NEEDED
#'
#' @export
#' @rdname logistic
logistic4p <- function(x, par) {
  par[1L] + (par[2L] - par[1L]) / (1 + x^(-par[3L]))^par[4L]
}

#' @export
#' @rdname logistic
logistic5p <- function(x, par) {
  par[1L] + (par[2L] - par[1L]) / (1 + (x/par[3L])^(-par[4L])) ^ par[5L]
}


#' @export
#' @rdname logistic
logistic3p <- function(x, par, par1 = 0.1) {
  par1 + (par[1L] - par1) / (1 + x^(-par[2L]))^par[3L]
}

#' Replace \code{NA}s in a \code{data.table} with zeros
#'
#' @param DT DESCRIPTION NEEDED
#' @param colsToUse DESCRIPTION NEEDED
#'
#' @return DESCRIPTION NEEDED
#'
#' @export
#' @importFrom data.table set
dtReplaceNAwith0 <- function(DT, colsToUse = NULL) {
  if (is.null(colsToUse)) {
    colsToUse <- names(DT)
  }
  for (i in colsToUse) {
    nas <- which(is.na(DT[[i]]))
    if (length(nas)) {
      set(DT, nas, i, 0)
    }
  }
  DT
}

#' Convert annual raster stacks to \code{data.table}
#'
#' @param annualStacks DESCRIPTION NEEDED
#' @param whNotNA DESCRIPTION NEEDED
#' @param ... DESCRIPTION NEEDED
#'
#' @return DESCRIPTION NEEDED
#'
#' @export
#' @importFrom data.table as.data.table set
#' @importFrom LandR asInteger
annualStacksToDTx1000 <- function(annualStacks, whNotNA, ...) {
  # whNotNA <- which(!is.na(rasterToMatch[]))
  rastersDT <- lapply(annualStacks, whNotNA = whNotNA, function(x, whNotNA) {
    a <- as.data.table(x[])[whNotNA]
    a <- dtReplaceNAwith0(a)
    a
  })
  lapply(rastersDT, function(x) {
    for (col in colnames(x)) {
      set(x, NULL, col, asInteger(x[[col]]*1000))
    }
  })

  rastersDT
}

#' Generate random beta variates between 2 values and a mean
#' @inheritParams stats::rbeta
#' @param shape2 If provided, passed to rbeta. If not, \code{m} must be (i.e., the mean)
#' @param l scalar numeric for the lower bound
#' @param u scalar numeric for the upper bound
#' @param m scalar numeric for the mean
#'
#' @export
rbetaBetween <- function(n, l, u, m, shape1, shape2 = NULL) {
  if (is.null(shape2)) {
    m1 <- ((1)/(u - l) * (m - l))
    shape2 <- (shape1 - shape1*m1)/m1
  }
  out <- rbeta(n, shape1, shape2)
  out * (u - l) + (l)
}
