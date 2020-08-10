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
  par[1L] + (par[2L] - par[1L]) / (1 + exp(x)^(-par[3L]))^par[4L]
}

#' @export
#' @rdname logistic
logistic5p <- function(x, par) {
  par[1L] + (par[2L] - par[1L]) / (1 + (exp(x)/par[3L])^(-par[4L])) ^ par[5L]
}


#' @param par1 DESCRIPTION NEEDED
#' @export
#' @rdname logistic
logistic3p <- function(x, par, par1 = 0.1) {
  par1 + (par[1L] - par1) / (1 + exp(x)^(-par[2L]))^par[3L]
}

#' @export
#' @rdname logistic
logistic2p <- function(x, par, par1 = 0.1, par4 = 0.5) {
  par1 + (par[1L] - par1) / (1 + exp(x)^(-par[2L]))^par4
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
#' @importFrom raster stack
#' @importFrom usefulFuns cbindFromList
annualStacksToDTx1000 <- function(annualStacks, whNotNA, ...) {
  stkName <- names(annualStacks)
  rastersDT <- lapply(names(annualStacks), whNotNA = whNotNA, function(x, whNotNA) {
    if (any(is(annualStacks, "list"),
            is(annualStacks, "RasterStack"))){
      lay <- annualStacks[[x]]
      } else {
        stop("annualStacks must be either a list or a RasterStack")
      }
    layDT <- as.data.table(lay[])[whNotNA]
    layDT <- dtReplaceNAwith0(layDT)
    names(layDT) <- names(lay)
    message("Layer ", names(layDT), " converted to data.table")
    return(layDT)
  })
  if (is(annualStacks, "RasterStack")){ # Should be the prediction, raster stack
    rastersDT <- cbindFromList(rastersDT)
    rastersDT[ , (names(rastersDT)) := lapply(X = .SD, FUN = function(column){
      column <- asInteger(column*1000)
      return(column)
    }), .SDcols = names(rastersDT)]
  } else {
    lapply(rastersDT, function(x) { # Should be the fitting, list of years. Doesn't work for stack
      for (col in colnames(x)) {
        set(x, NULL, col, asInteger(x[[col]]*1000))
        message("Layer ", col, " converted to integer")
      }
    })
  }
  return(rastersDT)
}

#' Convert annual raster stack to \code{data.table}
#'
#' @param annualStack DESCRIPTION NEEDED
#' @param whNotNA DESCRIPTION NEEDED
#' @param ... DESCRIPTION NEEDED
#'
#' @return DESCRIPTION NEEDED
#'
#' @export
#' @importFrom data.table as.data.table set
#' @importFrom LandR asInteger
#' @importFrom raster stack
#' @importFrom usefulFuns cbindFromList
annualStackToDTx1000 <- function(annualStack, whNotNA, ...) {
  stkName <- names(annualStack)
  rastersDT <- lapply(names(annualStack), whNotNA = whNotNA, function(x, whNotNA) {
    if (any(is(annualStack, "list"),
            is(annualStack, "RasterStack"))){
      lay <- annualStack[[x]]
    } else {
      stop("annualStack must be either a list or a RasterStack")
    }
    layDT <- as.data.table(lay[])[whNotNA]
    layDT <- dtReplaceNAwith0(layDT)
    names(layDT) <- names(lay)
    return(layDT)
  })
  lapply(rastersDT, function(x) { # Should be the fitting and predicting
    for (col in colnames(x)) {
      set(x, NULL, col, asInteger(x[[col]]*1000))
      message("Layer ", col, " converted to integer")
    }
  })
  if (is(annualStack, "RasterStack")){ # Should be the prediction, raster stack
    bindedList <- do.call(cbind, args = rastersDT)
    bindedList <- data.table::data.table(bindedList)
    if (any(duplicated(names(bindedList)))){
      warning("There are duplicated names in the raster layer. Deleting...", 
              immediate. = TRUE)
      bindedList[, `:=`(which(duplicated(names(bindedList))), NULL)]
    }
    rastersDT <- bindedList
  }
  return(rastersDT)
}

#' Generate random beta variates between 2 values and a mean
#' @inheritParams stats::rbeta
#' @param shape2 If provided, passed to \code{rbeta}. If not, \code{m} must be (i.e., the mean)
#' @param l scalar numeric for the lower bound
#' @param u scalar numeric for the upper bound
#' @param m scalar numeric for the mean
#'
#' @export
#' @importFrom stats rbeta
rbetaBetween <- function(n, l, u, m, shape1, shape2 = NULL) {
  if (is.null(shape2)) {
    m1 <- ((1)/(u - l) * (m - l))
    shape2 <- (shape1 - shape1*m1)/m1
  }
  out <- rbeta(n, shape1, shape2)
  out * (u - l) + (l)
}
