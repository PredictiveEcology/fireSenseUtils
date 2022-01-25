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
  par[1L] + (par[2L] - par[1L]) / (1 + (exp(x) / par[3L])^(-par[4L]))^par[5L]
}

#' @param par1 DESCRIPTION NEEDED
#' @export
#' @rdname logistic
logistic3p <- function(x, par, par1 = 0.1) {
  par1 + (par[1L] - par1) / (1 + exp(x)^(-par[2L]))^par[3L]
}

#' @param par4 DESCRIPTION NEEDED
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

#' Convert list of annual raster stack to \code{data.table}
#'
#' @param x \code{RasterStack} or list of rasters to convert to \code{data.table}
#'          and multiply by 1000 to save space
#' @param whNotNA Pixel indexes that should go through this process (i.e. not NA)
#' @param ... Not currently used
#'
#' @return \code{data.table} of the raster stack or the list
#' @export
#' @rdname annualStackToDTx1000
#'
#' @examples
#' library(raster)
#' r1 <- raster(extent(0, 10, 0, 10), vals = 1:100)
#' r2 <- raster(extent(0, 10, 0, 10), vals = 100:1)
#' r3 <- raster(extent(0, 10, 0, 10), vals = 200:101)
#' r4 <- raster(extent(0, 10, 0, 10), vals = 300:201)
#'
#' # list of Rasters
#' lRast <- list(r1, r2, r3)
#' lRast[[1]][5] <- NA
#' whNotNA <- setdiff(1:ncell(r1), 5)
#'
#' # unnamed -- should error
#' try(out1 <- annualStackToDTx1000(lRast, whNotNA))
#'
#' # named
#' names(lRast) <- c("OneToHun", "HunToOne", "TwoHunToOneHun")
#' out1 <- annualStackToDTx1000(lRast, whNotNA)
#'
#' # RasterStack
#' out2 <- annualStackToDTx1000(raster::stack(lRast), whNotNA)
#'
#' # List of RasterStacks
#' s1 <- raster::stack(r1, r2)
#' names(s1) <- names(lRast)[1:2]
#' s2 <- raster::stack(r4, r3)
#' names(s2) <- c(names(lRast)[3], "ThreeHunToTwoHun")
#' out3 <- annualStackToDTx1000(list(s1 = s1, s2 = s2), whNotNA) ## named list required
#'
#' # With duplicated names -- to remove duplicates;
#' #  actually, this doesn't make sense: RasterStack can't have duplicated names
#' names(lRast) <- c("OneToHun", "OneToHun", "TwoHunToOneHun")
#' out4 <- annualStackToDTx1000(raster::stack(lRast), whNotNA)
annualStackToDTx1000 <- function(x, whNotNA, ...) {
  UseMethod("annualStackToDTx1000")
}

#' @export
#' @importFrom data.table as.data.table set
#' @rdname annualStackToDTx1000
annualStackToDTx1000.RasterLayer <- function(x, whNotNA, ...) {
  layDT <- as.data.table(x[])[whNotNA]
  layDT <- dtReplaceNAwith0(layDT)
  set(layDT, NULL, 1L, asInteger(layDT[[1L]] * 1000))
  names(layDT) <- names(x)
  message("Layer ", names(layDT), " converted to data.table")
  layDT
}

#' @export
#' @importFrom data.table as.data.table
#' @rdname annualStackToDTx1000
annualStackToDTx1000.list <- function(x, whNotNA, ...) {
  # check for names
  # check for rasters
  if (is.null(names(x))) {
    stop("x must be a named list (or stack)")
  }
  out <- lapply(x, whNotNA = whNotNA, annualStackToDTx1000, ...)
  rastersDT <- as.data.table(out)
  # rastersDT <- lapply(names(x), whNotNA = whNotNA, function(x, whNotNA) {
  #   if (any(is(x, "list"), is(x, "RasterStack"))) {
  #     lay <- x[[x]]
  #     } else {
  #       stop("x must be either a list or a RasterStack")
  #     }
  #   layDT <- as.data.table(lay[])[whNotNA]
  #   layDT <- dtReplaceNAwith0(layDT)
  #   names(layDT) <- names(lay)
  #   message("Layer ", names(layDT), " converted to data.table")
  #   return(layDT)
  # })
  # if (is(x, "RasterStack")) { # Should be the prediction, raster stack
  #   rastersDT <- cbindFromList(rastersDT)
  #   rastersDT[ , (names(rastersDT)) := lapply(X = .SD, FUN = function(column){
  #     column <- asInteger(column*1000)
  #     return(column)
  #   }), .SDcols = names(rastersDT)]
  # } else {
  #   lapply(rastersDT, function(x) { # Should be the fitting, list of years. Doesn't work for stack
  #     for (col in colnames(x)) {
  #       set(x, NULL, col, asInteger(x[[col]]*1000))
  #       message("Layer ", col, " converted to integer")
  #     }
  #   })
  # }
  return(rastersDT)
}

#' @export
#' @importFrom data.table as.data.table set
#' @importFrom raster stack unstack
#' @rdname annualStackToDTx1000
annualStackToDTx1000.RasterStack <- function(x, whNotNA, ...) {
  annualList <- raster::unstack(x)
  names(annualList) <- names(x)
  rastersDT <- annualStackToDTx1000(annualList, whNotNA, ...)
  # stkName <- names(x)
  # rastersDT <- lapply(names(x), whNotNA = whNotNA, function(x, whNotNA) {
  #   if (any(is(x, "list"), is(x, "RasterStack"))) {
  #     lay <- x[[x]]
  #   } else {
  #     stop("x must be either a list or a RasterStack")
  #   }
  #   layDT <- as.data.table(lay[])[whNotNA]
  #   layDT <- dtReplaceNAwith0(layDT)
  #   names(layDT) <- names(lay)
  #   return(layDT)
  # })
  # lapply(rastersDT, function(x) { # Should be the fitting and predicting
  #   for (col in colnames(x)) {
  #     set(x, NULL, col, asInteger(x[[col]]*1000))
  #     message("Layer ", col, " converted to integer")
  #   }
  # })
  # if (is(x, "RasterStack")) { # Should be the prediction, raster stack
  #   bindedList <- do.call(cbind, args = rastersDT)
  #   bindedList <- data.table::data.table(bindedList)
  #   if (any(duplicated(names(bindedList)))) {
  #     warning("There are duplicated names in the raster layer. Deleting...",
  #             immediate. = TRUE)
  #     bindedList[, `:=`(which(duplicated(names(bindedList))), NULL)]
  #   }
  #   rastersDT <- bindedList
  # }
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
    m1 <- ((1) / (u - l) * (m - l))
    shape2 <- (shape1 - shape1 * m1) / m1
  }
  out <- rbeta(n, shape1, shape2)
  out * (u - l) + (l)
}

asInteger <- function(x) as.integer(floor(x + 0.5))
