#' Handling piecewise terms in a formula
#'
#' @param variable DESCRIPTION NEEDED
#' @param knot DESCRIPTION NEEDED
#' @return DESCRIPTION NEEDED
#'
#' @export
pw <- function(variable, knot) pmax(variable - knot, 0)

#' Order of Magnitude
#'
#' @param x a numeric
#' @return the order of magnitude
#' @export
oom <- function(x) 10^(ceiling(log10(abs(x))))

#' Extract the elements of the special terms, i.e. the variable and the knot value
#'
#' @param v DESCRIPTION NEEDED
#' @param k DESCRIPTION NEEDED
#'
#' @return DESCRIPTION NEEDED
#'
#' @export
extractSpecial <- function(v, k) {
  cl <- match.call()
  if (missing(k)) {
    stop("> Argument 'knotName' is missing (variable '", as.character(cl$v), "')")
  } else {
    list(variable = if (is(cl$v, "AsIs")) cl$v else as.character(cl$v), knot = as.character(cl$k))
  }
}

#' Objective function when no piecewise model is used
#'
#' @param params DESCRIPTION NEEDED
#' @param linkinv the link function
#' @param nll the log-likelihood function
#' @param sm scaling matrix
#' @param nx number of covariates
#' @param mm model matrix containing data
#' @param mod_env the environment containing params - can be a data.frame
#' @param offset DESCRIPTION NEEDED
#'
#' @return DESCRIPTION NEEDED
#'
#' @export
#' @rdname objFunIgnition
.objFunIgnition <- function(params, linkinv, nll, sm, nx, mm, mod_env, offset) {
  ## Parameters scaling
  params <- drop(params %*% sm)

  mu <- drop(mm %*% params[1L:nx]) + offset

  ## link implementation
  mu <- linkinv(mu)

  if (any(mu <= 0) || anyNA(mu) || any(is.infinite(mu)) || length(mu) == 0) {
    return(1e20)
  } else {
    return(eval(nll, envir = mod_env))
  }
}

#' Function to pass to the optimizer - Piece-wise version
#'
#' @param params DESCRIPTION NEEDED
#' @param formula DESCRIPTION NEEDED
#' @param linkinv DESCRIPTION NEEDED
#' @param nll DESCRIPTION NEEDED
#' @param sm DESCRIPTION NEEDED
#' @param updateKnotExpr DESCRIPTION NEEDED
#' @param nx DESCRIPTION NEEDED
#' @param mod_env the environment containing params - can be a data.frame
#' @param offset DESCRIPTION NEEDED
#'
#' @return DESCRIPTION NEEDED
#'
#' @export
#' @importFrom stats model.matrix as.formula
#' @rdname objFunIgnitionPW
.objFunIgnitionPW <- function(params, formula, linkinv, nll, sm, updateKnotExpr, nx, mod_env, offset) {
  ## Parameters scaling
  params <- drop(params %*% sm)
  formula <- as.formula(formula)
  eval(updateKnotExpr, envir = mod_env) ## update knot's values

  mu <- drop(model.matrix(formula, mod_env) %*% params[1:nx]) + offset

  ## link implementation
  mu <- linkinv(mu)

  if (any(mu <= 0) || anyNA(mu) || any(is.infinite(mu)) || length(mu) == 0) {
    return(1e20)
  } else {
    return(eval(nll, envir = mod_env))
  }
}

#' \code{objNlminb}
#'
#' Wrapper around \code{stats::nlminb}
#'
#' @param start DESCRIPTION NEEDED
#' @param objective objective function
#' @param lower lower bounds on coefficients
#' @param upper upper bounds on coefficients
#' @param control DESCRIPTION NEEDED
#' @param hvPW logical indicating whether the formula is piece-wise #IE added
#' @param ... additional arguments passed to objective function
#'
#' @return DESCRIPTION NEEDED
#'
#' @export
#' @importFrom stats nlminb
objNlminb <- function(start, objective, lower, upper, control, hvPW, ...) {
  dots <- list(...)
  nlminb.call <- quote(nlminb(start = start, objective = objective, lower = lower, upper = upper, control = control,
                              linkinv = dots$linkinv, nll = dots$nll, sm = dots$sm, nx = dots$nx,
                              mm = dots$mm, updateKnotExpr = dots$updateKnotExpr, #only one of these needed depending on hvPW
                              mod_env = dots$mod_env, offset = dots$offset, formula = dots$formula))

  if (hvPW){
    nlminb.call$mm <- NULL
  } else {
    nlminb.call$updateKnotExpr <- NULL
  }
  o <- eval(nlminb.call)

  i <- 1L

  ## When there is no convergence and restart is possible, run nlminb() again
  while (as.integer(gsub("[\\(\\)]", "", regmatches(o$message, gregexpr("\\(.*?\\)", o$message))[[1L]])) %in% 7:14 & i < 3L) {
    i <- i + 1L
    o <- eval(nlminb.call)
  }
  o
}
