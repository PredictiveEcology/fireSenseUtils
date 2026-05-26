#' Handling piecewise terms in a formula
#'
#' Implements the right-hand half of a piecewise-linear (hinge) term:
#' `pmax(variable - knot, 0)`. Used in formulas passed to the ignition
#' fitting routines so that a covariate can have zero effect below a
#' threshold and a linear effect above it.
#'
#' @param variable Numeric vector of covariate values.
#' @param knot Numeric scalar threshold (knot) below which the term is `0`.
#' @return Numeric vector the same length as `variable`, with each element
#'   equal to `max(variable[i] - knot, 0)`.
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
#' Used by the ignition formula machinery to parse a `pw(v, k)` term and
#' return its components without evaluating them. The variable and knot
#' names are returned as character strings (the captured symbols), or kept
#' as `AsIs` if `v` was wrapped with `I()`.
#'
#' @param v The variable symbol or expression appearing inside `pw()`.
#' @param k The knot value or symbol; **must** be supplied.
#'
#' @return A list with elements `variable` (character, or `AsIs` expression)
#'   and `knot` (character form of `k`). Errors if `k` is missing.
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
#' Internal callback passed to the optimizer (e.g. [stats::nlminb()]) when
#' the ignition formula contains no [pw()] terms. Computes the linear
#' predictor `mm %*% params + offset`, applies `linkinv`, and returns the
#' negative log-likelihood evaluated in `mod_env`. A large sentinel
#' (`1e20`) is returned whenever the link inverse produces invalid values
#' (non-positive, `NA`, infinite, or empty), preventing the optimizer from
#' wandering into infeasible regions.
#'
#' @param params Numeric vector of parameter values supplied by the optimizer.
#' @param linkinv Function. Inverse link function (e.g. [stats::plogis()]
#'   for a logit link).
#' @param nll Quoted/unevaluated R expression that computes the negative
#'   log-likelihood when evaluated in `mod_env`.
#' @param sm Square numeric matrix. Parameter scaling matrix applied as
#'   `params %*% sm` before the linear predictor is formed.
#' @param nx Integer. Number of covariate parameters (i.e. number of
#'   columns in `mm`).
#' @param mm Numeric model matrix containing covariate values.
#' @param mod_env An environment (or data-frame coerced via [base::list2env()])
#'   containing the variables referenced by `nll`.
#' @param offset Numeric vector the same length as `nrow(mm)`. Model offset
#'   added to the linear predictor.
#'
#' @return Numeric scalar: the negative log-likelihood, or `1e20` if the
#'   inverse-link produces non-finite or non-positive values.
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
#' Internal callback for the optimizer when the ignition formula contains
#' [pw()] terms. Differs from [.objFunIgnition] in that the model matrix is
#' rebuilt at every call (via [stats::model.matrix()]) after running
#' `updateKnotExpr`, which updates the knot values in `mod_env`. This lets
#' the optimizer jointly fit covariate coefficients and knot positions.
#'
#' @param params Numeric vector of parameter values supplied by the optimizer.
#'   Includes both covariate coefficients and knot values.
#' @param formula Model formula (or string coercible by [stats::as.formula()])
#'   containing one or more `pw()` calls.
#' @param linkinv Function. Inverse link function.
#' @param nll Quoted R expression evaluating to the negative log-likelihood
#'   when evaluated in `mod_env`.
#' @param sm Square numeric matrix. Parameter scaling matrix.
#' @param updateKnotExpr Quoted R expression that copies the latest knot
#'   values from `params` into the variables in `mod_env` referenced by
#'   `formula`.
#' @param nx Integer. Number of covariate parameters.
#' @param mod_env An environment (or data-frame) containing the covariates
#'   and the variables updated by `updateKnotExpr`.
#' @param offset Numeric vector. Model offset added to the linear predictor.
#'
#' @return Numeric scalar: the negative log-likelihood, or `1e20` if the
#'   inverse-link produces negative, `NA`, or infinite values.
#'
#' @export
#' @importFrom stats model.matrix as.formula
#' @rdname objFunIgnitionPW
.objFunIgnitionPW <- function(params, formula,
                              linkinv, nll, sm, updateKnotExpr, nx, mod_env, offset) {
  ## Parameters scaling

  params <- drop(params %*% sm)
  formula <- as.formula(formula)

  eval(updateKnotExpr) ## update knot's values

  mm <- model.matrix(formula, mod_env)
  mu <- drop(mm %*% params[1:nx]) + offset

  ## link implementation
  mu <- linkinv(mu)

  if (any(mu < 0) || anyNA(mu) || any(is.infinite(mu)) || length(mu) == 0) {
    return(1e20)
  } else {
    return(eval(nll, envir = mod_env))
  }
}

#' `objNlminb`
#'
#' Wrapper around [stats::nlminb()] that calls one of the ignition objective
#' functions (`.objFunIgnition` or `.objFunIgnitionPW`, dispatched by `hvPW`)
#' and falls back to [stats::optim()] (L-BFGS-B) when `nlminb` reports a
#' "false convergence" or similar non-convergence code (codes 7-14).
#'
#' @param x Numeric vector. Starting parameter values for the optimizer.
#' @param objective Function. Objective function to minimize, typically
#'   `.objFunIgnition` or `.objFunIgnitionPW`.
#' @param lower Numeric vector. Lower bounds on each parameter, same length
#'   as `x`.
#' @param upper Numeric vector. Upper bounds on each parameter, same length
#'   as `x`.
#' @param control List of control options passed to [stats::nlminb()] (or
#'   [stats::optim()] on fallback, with `iter.max` and `eval.max` dropped).
#' @param hvPW Logical. If `TRUE`, dispatches to the piecewise objective
#'   (formula is rebuilt per call); if `FALSE`, dispatches to the
#'   non-piecewise objective (pre-built `mm` is supplied).
#' @param ... Additional arguments passed through to `objective`. Expected
#'   names are: `linkinv`, `nll`, `sm`, `nx`, `mm`, `updateKnotExpr`,
#'   `mod_env`, `offset`, `formula`.
#'
#' @return The list returned by [stats::nlminb()] on success, or the result
#'   of [stats::optim()] (with `value` renamed to `objective`) on fallback.
#'
#' @export
#' @importFrom stats nlminb
objNlminb <- function(x, objective, lower, upper, control, hvPW, ...) {
  dots <- list(...)

  controlOptim <- control[setdiff(names(control), c("iter.max", "eval.max"))]
  optim.call <- quote(optim(
    par = x, fn = objective, lower = lower, upper = upper, control = controlOptim,
    linkinv = dots$linkinv, nll = dots$nll, sm = dots$sm, nx = dots$nx,
    mm = dots$mm, updateKnotExpr = dots$updateKnotExpr, # only one of these needed depending on hvPW
    mod_env = dots$mod_env, method = "L-BFGS-B",
    offset = dots$offset, formula = dots$formula
  ))

  nlminb.call <- quote(nlminb(
    start = x, objective = objective, lower = lower, upper = upper, control = control,
    linkinv = dots$linkinv, nll = dots$nll, sm = dots$sm, nx = dots$nx,
    mm = dots$mm, updateKnotExpr = dots$updateKnotExpr, # only one of these needed depending on hvPW
    mod_env = dots$mod_env,
    offset = dots$offset, formula = dots$formula
  ))

  if (hvPW) {
    nlminb.call$mm <- NULL
    optim.call$mm <- NULL
  } else {
    nlminb.call$formula <- NULL
    optim.call$formula <- NULL
    nlminb.call$updateKnotExpr <- NULL
    optim.call$updateKnotExpr <- NULL
  }

  o <- eval(nlminb.call)

  # i <- 1L

  ## When there is no convergence and restart is possible, run nlminb() again (Jean)
  ## Eliot commented this out, because fitting seems to be fairly deterministic, so restarting
  ##   with same starts gives exactly same answer
  if (as.integer(gsub("[\\(\\)]", "", regmatches(o$message, gregexpr("\\(.*?\\)", o$message))[[1L]])) %in% 7:14) {
    #   i <- i + 1L
    message("nlminb did not converge; trying optim")
    o2 <- try(eval(optim.call))
    if (!is(o2, "try-error")) {
      o <- o2
      o$objective <- o$value
      o$value <- NULL
    }
    # break
  }
  o
}
