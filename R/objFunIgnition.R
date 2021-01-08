
#' Handling piecewise terms in a formula
#' @param variable DESCRIPTION NEEDED
#' @param knot DESCRIPTION NEEDED
#' @return DESCRIPTION NEEDED

#' @export
pw <- function(variable, knot) pmax(variable - knot, 0)


#' order of magnitude
#'
#' @param x DESCRIPTION NEEDED
#' @return DESCRIPTION NEEDED
#' @export
oom <- function(x) 10^(ceiling(log10(abs(x))))

#' Extract the elements of the special terms, i.e. the variable and the knot value
#' @param v DESCRIPTION NEEDED
#' @param k DESCRIPTION NEEDED
#' @return DESCRIPTION NEEDED
#' @export
extractSpecial <- function(v, k){
  cl <- match.call()

  if(missing(k)) stop("> Argument 'knotName' is missing (variable '", as.character(cl$v), "')")
  else list(variable = if (is(cl$v, "AsIs")) cl$v else as.character(cl$v), knot = as.character(cl$k))
}

#' Objective function when no piece-wise model is used
#' @param params DESCRIPTION NEEDED
#' @param linkinv DESCRIPTION NEEDED
#' @param nll DESCRIPTION NEEDED
#' @param sm DESCRIPTION NEEDED
#' @param nx DESCRIPTION NEEDED
#' @param mm DESCRIPTION NEEDED
#' @param mod_env DESCRIPTION NEEDED
#' @param offset DESCRIPTION NEEDED
#' @return DESCRIPTION NEEDED
objIgnition <- function(params, linkinv, nll, sm, nx, mm, mod_env, offset) {
  ## Parameters scaling
  params <- drop(params %*% sm)

  mu <- drop(mm %*% params[1L:nx]) + offset

  ## link implementation
  mu <- linkinv(mu)

  if(any(mu <= 0) || anyNA(mu) || any(is.infinite(mu)) || length(mu) == 0) return(1e20)
  else return(eval(nll, envir = as.list(environment()), enclos = mod_env))
}


#' Function to pass to the optimizer - Piece-wise version
#' @param params DESCRIPTION NEEDED
#' @param formula DESCRIPTION NEEDED
#' @param linkinv DESCRIPTION NEEDED
#' @param nll DESCRIPTION NEEDED
#' @param sm DESCRIPTION NEEDED
#' @param updateKnotExpr DESCRIPTION NEEDED
#' @param nx DESCRIPTION NEEDED
#' @param mod_env DESCRIPTION NEEDED
#' @param offset DESCRIPTION NEEDED
#'
#' @importFrom stats model.matrix
#' @return DESCRIPTION NEEDED
objIgnitionPW <- function(params, formula, linkinv, nll, sm, updateKnotExpr, nx, mod_env, offset) {
  ## Parameters scaling
  params <- drop(params %*% sm)

  eval(updateKnotExpr, envir = environment(), enclos = mod_env) ## update knot's values

  mu <- drop(model.matrix(formula, mod_env) %*% params[1:nx]) + offset

  ## link implementation
  mu <- linkinv(mu)

  if(any(mu <= 0) || anyNA(mu) || any(is.infinite(mu)) || length(mu) == 0) return(1e20)
  else return(eval(nll, envir = as.list(environment()), enclos = mod_env))
}


#' Nlminb wrapper
#' @param start DESCRIPTION NEEDED
#' @param objective DESCRIPTION NEEDED
#' @param lower DESCRIPTION NEEDED
#' @param upper DESCRIPTION NEEDED
#' @param control DESCRIPTION NEEDED
#' @return DESCRIPTION NEEDED
objNlminb <- function(start, objective, lower, upper, control) {
  nlminb.call <- quote(nlminb(start = start, objective = objective, lower = lower, upper = upper, control = control))
  nlminb.call[names(formals(objective)[-1L])] <- parse(text = formalArgs(objective)[-1L])

  o <- eval(nlminb.call)

  i <- 1L

  ## When there is no convergence and restart is possible, run nlminb() again
  while(as.integer(gsub("[\\(\\)]", "", regmatches(o$message, gregexpr("\\(.*?\\)", o$message))[[1L]])) %in% 7:14 & i < 3L)
  {
    i <- i + 1L
    o <- eval(nlminb.call)
  }
  o
}
