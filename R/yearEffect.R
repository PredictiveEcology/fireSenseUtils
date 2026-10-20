#' The year effect on an expected count: a mean-one lognormal multiplier
#'
#' FireSense gives every simulated year one draw `z ~ N(0, 1)`, shared by spread, ignitions and
#' escapes (drawn in `fireSense_burn`): a bad spread year is also a bad ignition and escape year.
#' Spread adds `z * yearSpreadSD` to its linear predictor (see `yearSpreadSD` in
#' [.objfunSpreadFit()]). An expected count (ignitions or escapes per pixel) is multiplied by
#' `u = exp(sd * z - sd^2 / 2)`, so `E(u) = 1` over `z`: the long-run mean is kept and the result
#' is a valid Poisson mean (`>= 0`). `yearEffectSD()` is the matching estimator of `sd`.
#'
#' @param sd Numeric: the sd of the year effect, one number or one per element of `z`; `NA` is 0.
#' @param z Numeric: the year's draw(s) from `N(0, 1)`.
#'
#' @return `yearEffectMultiplier()`: numeric, `exp(sd * z - sd^2 / 2)`, recycled to the longer of
#'   `sd` and `z`; exactly 1 where `sd` is 0.
#' @export
#' @rdname yearEffect
#' @examples
#' yearEffectMultiplier(0.5, c(-2, 0, 2))
#' mean(yearEffectMultiplier(0.5, rnorm(1e5)))   # about 1
yearEffectMultiplier <- function(sd, z) {
  sd[is.na(sd)] <- 0
  exp(sd * z - sd^2 / 2)
}

#' @description
#' `yearEffectSD()` estimates the sd of the year effect from yearly totals, by the method of moments.
#' Model: `O_j ~ Poisson(mu_j * u_j)` with `u_j = yearEffectMultiplier(sd, z_j)`, so
#' `Var(u_j) = exp(sd^2) - 1` and `E((O_j - mu_j)^2 - O_j) = mu_j^2 * (exp(sd^2) - 1)`: subtracting
#' `O_j` removes the Poisson sampling variance. The estimate is
#' `exp(sd^2) - 1 = sum((O_j - mu_j)^2 - O_j) / sum(mu_j^2)`, floored at 0 (less than Poisson spread
#' means no year effect). `mu_j` is `predicted_j * sum(O) / sum(predicted)`, so a bias in the overall
#' level is not counted as year-to-year variation. It works on counts, not on `log(O_j / mu_j)`, so
#' years with no event (common for escapes) need no special treatment.
#'
#' @param observed,predicted Numeric vectors, one element per year: the observed and the predicted
#'   totals. Years with an `NA` in either are dropped.
#'
#' @return `yearEffectSD()`: `sd`, one number `>= 0`; 0 with fewer than 2 years, or when nothing
#'   was observed or predicted.
#' @export
#' @rdname yearEffect
#' @examples
#' set.seed(1)
#' mu <- rep(c(30, 80, 150, 300), 10)
#' obs <- rpois(length(mu), mu * yearEffectMultiplier(0.5, rnorm(length(mu))))
#' yearEffectSD(obs, mu)                          # about 0.5
#' yearEffectSD(rpois(40, 100), rep(100, 40))     # about 0: Poisson noise only
yearEffectSD <- function(observed, predicted) {
  ok <- !is.na(observed) & !is.na(predicted)
  observed <- observed[ok]
  predicted <- predicted[ok]
  if (length(observed) < 2L || sum(predicted) <= 0 || sum(observed) <= 0) return(0)
  mu <- predicted * sum(observed) / sum(predicted)
  extra <- sum((observed - mu)^2 - observed) / sum(mu^2)
  sqrt(log1p(max(extra, 0)))
}
