#' Predictions from ignition model
#'
#' @param model the model used to fit ingition data - must have a predict method
#' @param data data for prediction
#' @param rescaleFactor spatial rescaling factor when predicted and fitted data are at different scales.
#'     Calculated as: (predResolution/fitResolution)^2
#' @param lambdaRescaleFactor  If the data for fitting has been sampled for pseudo-absences,
#'    this imposes a new baseline probability of fire occurrences, hence predictions need to be adjusted.
#'    If the  original fire prob. is (total no. fires)/(total no. fires + total no. absences), and the fire
#'    probability imposed by sampling is (total no. fires)/(total no. fires + no. sampled pseudo-absences), to
#'    adjust predicted values, one needs to multiply them by (total no. fires + no. sampled pseudo-absences/(total no. fires + total no. absences)
#'
#' @return vector of predicted values.
#'
#' @export

predictIgnition <- function(model, data, rescaleFactor, lambdaRescaleFactor) {

  ignitions <- predict(model, newdata = data, se.fit = FALSE, re.form = NA,
                       type = "response")

  ignitions <- ignitions * rescaleFactor * lambdaRescaleFactor

  #old method
  # mm <- model.matrix(model, data)
  # pred <- mm %*% coefs
  # pred <- drop(pred)
  # pred <- linkinv(pred)
  # pred <- pred * rescaleFactor
  # pred * lambdaRescaleFactor

  return(ignitions)
}
