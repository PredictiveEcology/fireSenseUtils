#' Generate, Summarize, Predict Burn Classes from Covariates
#'
#' @param df A data.frame (or data.table), with covarites, including "burned" (a binary 0, 1
#' notburned = 0, burned = 1), e.g., timeSinceFire, biomassJackPine, etc. that will
#' be used to find fuel classes. This set of covariates must be available both during
#' fitting and for prediction. These must be quantitative.
#' @param numClasses A vector indicating how many classes should be attempted. The function
#' will return the number of classes that best classify the data into homogeneous groups.
#'
#' @return
#' A list with 2 elements, first the \code{model}, which comes from \code{mclust::Mclust},
#' and second the Area Under the Curve or AUC as an indicator of the overal goodness of fit.
#'
#' @rdname burnClass
#' @importFrom mclust Mclust
#' @importFrom pROC roc
#' @importFrom utils modifyList
#' @export
#' @examples
#' #################################
#' # Use own data; here is a generated set for reprex
#' library("data.table")
#' N <- 1e5
#' DT <- data.table(burned = sample(c(0, 0, 0, 1), replace = TRUE, size = N))
#' set(DT, NULL, "jp", rlnorm(N, mean = 4 + 0.5 * DT$burned, sd = 0.25))
#' set(DT, NULL, "bs", rlnorm(N, mean = 4 + 0.3 * DT$burned, sd = 0.25))
#' set(DT, NULL, "ws", rlnorm(N, mean = 4 + 0.2 * DT$burned, sd = 0.25))
#' set(DT, NULL, "age", rlnorm(N, mean = 4 - 0.2* DT$burned, sd = 0.25))
#' DT[, c("jp", "bs", "ws", "age") := lapply(.SD, function(x) x/max(x) * 1000), .SDcols = c("jp", "bs", "ws", "age")]
#' DT[, c("age") := lapply(.SD, function(x) x/max(x) * 200), .SDcols = c("age")]
#' DT2 <- copy(DT)
#' summary(DT)
#' boxplot(DT$age ~ DT$burned)
#'
#'
#' bc <- burnClassGenerator(DT, 4:8)
#' (summ <- burnClassSummary(bc$model))
#' set(DT, NULL, "burnClass", burnClassPredict(bc$model, df = DT))
#' prob <- burnProbFromClass(bc$model, DT)
burnClassGenerator <- function(df, numClasses = 4:9, AUC = TRUE, plotAUC = FALSE) {
  df <- as.data.table(df)
  sam <- sample(NROW(df), NROW(df)*0.75)
  dftrain <- df[sam]
  dftest <- df[-sam]

  mod <- Mclust(dftrain, G = numClasses)

  summ <- burnClassSummary(mod)
  set(dftest, NULL, "burnClass", burnClassPredict(mod, dftest))
  dftest[summ[,c("burnClass", "burned")], burnProb := i.burned, on = "burnClass"]
  rocDF <- data.frame(predictions = dftest$burnProb, labels = dftest$burned)
  args <- list(rocDF$labels, rocDF$predictions, ci=TRUE, ci.alpha=0.9, stratified=FALSE, smoothed = TRUE)
  if (plotAUC) args <- modifyList(args, list(plot=TRUE, auc.polygon=TRUE, max.auc.polygon=TRUE, grid=TRUE,
                                             print.auc=TRUE, show.thres=TRUE))
  pROC_obj <- do.call(roc, args)
  list(model = mod, AUC = pROC_obj$auc)
}

#' @rdname burnClass
#' @export
#' @param mod A model of class \code{Mclust}, e.g,. coming from \code{Mclust} or
#'   \code{burnClassGenerator}
burnClassSummary <- function(mod) {
  df <- as.data.table(mod$data)
  set(df, NULL, "burnClass", mod$classification)
  df[, lapply(.SD, mean), by = "burnClass"]
}

#' @rdname burnClass
#' @export
burnClassPredict <- function(mod, df) {
  pred <- predict(mod, newdata = df)
  pred$classification
}

#' @rdname burnClass
#' @export
burnProbFromClass <- function(mod, df) {
  df <- as.data.table(df)
  summ <- burnClassSummary(mod)
  colsToUse <- intersect(colnames(mod$data), colnames(df))
  set(df, NULL, "burnClass", burnClassPredict(mod, df[, ..colsToUse]))
  df[summ[,c("burnClass", "burned")], burnProb := i.burned, on = "burnClass"]
  df$burnProb
}

