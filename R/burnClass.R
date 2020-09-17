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
#' @importFrom mclust Mclust mclustBIC
#' @importFrom pROC roc
#' @importFrom utils modifyList
#' @export
#' @examples
#' \dontrun{
#' #################################
#' # Use own data; here is a generated set for reprex
#' library("data.table")
#' N <- 1e5
#' DT <- list()
#' for (i in c("train", "test")) {
#'   DT[[i]] <- data.table(burned = sample(c(0, 0, 0, 1), replace = TRUE, size = N))
#'   set(DT[[i]], NULL, "jp", rlnorm(N, mean = 4 + 0.5 * DT[[i]]$burned, sd = 0.25))
#'   set(DT[[i]], NULL, "bs", rlnorm(N, mean = 4 + 0.3 * DT[[i]]$burned, sd = 0.25))
#'   set(DT[[i]], NULL, "ws", rlnorm(N, mean = 4 + 0.2 * DT[[i]]$burned, sd = 0.25))
#'   set(DT[[i]], NULL, "age", rlnorm(N, mean = 4 - 0.2* DT[[i]]$burned, sd = 0.25))
#'   DT[[i]][, c("jp", "bs", "ws", "age") := lapply(.SD, function(x) x/max(x) * 1000), .SDcols = c("jp", "bs", "ws", "age")]
#'   DT[[i]][, c("age") := lapply(.SD, function(x) x/max(x) * 200), .SDcols = c("age")]
#'   summary(DT[[i]])
#'   boxplot(DT[[i]]$age ~ DT[[i]]$burned)
#' }
#'
#' bc <- burnClassGenerator(DT[["train"]], 4:8)
#' # Show if the model is good at predicting burn state
#' (bc$AUC) # area under the curve
#'
#' # print summary of mean values of each burn class
#' (summ <- burnClassSummary(bc$model))
#'
#' # predict -- add Burn Class to object
#' set(DT[["test"]], NULL, "burnClass", burnClassPredict(bc$model, df = DT[["test"]]))
#' prob <- burnProbFromClass(bc$model, DT[["test"]])
#'
#' }
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

