utils::globalVariables(c(
  "burnProb", "i.burned"))

#' Generate, Summarize, Predict Burn Classes from Covariates
#'
#' @param df A data.frame (or data.table), with covarites, including "burned" (a binary 0, 1
#' notburned = 0, burned = 1), e.g., timeSinceFire, biomassJackPine, etc. that will
#' be used to find fuel classes. This set of covariates must be available both during
#' fitting and for prediction. These must be quantitative.
#' @param numClasses A vector indicating how many classes should be attempted. The function
#' will return the number of classes that best classify the data into homogeneous groups.
#' @param AUC Logical. Should the Area Under the receiver operating Curve be returned?
#' @param plotAUC Logical. Should the plot of the AUC be made.
#'
#' @details
#' This was inspired by reading here:
#' \url{https://www.datanovia.com/en/blog/types-of-clustering-methods-overview-and-quick-start-r-code/}
#' and here:
#' \url{https://www.datanovia.com/en/lessons/model-based-clustering-essentials/}, with
#' citation here:
#' Scrucca L., Fop M., Murphy T. B. and Raftery A. E. (2016) mclust 5: clustering,
#' classification and density estimation using Gaussian finite mixture models,
#' The R Journal, 8/1, pp. 205-233.
#' \url{https://journal.r-project.org/archive/2016/RJ-2016-021/RJ-2016-021.pdf}
#'
#' @section The algorithm:
#' The basic solution is to take all covariates, including the binary "not burned", "burned"
#' (coded as 0 and 1, respectively), and do model-based clustering with the \code{mclust}
#' R package. We can choose a fixed number of burn classes, or a finite range (see
#' \code{numClasses} argument.
#' This will make \code{numClasses} "homogeneous" groups,
#' including whether they burned or not.
#' From this, we can identify groups by looking at the mean values of "burned" to see
#' what their burn tendency is as a "homogeneous" group.
#'
#' @section Categorical data:
#' For now, it is recommended to convert categorical data to dummy variables, 0 and 1.
#' E.g., For land cover, wetland class can be converted to a column "wetland" with
#' 1 for data points that are wetlands and 0 for non-wetland.
#'
#' @section How much data to include:
#' This has not been tested yet; however, I believe that having a relatively similar
#' number of "burned" and "unburned" pixels (within 3x either way), is probably a good idea.
#' In other words, if there are 100,000 burned data points, there should be between 30,000 and
#' 300,000 unburned data points. If there are already buffers around the burned polygons that
#' include unburned pixels, then these buffers can be used as part of the unburned
#' content.
#'
#' @return
#' A list with 2 elements, first the \code{model}, which comes from \code{mclust::Mclust},
#' and second the Area Under the Curve or AUC as an indicator of the overal goodness of fit.
#'
#' @rdname burnClass
#' @importFrom mclust Mclust mclustBIC
#' @importFrom pROC roc
#' @importFrom utils modifyList
#' @author Eliot McIntire
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
#'   DT[[i]][, c("jp", "bs", "ws", "age") := lapply(.SD, function(x) x/max(x) * 1000),
#'             .SDcols = c("jp", "bs", "ws", "age")]
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

  AUC <- NA
  if (isTRUE(AUC)) {
    summ <- burnClassSummary(mod)
    set(dftest, NULL, "burnClass", burnClassPredict(mod, dftest))
    dftest[summ[,c("burnClass", "burned")], burnProb := i.burned, on = "burnClass"]
    rocDF <- data.frame(predictions = dftest$burnProb, labels = dftest$burned)
    args <- list(rocDF$labels, rocDF$predictions, ci=TRUE, ci.alpha=0.9, stratified=FALSE, smoothed = TRUE)
    if (plotAUC) args <- modifyList(args, list(plot=TRUE, auc.polygon=TRUE, max.auc.polygon=TRUE, grid=TRUE,
                                               print.auc=TRUE, show.thres=TRUE))
    pROC_obj <- do.call(roc, args)
    AUC <- pROC_obj$auc
  }
  list(model = mod, AUC = AUC)
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
#' @importFrom stats predict
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

