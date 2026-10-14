#' bHIVE classification model for caret
#'
#' A custom model definition that lets \code{caret::train()} tune, resample
#' and evaluate \code{\link{bHIVE}} classifiers.
#'
#' @name bHIVEmodel
#' @format A list with the components \pkg{caret} expects of a custom model.
#'
#' @details
#' Pass \code{bHIVEmodel} as the \code{method} of \code{caret::train()}.
#' Training fits \code{\link{bHIVE}} with the tuned \code{nAntibodies},
#' \code{beta} and \code{epsilon}. Further arguments given to
#' \code{caret::train()} (for example \code{maxIter}, \code{scale} or
#' immunology modules) are passed on to \code{bHIVE()}. Progress output is
#' off unless \code{verbose = TRUE} is passed. Prediction labels each new
#' observation with the class of its best-matching antibody, using the fitted
#' \code{\link{AINet}} model, so new data gets the same scaling and kernel as
#' the training data.
#'
#' caret requires column names on \code{x}.
#'
#' @section Components:
#' \describe{
#'   \item{\code{label}}{"Artificial Immune Network (bHIVE)".}
#'   \item{\code{library}}{"bHIVE".}
#'   \item{\code{type}}{"Classification".}
#'   \item{\code{parameters}}{The tunable parameters: \code{nAntibodies}
#'   (initial antibody count), \code{beta} (clone multiplier) and
#'   \code{epsilon} (suppression threshold).}
#'   \item{\code{grid}}{\code{function(x, y, len, search)}. A default tuning
#'   grid, regular (\code{search = "grid"}) or random.}
#'   \item{\code{fit}}{\code{function(x, y, wts, param, lev, last, classProbs,
#'   ...)}. Fits \code{bHIVE()} and returns its result.}
#'   \item{\code{predict}}{\code{function(modelFit, newdata, submodels)}.
#'   Returns a factor of predicted classes.}
#'   \item{\code{prob}}{\code{NULL}. Class probabilities are not implemented,
#'   so \code{caret::train()} warns and continues if \code{classProbs = TRUE}.}
#' }
#'
#' @return A list defining a custom \pkg{caret} model, for use as
#'   \code{caret::train(method = bHIVEmodel)}.
#'
#' @examples
#' bHIVEmodel$label
#' bHIVEmodel$parameters
#'
#' if (requireNamespace("caret", quietly = TRUE)) {
#'     X <- as.matrix(iris[, 1:4])
#'     y <- iris$Species
#'     set.seed(42)
#'     fit <- caret::train(
#'         x = X, y = y, method = bHIVEmodel,
#'         trControl = caret::trainControl(method = "cv", number = 3),
#'         tuneGrid = data.frame(nAntibodies = 40, beta = 5, epsilon = 0.05),
#'         maxIter = 10
#'     )
#'     fit
#'     table(Predicted = predict(fit, X), Actual = y)
#' }
#'
#' @seealso \code{\link{bHIVE}}, \code{\link{honeycombHIVEmodel}},
#'   \code{\link[caret]{train}}
#'
#' @importFrom stats runif
#' @export
bHIVEmodel <- list(
  label = "Artificial Immune Network (bHIVE)",
  library = "bHIVE",
  type = "Classification",
  parameters = data.frame(
    parameter = c("nAntibodies", "beta", "epsilon"),
    class = c("numeric", "numeric", "numeric"),
    label = c("Number of Antibodies", "Clone Multiplier", "Epsilon Threshold")
  ),
  grid = function(x, y, len = NULL, search = "grid") {
    if (search == "grid") {
      expand.grid(nAntibodies = round(seq(10, 50, length.out = len)),
                  beta        = seq(1, 10, length.out = len),
                  epsilon     = seq(0.01, 0.1, length.out = len))
    } else {
      data.frame(nAntibodies = sample(10:50, len, replace = TRUE),
                 beta        = runif(len, 1, 10),
                 epsilon     = runif(len, 0.01, 0.1))
    }
  },
  fit = function(x, y, wts, param, lev, last, classProbs, ...) {
    if (!is.factor(y)) {
      stop("bHIVEmodel supports classification only: y must be a factor.",
           call. = FALSE)
    }
    args <- c(list(X = as.matrix(x), y = y, task = "classification",
                   nAntibodies = param$nAntibodies, beta = param$beta,
                   epsilon = param$epsilon),
              list(...))
    if (is.null(args$verbose)) args$verbose <- FALSE
    res <- do.call(.bHIVE_matrix, args)
    # caret sets obsLevels itself; keep a copy for direct fit/predict calls.
    res$obsLevels <- if (is.null(lev)) levels(y) else lev
    res
  },
  predict = function(modelFit, newdata, submodels = NULL) {
    factor(modelFit$model$predict(as.matrix(newdata)),
           levels = modelFit$obsLevels)
  },
  prob = NULL
)
