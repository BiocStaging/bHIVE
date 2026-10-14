#' Multilayered honeycombHIVE classification model for caret
#'
#' A custom model definition that lets \code{caret::train()} tune, resample
#' and evaluate \code{\link{honeycombHIVE}} classifiers with gradient-based
#' refinement.
#'
#' @name honeycombHIVEmodel
#' @format A list with the components \pkg{caret} expects of a custom model.
#'
#' @details
#' Training runs \code{honeycombHIVE(refine = TRUE)} with the tuned
#' parameters. Further arguments given to \code{caret::train()} are passed on
#' to \code{honeycombHIVE()}. Progress output is off unless
#' \code{verbose = TRUE} is passed. Prediction uses the fitted
#' \code{\link{AINet}} model of the final layer. caret requires column names
#' on \code{x}. Class probabilities are not implemented (\code{prob = NULL}).
#'
#' @section Parameters:
#'   \itemize{
#'     \item \code{nAntibodies}: Number of initial antibodies per layer.
#'     \item \code{beta}: Clone multiplier controlling the number of clones per
#'           antibody.
#'     \item \code{epsilon}: Threshold for network suppression to remove
#'           redundant antibodies.
#'     \item \code{layers}: Number of hierarchical layers.
#'     \item \code{refineOptimizer}: Optimizer for gradient-based refinement
#'           ("sgd", "momentum", "adagrad", "adam" or "rmsprop").
#'     \item \code{refineSteps}: Number of gradient update steps in refinement.
#'     \item \code{refineLR}: Learning rate for refinement.
#'   }
#'
#' @return A list defining a custom \pkg{caret} model, for use as
#'   \code{caret::train(method = honeycombHIVEmodel)}.
#'
#' @examples
#' if (requireNamespace("caret", quietly = TRUE)) {
#'     X <- as.matrix(iris[, 1:4])
#'     y <- iris$Species
#'     set.seed(42)
#'     fit <- caret::train(
#'         x = X, y = y, method = honeycombHIVEmodel,
#'         trControl = caret::trainControl(method = "cv", number = 3),
#'         tuneGrid = data.frame(nAntibodies = 30, beta = 5, epsilon = 0.05,
#'                               layers = 2, refineOptimizer = "adam",
#'                               refineSteps = 3, refineLR = 0.01),
#'         maxIter = 10
#'     )
#'     fit
#'     table(Predicted = predict(fit, X), Actual = y)
#' }
#'
#' @seealso \code{\link{honeycombHIVE}}, \code{\link{bHIVEmodel}},
#'   \code{\link[caret]{train}}
#' @export
honeycombHIVEmodel <- list(
  label = "Multilayered Artificial Immune Network (honeycombHIVE)",
  library = "bHIVE",
  type = "Classification",
  parameters = data.frame(
    parameter = c("nAntibodies", "beta", "epsilon", "layers",
                  "refineOptimizer", "refineSteps", "refineLR"),
    class = c("numeric", "numeric", "numeric", "numeric",
              "character", "numeric", "numeric"),
    label = c("Number of Antibodies", "Clone Multiplier (beta)",
              "Epsilon Threshold", "Number of Layers", "Refinement Optimizer",
              "Refinement Steps", "Refinement Learning Rate")
  ),
  # Refinement settings are fixed in the default grid. Tune them with
  # tuneGrid.
  grid = function(x, y, len = NULL, search = "grid") {
    if (search == "grid") {
      g <- expand.grid(nAntibodies = round(seq(10, 50, length.out = len)),
                       beta        = seq(1, 10, length.out = len),
                       epsilon     = seq(0.01, 0.1, length.out = len),
                       layers      = unique(round(seq(1, 3, length.out = len))))
    } else {
      g <- data.frame(nAntibodies = sample(10:50, len, replace = TRUE),
                      beta        = runif(len, 1, 10),
                      epsilon     = runif(len, 0.01, 0.1),
                      layers      = sample(1:3, len, replace = TRUE))
    }
    g$refineOptimizer <- "adam"
    g$refineSteps <- 5
    g$refineLR <- 0.01
    g
  },
  fit = function(x, y, wts, param, lev, last, classProbs, ...) {
    if (!is.factor(y)) {
      stop("honeycombHIVEmodel supports classification only: y must be a factor.",
           call. = FALSE)
    }
    args <- c(list(X = as.matrix(x), y = y, task = "classification",
                   nAntibodies = param$nAntibodies, beta = param$beta,
                   epsilon = param$epsilon,
                   layers = as.integer(round(param$layers)), refine = TRUE,
                   refineSteps = param$refineSteps, refineLR = param$refineLR,
                   # tuning grids store strings as factors
                   refineOptimizer = as.character(param$refineOptimizer)),
              list(...))
    if (is.null(args$verbose)) args$verbose <- FALSE
    results <- do.call(honeycombHIVE, args)
    structure(list(task = "classification", results = results,
                   finalLayer = results[[length(results)]],
                   obsLevels = if (is.null(lev)) levels(y) else lev),
              class = "honeycombHIVEmodel")
  },
  predict = function(modelFit, newdata, submodels = NULL) {
    factor(modelFit$finalLayer$model$predict(as.matrix(newdata)),
           levels = modelFit$obsLevels)
  },
  prob = NULL
)
