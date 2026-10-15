#' Tune Hyperparameters for bHIVE (Swarm/Grid Search)
#'
#' Performs hyperparameter tuning for the bHIVE algorithm over a grid of
#' hyperparameter values or an externally provided data frame of parameter
#' combinations. Evaluates each combination using different metrics:
#'
#' - **Classification**: "accuracy", "balanced_accuracy", "f1", "kappa"
#' - **Clustering**: "silhouette", "davies_bouldin", or "calinski_harabasz"
#'
#' **Note**: Some metrics require additional packages or assumptions
#' (e.g., multi-class classification for "f1" is calculated as a macro-average).
#'
#' @param X A numeric matrix or data frame of input features (rows =
#' observations, columns = features).
#' @param y Optional. A factor target vector for classification.
#'   If \code{NULL}, clustering is performed.
#' @param task Character. One of \code{"clustering"} or \code{"classification"}.
#' @param grid A data frame specifying the hyperparameter combinations, one
#' row per combination. Any column named after a \code{\link{bHIVE}} or
#' \code{\link{AINet}} argument is tuned (for example \code{nAntibodies},
#' \code{beta}, \code{epsilon}, \code{targetK}, \code{scale},
#' \code{affinityFunc} or \code{coverageBoost}). Other columns are ignored.
#' Factor columns, as produced by \code{expand.grid()}, are converted to
#' character.
#' @param metric Character. Name of the evaluation metric. Options:
#'   \itemize{
#'     \item \strong{Classification}: "accuracy", "balanced_accuracy", "f1", "kappa"
#'     \item \strong{Clustering}: "silhouette", "davies_bouldin", "calinski_harabasz"
#'   }
#' @param maxIter Integer. Maximum iterations for each \code{bHIVE} run
#' (default 50). A \code{maxIter} column in \code{grid} takes precedence.
#' @param BPPARAM Character. A BiocParallel::bpparam() object that can be used 
#' for parallelization. The function supports \code{SerialParam}, \code{MulticoreParam}, 
#' \code{BatchtoolsParam}, and \code{SnowParam}.
#' @param verbose Logical. If \code{TRUE}, reports progress with
#' \code{message()}.
#'
#' @return A list:
#'   \itemize{
#'     \item \code{best_params}: A list (row) of the best hyperparameters.
#'     \item \code{results}: A data frame with the full grid search results, 
#'       including the \code{metric_value} for each combination.
#'   }
#'
#' @examples
#' data(iris)
#' X <- as.matrix(iris[, 1:4])
#' y <- iris$Species  # classification
#'
#' # Define hyperparameter grid
#' grid <- expand.grid(
#'   nAntibodies = c(10, 20),
#'   beta        = c(3, 5),
#'   epsilon     = c(0.01, 0.05)
#' )
#'
#' # Perform hyperparameter tuning for classification
#' tuning_results <- swarmbHIVE(X = X, 
#'                              y = y, 
#'                              task = "classification", 
#'                              grid = grid, 
#'                              metric = "balanced_accuracy",
#'                              maxIter = 10)
#'
#' # For clustering with silhouette
#' set.seed(42)
#' X_clust <- matrix(rnorm(100 * 5), ncol = 5)
#' grid_clust <- expand.grid(nAntibodies = c(5, 10), 
#'                           beta = c(3, 5), 
#'                           epsilon = c(0.01, 0.05))
#' res_clust <- swarmbHIVE(X_clust, 
#'                         task = "clustering", 
#'                         grid = grid_clust, 
#'                         metric = "silhouette")
#' res_clust$best_params
#' 
#'
#' @importFrom stats dist
#' @importFrom cluster silhouette
#' @importFrom BiocParallel bplapply SerialParam
#' @importFrom clusterCrit intCriteria
#' @export
swarmbHIVE <- function(X,
                       y = NULL,
                       task = c("clustering","classification"),
                       grid,
                       metric = NULL,
                       maxIter = 50,
                       BPPARAM = SerialParam(),
                       verbose = TRUE) {
  task <- match.arg(task)
  .validate_bHIVE_input(X, y)
  metric <- .swarm_check_metric(task, metric)

  # Silhouette needs all pairwise distances. Compute them once, not per combo.
  dist_mat <- NULL
  if (metric == "silhouette") {
    if (verbose) message("Precomputing distance matrix for the silhouette metric.")
    dist_mat <- dist(X)
  }

  n_combos <- nrow(grid)
  if (verbose) {
    message(sprintf(
      "Starting swarmbHIVE with %d parameter combinations (task=%s, metric=%s).",
      n_combos, task, metric))
  }
  results_list <- bplapply(seq_len(n_combos), function(i) {
    row <- grid[i, , drop = FALSE]
    if (verbose) {
      message(sprintf("Evaluating combo %d/%d: %s", i, n_combos,
                      .swarm_format_row(row)))
    }
    .swarm_eval_combo(row, X, y, task, maxIter, metric, dist_mat)
  }, BPPARAM = BPPARAM)

  results_df <- do.call(rbind, results_list)
  rownames(results_df) <- NULL
  best_params <- results_df[.swarm_best_index(results_df$metric_value, metric),
                            , drop = FALSE]
  if (verbose) message("Best parameters: ", .swarm_format_row(best_params))
  list(best_params = best_params, results = results_df)
}


# ============================================================================
# Internal helpers
# ============================================================================

# Valid metrics per task. The first is the default.
.SWARM_METRICS <- list(
  classification = c("accuracy", "balanced_accuracy", "f1", "kappa"),
  clustering     = c("silhouette", "davies_bouldin", "calinski_harabasz"))

#' Resolve and validate the tuning metric for a task
#' @keywords internal
#' @noRd
.swarm_check_metric <- function(task, metric) {
  valid <- .SWARM_METRICS[[task]]
  if (is.null(metric)) return(valid[1])
  if (!metric %in% valid) {
    stop(sprintf("Invalid metric '%s' for %s. Valid: %s",
                 metric, task, paste(valid, collapse = ", ")))
  }
  metric
}

#' Fit bHIVE for one grid row and score it
#'
#' Every grid column named after a bHIVE/AINet argument is forwarded, so the
#' swarm can tune any of them. Factor columns (expand.grid's default) become
#' character so match.arg() accepts them.
#' @param row One-row data frame.
#' @return The row with a metric_value column appended.
#' @keywords internal
#' @noRd
.swarm_eval_combo <- function(row, X, y, task, maxIter, metric, dist_mat) {
  args <- lapply(as.list(row), function(v) {
    if (is.factor(v)) as.character(v) else v
  })
  args <- args[intersect(names(args), .bhive_tunable_args())]
  if (is.null(args$maxIter)) args$maxIter <- maxIter
  args <- c(args, list(X = X, y = y, task = task, verbose = FALSE))
  fit <- do.call(.bHIVE_matrix, args)
  cbind(row, metric_value = .swarm_metric(fit$assignments, X, y, task,
                                          metric, dist_mat))
}

#' Score a set of labels with a clustering or classification metric
#' @keywords internal
#' @noRd
.swarm_metric <- function(labels, X, y, task, metric, dist_mat) {
  if (task == "clustering") {
    .swarm_cluster_metric(labels, X, metric, dist_mat)
  } else {
    .swarm_class_metric(labels, y, metric)
  }
}

#' Internal clustering metrics (silhouette, Davies-Bouldin, Calinski-Harabasz)
#' @keywords internal
#' @noRd
.swarm_cluster_metric <- function(labels, X, metric, dist_mat) {
  part <- as.integer(factor(labels))
  # Every internal index is undefined for a single cluster.
  if (length(unique(part)) < 2L) return(NA_real_)
  if (metric == "silhouette") {
    sil <- silhouette(part, dist_mat)
    return(mean(sil[, "sil_width"], na.rm = TRUE))
  }
  intCriteria(traj = as.matrix(X), part = part, crit = metric)[[metric]]
}

#' Classification metrics from predicted and true labels
#'
#' The confusion matrix is built over the union of true and predicted labels
#' so it is always square, even when the model never predicts some class.
#' @keywords internal
#' @noRd
.swarm_class_metric <- function(pred, y, metric) {
  lvls <- sort(unique(c(as.character(y), as.character(pred))))
  tb <- table(factor(as.character(y), levels = lvls),
              factor(as.character(pred), levels = lvls))  # rows = truth
  switch(metric,
         accuracy          = mean(as.character(pred) == as.character(y)),
         balanced_accuracy = .balanced_accuracy(tb),
         f1                = .macro_f1(tb),
         kappa             = .cohen_kappa(tb))
}

#' Mean per-class recall (rows = truth, columns = prediction)
#' @keywords internal
#' @noRd
.balanced_accuracy <- function(tb) {
  mean(diag(prop.table(tb, margin = 1)), na.rm = TRUE)
}

#' Macro-averaged F1. A class with undefined precision or recall scores 0.
#' @keywords internal
#' @noRd
.macro_f1 <- function(tb) {
  tp   <- diag(tb)
  prec <- tp / colSums(tb)
  rec  <- tp / rowSums(tb)
  f1 <- ifelse(is.na(prec) | is.na(rec) | (prec + rec) == 0, 0,
               2 * prec * rec / (prec + rec))
  mean(f1)
}

#' Cohen's kappa: (p_o - p_e) / (1 - p_e)
#' @keywords internal
#' @noRd
.cohen_kappa <- function(tb) {
  n  <- sum(tb)
  p0 <- sum(diag(tb)) / n
  pE <- sum(rowSums(tb) * colSums(tb)) / (n * n)
  (p0 - pE) / (1 - pE)
}

#' Row index of the best score. Davies-Bouldin is smaller-is-better.
#' @keywords internal
#' @noRd
.swarm_best_index <- function(values, metric) {
  if (metric == "davies_bouldin") which.min(values) else which.max(values)
}

#' Format a one-row data frame as "name=value, ..." for progress messages
#' @keywords internal
#' @noRd
.swarm_format_row <- function(row) {
  vals <- vapply(row, function(v) {
    format(if (is.factor(v)) as.character(v) else v, digits = 4)
  }, character(1))
  paste(names(row), vals, sep = "=", collapse = ", ")
}
