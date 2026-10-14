#' honeycombHIVE: Multilayer AIS with optional gradient-based fine-tuning
#'
#' The \code{honeycombHIVE} function implements a multilayer artificial immune
#' system that iteratively refines a set of prototypes - referred to as
#' antibodies - to model the structure of the input data. In each layer, the
#' function first uses the \code{\link{bHIVE}} algorithm to generate or update
#' antibodies based on the current data representation and task (clustering or
#' classification). Optionally, it applies gradient-based
#' fine-tuning (via \code{\link{refineB}}) to these antibodies, allowing for
#' advanced refinement through various optimizers (e.g., SGD, Adam, RMSProp) and
#' customizable loss functions. The final output is a hierarchical set of layers
#' that encapsulate both the refined prototypes and the corresponding cluster
#' assignments or predictions for the original observations, making
#' \code{honeycombHIVE} a versatile tool for adaptive learning and pattern
#' recognition.
#'
#' @param X A numeric matrix or data frame of input features
#' (rows = observations, columns = features).
#' @param y Optional target factor vector for classification.
#' @param task Character, one of "clustering" or "classification".
#' @param layers Integer, how many layers (AIS iterations) to run.
#' @param nAntibodies Integer, how many antibodies (prototypes) to generate 
#' initially in each layer.
#' @param minAntibodies Integer, minimal number of antibodies to keep in each 
#' layer.
#' @param epsilon Numeric, threshold param for \code{bHIVE} suppression.
#' @param beta Numeric, selection pressure param for \code{bHIVE}.
#' @param maxIter Integer, maximum iterations for \code{bHIVE} each layer.
#' @param collapseMethod One of "centroid","medoid","median","mode".
#' @param minClusterSize Minimum cluster size. Clusters with fewer members are
#' handled per \code{smallClusterAction}. \code{NULL} (default) keeps every
#' cluster regardless of size.
#' @param smallClusterAction One of \code{"merge"} (default), \code{"keep"}, or
#' \code{"drop"}, controlling clusters below \code{minClusterSize}.
#' \code{"merge"} reassigns their members to the nearest surviving prototype so
#' no observation is lost (the rare cluster is absorbed, not deleted).
#' \code{"keep"} retains the small cluster as its own prototype (use this to
#' \emph{preserve} rare populations). \code{"drop"} is the legacy behavior that
#' discards the cluster and its members from subsequent layers. Ignored when
#' \code{minClusterSize} is \code{NULL}.
#' @param distance Distance metric for medoid calculation, e.g. "euclidean".
#' @param verbose Logical, if TRUE reports progress at each layer with
#' \code{message()}.
#' @param refine Logical, if TRUE apply gradient-based refinement via
#' \code{refineB()} to each layer's antibodies, in the scaled space the layer
#' was trained in. Refinement polishes the reported antibodies and the
#' final layer's model. The collapse into the next layer's prototypes uses
#' the bHIVE assignments.
#' @param refineLoss Character specifying the loss for \code{refineB()}
#' (e.g. "categorical_crossentropy", "mae").
#' @param refineSteps Integer, number of gradient steps in \code{refineB()}.
#' @param refineLR Numeric, learning rate for gradient updates.
#' @param refinePushAway Logical, if TRUE and classification, push prototypes
#' away from differently labeled points.
#' @param refineOptimizer Character, one of \code{"sgd", "momentum", "adagrad",
#' "adam", "rmsprop"} to be passed to \code{refineB()}.
#' @param refineMomentumCoef Numeric, momentum coefficient (if using momentum).
#' @param refineBeta1 Numeric, first moment decay rate (if using Adam).
#' @param refineBeta2 Numeric, second moment decay rate (if using Adam).
#' @param refineRmspropDecay Numeric, decay rate for the moving average of 
#' squared gradients (if using RMSProp).
#' @param refineEpsilon Numeric, a small constant for numerical stability 
#' (used in adaptive optimizers).
#' @param ... Additional arguments passed to \code{bHIVE}.
#'
#' @examples
#' # Clustering
#' data(iris)
#' X_iris <- iris[, 1:4]
#' resC <- honeycombHIVE(
#'   X = X_iris,
#'   task = "clustering",
#'   layers = 3,
#'   nAntibodies = 15,
#'   beta = 5,
#'   maxIter = 10
#' )
#'
#' @return A list of length \code{layers}. Each element (layer) is the
#'   \code{\link{bHIVE}} result for that layer plus:
#'   \itemize{
#'     \item \code{antibodies}: The antibodies learned in that layer.
#'     \item \code{assignments}: For each row of the layer's input (the
#'     original rows for layer 1, the previous layer's prototypes after that),
#'     the cluster ID or, for classification, the class label.
#'     \item \code{prototypes}: The collapsed prototypes that feed the next
#'     layer, one row per cluster. Row \code{k} summarizes the original
#'     observations with \code{membership == k}.
#'     \item \code{membership}: For each \strong{original} row in \code{X},
#'     the row of \code{prototypes} it belongs to (\code{NA} if its cluster
#'     was dropped).
#'     \item \code{predictions}: For classification, the predicted label for
#'     each original row: the majority label of its cluster.
#'     \item \code{model}: Final layer only. The fitted \code{\link{AINet}}
#'     model, for labeling new observations with \code{$predict()}.
#'     \item \code{task}: The specified task.
#'   }
#'   For classification, each layer collapses observations by the antibody
#'   that labeled them, so a layer keeps several prototypes per class.
#'
#' @importFrom stats median setNames
#' @export
honeycombHIVE <- function(X, y = NULL,
                          task = c("clustering", "classification"),
                          layers = 3, nAntibodies = 20, minAntibodies = 5,
                          epsilon = 0.05, beta = 5, maxIter = 10,
                          collapseMethod = c("centroid", "medoid", "median",
                                             "mode"),
                          minClusterSize = NULL,
                          smallClusterAction = c("merge", "keep", "drop"),
                          distance = "euclidean", verbose = TRUE,
                          refine = FALSE,
                          refineLoss = "categorical_crossentropy",
                          refineSteps = 5, refineLR = 0.01,
                          refinePushAway = TRUE, refineOptimizer = "sgd",
                          refineMomentumCoef = 0.9, refineBeta1 = 0.9,
                          refineBeta2 = 0.999, refineRmspropDecay = 0.9,
                          refineEpsilon = 1e-8, ...) {
  task <- match.arg(task)
  collapse <- list(method = match.arg(collapseMethod), distance = distance,
                   minClusterSize = minClusterSize,
                   action = match.arg(smallClusterAction))
  .validate_bHIVE_input(X, y)
  refine_args <- if (refine) {
    .honeycomb_refine_args(task, refineLoss, refineSteps, refineLR,
                           refinePushAway, refineOptimizer, refineMomentumCoef,
                           refineBeta1, refineBeta2, refineRmspropDecay,
                           refineEpsilon, verbose)
  }
  st <- .honeycomb_prepare(X, y)
  layer_args <- list(task = task, nAntibodies = nAntibodies,
                     minAntibodies = minAntibodies, epsilon = epsilon,
                     beta = beta, maxIter = maxIter, verbose = verbose,
                     dots = list(...))
  current <- st
  results <- vector("list", length = layers)
  for (layer_i in seq_len(layers)) {
    out <- .honeycomb_layer(current, st, layer_i, layers, layer_args,
                            refine_args, collapse)
    results[[layer_i]] <- out$result
    current <- out$next_input
  }
  results
}


# ============================================================================
# Internal helpers
# ============================================================================

#' Fit, optionally refine, and collapse one layer
#'
#' @param current list(X, y, rowIndices): this layer's input and its map back
#'   to the original rows.
#' @param st The prepared original data (.honeycomb_prepare()).
#' @return list(result = the layer's output, next_input = the next layer's
#'   input in the same shape as `current`).
#' @keywords internal
#' @noRd
.honeycomb_layer <- function(current, st, layer_i, layers, layer_args,
                             refine_args, collapse) {
  task    <- layer_args$task
  verbose <- layer_args$verbose
  if (verbose) {
    message(sprintf("=== honeycombHIVE: Layer %d / %d (task=%s) ===",
                    layer_i, layers, task))
  }
  res_layer <- .honeycomb_fit_layer(current$X, current$y, layer_args)
  if (!is.null(refine_args)) {
    res_layer <- .honeycomb_refine(res_layer, current$X, current$y, task,
                                   refine_args)
  }
  col <- .honeycomb_collapse(res_layer, current$X, current$rowIndices, st$y,
                             task, collapse, layer_i, nrow(st$X), verbose)
  res_layer$prototypes  <- col$prototypes
  res_layer$membership  <- col$membership
  res_layer$predictions <- col$predictions
  # Keep the fitted engine only on the final layer (for predict()); one R6
  # object per layer would bloat the returned hierarchy.
  if (layer_i < layers) res_layer$model <- NULL
  if (verbose) {
    message(sprintf("Layer %d completed. Next layer will use %d prototypes.",
                    layer_i, nrow(col$prototypes)))
  }
  list(result = res_layer,
       next_input = list(X = as.data.frame(col$prototypes),
                         y = .honeycomb_next_y(col$rowIndices, st$y, task),
                         rowIndices = col$rowIndices))
}

#' Validate the refine settings and bundle them for refineB()
#' @keywords internal
#' @noRd
.honeycomb_refine_args <- function(task, loss, steps, lr, push_away,
                                   optimizer, momentum_coef, beta1, beta2,
                                   rmsprop_decay, epsilon, verbose) {
  valid <- c("categorical_crossentropy", "binary_crossentropy",
             "kullback_leibler", "cosine", "mae")
  if (!loss %in% valid) {
    stop(sprintf("Invalid refineLoss '%s' for task '%s'. Supported losses: %s",
                 loss, task, paste(valid, collapse = ", ")))
  }
  list(loss = loss, steps = steps, lr = lr,
       push_away = task == "classification" && push_away,
       verbose = verbose, optimizer = as.character(optimizer),
       momentum_coef = momentum_coef, beta1 = beta1, beta2 = beta2,
       rmsprop_decay = rmsprop_decay, epsilon = epsilon)
}

#' Coerce X to a data frame with row names and map rows to original indices
#' @keywords internal
#' @noRd
.honeycomb_prepare <- function(X, y) {
  X <- as.data.frame(X)
  n <- nrow(X)
  if (is.null(rownames(X))) rownames(X) <- paste0("row_", seq_len(n))
  if (!is.null(y) && is.null(names(y))) names(y) <- rownames(X)
  list(X = X, y = y,
       rowIndices = setNames(as.list(seq_len(n)), rownames(X)))
}

#' Run bHIVE on one layer's input
#' @keywords internal
#' @noRd
.honeycomb_fit_layer <- function(current_X, current_y, a) {
  n_current <- nrow(current_X)
  args <- list(X = current_X, y = current_y, task = a$task,
               nAntibodies = min(a$nAntibodies, max(a$minAntibodies, n_current)),
               epsilon = a$epsilon, beta = a$beta, maxIter = a$maxIter,
               verbose = a$verbose)
  args <- c(args, a$dots[setdiff(names(a$dots), names(args))])
  do.call(.bHIVE_matrix, args)
}

#' Refine a layer's antibodies with refineB() in the model's scaled space
#'
#' Each sample pulls (or pushes) the antibody that actually labeled it:
#' the cluster ID for clustering, antibody_index for classification.
#' Clustering uses a single dummy class so every sample pulls.
#' @keywords internal
#' @noRd
.honeycomb_refine <- function(res_layer, current_X, current_y, task,
                              refine_args) {
  model   <- res_layer$model
  scaling <- model$config$scaling
  Xs <- .bhive_apply_scaling(as.matrix(current_X), scaling)
  if (task == "classification") {
    idx <- res_layer$antibody_index
    y_ref <- current_y
  } else {
    idx <- as.integer(res_layer$assignments)
    y_ref <- factor(rep("c", nrow(Xs)))
  }
  A <- do.call(refineB, c(list(A = res_layer$antibodies, X = Xs, y = y_ref,
                               assignments = idx, task = "classification"),
                          refine_args))
  res_layer$antibodies <- A
  res_layer$antibodies_unscaled <- .bhive_invert_scaling(A, scaling)
  model$result$antibodies <- A
  model$result$antibodies_unscaled <- res_layer$antibodies_unscaled
  res_layer
}

#' Collapse a layer's clusters into prototypes for the next layer
#'
#' Groups the layer's input rows by cluster (clustering) or by the antibody
#' that labeled them (classification), summarizes each group into one
#' prototype, handles small clusters, and maps the groups back to the
#' original observations.
#' @return list(prototypes, rowIndices, membership, predictions)
#' @keywords internal
#' @noRd
.honeycomb_collapse <- function(res_layer, current_X, rowIndices, y, task,
                                collapse, layer_i, n_original, verbose) {
  ids <- if (task == "classification") {
    res_layer$antibody_index
  } else {
    res_layer$assignments
  }
  groups  <- split(rownames(current_X), f = factor(ids, levels = unique(ids)))
  members <- lapply(groups, function(rn) unlist(rowIndices[rn], use.names = FALSE))
  protos  <- do.call(rbind, lapply(groups, function(rn) {
    .honeycomb_collapse_proto(current_X[rn, , drop = FALSE], collapse$method,
                              collapse$distance)
  }))
  colnames(protos) <- colnames(current_X)
  is_small <- if (is.null(collapse$minClusterSize)) {
    rep(FALSE, length(members))
  } else {
    lengths(members) < collapse$minClusterSize
  }
  kept <- .honeycomb_small_clusters(protos, members, is_small,
                                    collapse$action, verbose)
  if (nrow(kept$protos) == 0) {
    stop("No valid prototypes generated. Check parameters & input data.")
  }
  row_names <- paste0("Layer", layer_i, "_Cluster", seq_len(nrow(kept$protos)))
  rownames(kept$protos) <- row_names
  names(kept$members)   <- row_names
  .honeycomb_finalize_layer(kept, y, task, n_original)
}

#' Summarize one cluster's rows into a single prototype vector
#' @keywords internal
#' @noRd
.honeycomb_collapse_proto <- function(sub_data, method, distance) {
  switch(
    method,
    "centroid" = colMeans(sub_data, na.rm = TRUE),
    "medoid" = {
      dmat <- as.matrix(dist(sub_data, method = distance))
      as.numeric(sub_data[which.min(rowSums(dmat)), , drop = FALSE])
    },
    "median" = apply(sub_data, 2, median, na.rm = TRUE),
    "mode" = apply(sub_data, 2, function(vec) {
      tb <- table(vec)
      as.numeric(names(tb)[which.max(tb)])
    })
  )
}

#' Merge, keep or drop clusters below minClusterSize
#'
#' "merge" hands each small cluster's members to the nearest surviving
#' prototype so no observation is lost. "keep" retains small clusters (rare
#' populations). "drop" discards them and their members.
#' @keywords internal
#' @noRd
.honeycomb_small_clusters <- function(protos, members, is_small, action,
                                      verbose) {
  if (!any(is_small)) return(list(protos = protos, members = members))
  big   <- which(!is_small)
  small <- which(is_small)
  keep  <- rep(TRUE, length(members))
  if (action == "merge" && length(big) > 0) {
    big_protos <- protos[big, , drop = FALSE]
    for (si in small) {
      tgt <- big[which.min(colSums((t(big_protos) - protos[si, ])^2))]
      members[[tgt]] <- c(members[[tgt]], members[[si]])
    }
    keep <- !is_small
    if (verbose) {
      message(sprintf("Merging %d small cluster(s) into nearest prototype.",
                      length(small)))
    }
  } else if (action == "drop") {
    keep <- !is_small
    if (verbose) {
      message(sprintf("Dropping %d small cluster(s) and their members.",
                      length(small)))
    }
  } else if (verbose && action == "keep") {
    message(sprintf("Keeping %d small cluster(s) as rare prototypes.",
                    length(small)))
  }
  list(protos = protos[keep, , drop = FALSE], members = members[keep])
}

#' Map the final clusters back to original rows
#'
#' membership[i] is the prototype row of original observation i. For
#' classification, every original row is predicted as the majority original
#' label of its cluster.
#' @keywords internal
#' @noRd
.honeycomb_finalize_layer <- function(kept, y, task, n_original) {
  membership  <- rep(NA_integer_, n_original)
  predictions <- if (task == "classification") rep(NA_character_, n_original)
  for (k in seq_along(kept$members)) {
    idx <- kept$members[[k]]
    membership[idx] <- k
    if (task == "classification") predictions[idx] <- .majority_label(y[idx])
  }
  list(prototypes = kept$protos, rowIndices = kept$members,
       membership = membership, predictions = predictions)
}

#' Class labels for the next layer's prototypes (majority of original labels)
#' @keywords internal
#' @noRd
.honeycomb_next_y <- function(rowIndices, y, task) {
  if (task != "classification" || is.null(y)) return(NULL)
  labs <- vapply(rowIndices, function(ix) .majority_label(y[ix]), character(1))
  stats::setNames(factor(labs, levels = levels(y)), names(rowIndices))
}

#' Most frequent label (first in level order on ties)
#' @keywords internal
#' @noRd
.majority_label <- function(y) {
  tb <- table(y)
  names(tb)[which.max(tb)]
}
