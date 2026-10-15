#' refineB: Gradient-based fine-tuning for bHIVE antibodies with multiple loss
#' functions and optimizers
#'
#' After running \code{bHIVE} (or within \code{honeycombHIVE}), you have a set of
#' final antibody positions (A) in feature space. This function refines those
#' prototypes by iterating over the data points assigned to each antibody and
#' applying gradient-based updates using a user-chosen loss function and optimizer.
#'
#' The user must provide:
#' - A numeric matrix \code{A} of antibody/prototype positions (nAb x nFeatures).
#' - A numeric matrix/data frame \code{X} of data (nSamples x nFeatures).
#' - Optional \code{y} for classification. If \code{task="clustering"},
#'   \code{y} can be NULL or ignored.
#' - An integer or character vector \code{assignments} (length=nSamples) giving the
#'   antibody index (or label) to which each data point is assigned.
#'
#' ## Available Losses
#'
#' ### Classification (factor y)
#' - **"categorical_crossentropy"**: Pull prototypes toward data points if they share
#'   the antibody's dominant label; push away otherwise.
#' - **"binary_crossentropy"**: Similar to categorical CE, but we interpret factor y
#'   as binary (two classes). Pull for same label, push for different label.
#' - **"kullback_leibler"**: Very rough approach that treats "dominant label vs. others"
#'   as p and q distributions, pushing/pulling prototypes.
#' - **"cosine"**: Interpreted as trying to maximize cosine similarity for same-label points
#'   and minimize for different-label points.
#' - **"mae"**: Sign-based pull/push.
#'
#' ## Available Optimizers
#' - **"sgd"**: Basic stochastic gradient descent.
#' - **"momentum"**: SGD with momentum.
#' - **"adagrad"**: Adaptive gradient descent.
#' - **"adam"**: Adaptive moment estimation.
#' - **"rmsprop"**: Root Mean Square Propagation.
#'
#' @param A Numeric matrix (nAb x d) of antibody/prototype positions.
#' @param X Matrix or data frame (nSamples x d) of feature data.
#' @param y Optional. Factor (classification) or NULL (clustering).
#' @param assignments Vector (length = nSamples) giving the antibody each
#'   sample belongs to. Numeric values are row indices of \code{A}. Character
#'   or factor values are matched to \code{rownames(A)} when every value is a
#'   row name, and otherwise mapped to rows in order of first appearance.
#' @param task One of \code{c("clustering", "classification")}.
#' @param loss One of \code{c("categorical_crossentropy", "binary_crossentropy",
#'   "kullback_leibler", "cosine", "mae")}.
#' @param steps Integer. How many gradient steps to run.
#' @param lr Numeric. Learning rate for each update.
#' @param push_away Logical (for classification). Whether to push prototypes away
#'   from differently-labeled samples.
#' @param verbose Logical. If \code{TRUE}, reports progress with
#'   \code{message()} each step.
#' @param optimizer One of \code{c("sgd", "momentum", "adagrad", "adam", "rmsprop")}. Specifies the
#'   gradient-based optimization approach.
#' @param momentum_coef Numeric. Momentum coefficient (used if \code{optimizer=="momentum"}).
#' @param beta1 Numeric. Exponential decay rate for the first moment estimates (used in Adam).
#' @param beta2 Numeric. Exponential decay rate for the second moment estimates (used in Adam).
#' @param rmsprop_decay Numeric. Decay factor for the squared gradient moving average (used in RMSProp).
#' @param epsilon Numeric. A small constant for numerical stability.
#'
#' @return Updated matrix \code{A} of shape (nAb x d).
#'
#' @details
#' This function performs gradient-based updates on each antibody using the selected loss function.
#' Depending on the chosen \code{optimizer}, it may use plain SGD, momentum-based updates, Adagrad,
#' Adam, or RMSProp.
#'
#' @examples
#' data(iris)
#' X <- as.matrix(iris[, 1:4])
#' y <- iris$Species
#' res <- bHIVE(X, y, task = "classification", nAntibodies = 10,
#'              maxIter = 5, verbose = FALSE)
#' # antibody_index gives the antibody (row of res$antibodies) that labeled
#' # each sample
#' A_refined <- refineB(res$antibodies, X, y = y,
#'                      assignments = res$antibody_index,
#'                      task = "classification",
#'                      loss = "categorical_crossentropy", optimizer = "adam",
#'                      steps = 3, lr = 0.01, verbose = FALSE)
#' dim(A_refined)
#'
#' @export
refineB <- function(A, X, y = NULL, assignments,
                    task = c("clustering", "classification"),
                    loss = c("categorical_crossentropy", "binary_crossentropy",
                             "kullback_leibler", "cosine", "mae"),
                    steps = 5, lr = 0.01, push_away = TRUE, verbose = TRUE,
                    optimizer = c("sgd", "momentum", "adagrad", "adam",
                                  "rmsprop"),
                    momentum_coef = 0.9, beta1 = 0.9, beta2 = 0.999,
                    rmsprop_decay = 0.9, epsilon = 1e-8) {
  task      <- match.arg(task)
  loss      <- match.arg(loss)
  optimizer <- match.arg(optimizer)

  inp <- .refineB_check(A, X, assignments, verbose)
  A <- inp$A
  X <- inp$X
  assignments <- inp$assignments

  ab_label <- rep(NA_character_, nrow(A))
  if (task == "classification") {
    if (!is.factor(y)) y <- as.factor(y)
    ab_label <- .refineB_ab_labels(y, assignments, nrow(A))
  }
  step <- .refineB_optimizer(optimizer, nrow(A), ncol(A), lr, momentum_coef,
                             beta1, beta2, rmsprop_decay, epsilon)

  for (st in seq_len(steps)) {
    if (verbose) message(sprintf("Refinement step %d/%d", st, steps))
    for (i in sample.int(nrow(X))) {
      j <- assignments[i]
      # For classification, pull toward samples that share the antibody's
      # dominant label and (optionally) push away from the rest.
      same_label <- NA
      if (task == "classification" && !is.na(ab_label[j])) {
        same_label <- (y[i] == ab_label[j])
      }
      grad <- .update_prototype(A[j, ], as.numeric(X[i, ]), same_label,
                                task = task, loss = loss,
                                push_away = push_away)
      A[j, ] <- A[j, ] + step(j, grad)
    }
  }
  A
}


# ============================================================================
# Internal helpers
# ============================================================================

#' Validate refineB inputs and resolve assignments to row indices of A
#' @return list(A, X, assignments) with numeric matrices and integer indices.
#' @keywords internal
#' @noRd
.refineB_check <- function(A, X, assignments, verbose) {
  if (!is.matrix(X)) X <- as.matrix(X)
  if (!is.numeric(X)) stop("'X' must be numeric.")
  if (!is.matrix(A) || !is.numeric(A)) A <- as.matrix(A)
  nAb <- nrow(A)
  if (ncol(X) != ncol(A)) {
    stop("Number of columns in 'X' must match the number of columns in 'A'.")
  }
  if (length(assignments) != nrow(X)) {
    stop("Length of 'assignments' must match nrow(X).")
  }
  if (is.numeric(assignments)) {
    idx <- as.integer(assignments)
  } else if (!is.null(rownames(A)) &&
             all(as.character(assignments) %in% rownames(A))) {
    idx <- match(as.character(assignments), rownames(A))
  } else {
    idx <- match(assignments, unique(assignments))
  }
  bad <- is.na(idx) | idx < 1L | idx > nAb
  if (any(bad)) {
    stop(sprintf("'assignments' contains invalid values: %s. Valid range is 1..%d.",
                 paste(unique(assignments[bad]), collapse = ", "), nAb))
  }
  n_used <- length(unique(idx))
  if (verbose && n_used < nAb) {
    message(sprintf("%d of %d prototypes have no assigned samples and will not move.",
                    nAb - n_used, nAb))
  }
  list(A = A, X = X, assignments = idx)
}

#' Dominant class label of the samples assigned to each antibody
#' @keywords internal
#' @noRd
.refineB_ab_labels <- function(y, assignments, nAb) {
  vapply(seq_len(nAb), function(j) {
    pts <- which(assignments == j)
    if (length(pts) == 0) return(NA_character_)
    tbl <- table(y[pts])
    names(tbl)[which.max(tbl)]
  }, character(1))
}

#' Build a per-antibody update function for the chosen optimizer
#'
#' Returns step(j, grad), which updates the optimizer state for antibody j and
#' returns the displacement to add to A[j, ]. State lives in an environment
#' shared by the returned function.
#' @keywords internal
#' @noRd
.refineB_optimizer <- function(optimizer, nAb, d, lr, momentum_coef, beta1,
                               beta2, rmsprop_decay, epsilon) {
  st <- new.env(parent = emptyenv())
  st$s1 <- matrix(0, nrow = nAb, ncol = d)  # velocity / accumulator / 1st moment
  st$s2 <- matrix(0, nrow = nAb, ncol = d)  # 2nd moment (Adam)
  st$t  <- rep(0, nAb)
  switch(optimizer,
    sgd = function(j, g) lr * g,
    momentum = function(j, g) {
      st$s1[j, ] <- momentum_coef * st$s1[j, ] + g
      lr * st$s1[j, ]
    },
    adagrad = function(j, g) {
      st$s1[j, ] <- st$s1[j, ] + g^2
      lr * g / (sqrt(st$s1[j, ]) + epsilon)
    },
    adam = function(j, g) {
      st$s1[j, ] <- beta1 * st$s1[j, ] + (1 - beta1) * g
      st$s2[j, ] <- beta2 * st$s2[j, ] + (1 - beta2) * (g^2)
      st$t[j] <- st$t[j] + 1
      m_hat <- st$s1[j, ] / (1 - beta1^(st$t[j]))
      v_hat <- st$s2[j, ] / (1 - beta2^(st$t[j]))
      lr * m_hat / (sqrt(v_hat) + epsilon)
    },
    rmsprop = function(j, g) {
      st$s1[j, ] <- rmsprop_decay * st$s1[j, ] + (1 - rmsprop_decay) * (g^2)
      lr * g / (sqrt(st$s1[j, ]) + epsilon)
    })
}
