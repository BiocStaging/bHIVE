#' @title ActivationGate
#' @description Two-signal activation gate implementing the immunological
#' principle that immune cell activation requires both antigen-specific
#' recognition (Signal 1) AND costimulatory context (Signal 2).
#'
#' @details
#' Prevents spurious activation on isolated outliers. An antibody is only
#' allowed to clone if both signals exceed their thresholds. This is
#' biologically-principled regularization.
#'
#' Signal 2 options:
#' \itemize{
#'   \item \code{"density"}: Local data density around the antibody
#'   \item \code{"danger"}: User-provided danger signal vector
#'   \item \code{"entropy"}: Local label entropy (classification only)
#' }
#'
#' @examples
#' # Two-signal activation gate
#' data(iris)
#' X <- as.matrix(iris[, 1:4])
#' A <- X[sample(150, 10), ]
#' rep <- ImmuneRepertoire$new(A)
#' gate <- ActivationGate$new(signal2_type = "density", threshold2 = 0.3)
#' aff <- rep$affinity_matrix(X, "gaussian")
#' activated <- gate$evaluate(aff, X, A)
#' sum(activated)  # number of activated interactions
#'
#' @param signal2_type Character. "density", "danger", or "entropy".
#' @param threshold1 Numeric. Minimum affinity threshold.
#' @param threshold2 Numeric. Minimum Signal 2 threshold.
#' @param danger_signals Numeric vector. Per-data-point danger scores.
#'
#' @return An \code{\link[R6]{R6Class}} generator. \code{ActivationGate$new()}
#'   returns an \code{ActivationGate} object.
#'
#' @importFrom R6 R6Class
#' @export
ActivationGate <- R6::R6Class(
  "ActivationGate",
  public = list(

    #' @field signal2_type Type of costimulatory signal.
    signal2_type = NULL,

    #' @field threshold1 Minimum affinity for Signal 1 (antigen recognition).
    threshold1 = NULL,

    #' @field threshold2 Minimum costimulatory signal for Signal 2.
    threshold2 = NULL,

    #' @field danger_signals User-provided danger signal vector (for "danger" type).
    danger_signals = NULL,

    #' @description Create a new ActivationGate.
    #' @param signal2_type Character. "density", "danger", or "entropy".
    #' @param threshold1 Numeric. Minimum affinity threshold.
    #' @param threshold2 Numeric. Minimum Signal 2 threshold.
    #' @param danger_signals Numeric vector. Per-data-point danger scores.
    initialize = function(signal2_type = "density",
                          threshold1 = 0.1,
                          threshold2 = 0.3,
                          danger_signals = NULL) {
      self$signal2_type   <- match.arg(signal2_type, c("density", "danger", "entropy"))
      self$threshold1     <- threshold1
      self$threshold2     <- threshold2
      self$danger_signals <- danger_signals
    },

    #' @description Evaluate which antibody-data interactions pass the two-signal gate.
    #' @param affinity_matrix Numeric matrix (n x m) of affinities.
    #' @param X Numeric matrix of data (n x d).
    #' @param A Numeric matrix of antibodies (m x d).
    #' @param y Target vector or NULL.
    #' @param task Character. Task type.
    #' @return Logical matrix (n x m) where TRUE means the interaction is activated.
    evaluate = function(affinity_matrix, X, A, y = NULL, task = "clustering") {
      n <- nrow(affinity_matrix)
      m <- ncol(affinity_matrix)

      # Signal 1: affinity exceeds threshold
      signal1 <- affinity_matrix > self$threshold1

      # Signal 2: costimulatory context
      signal2 <- switch(self$signal2_type,
        density = .gate_density_signal(affinity_matrix, self$threshold2),
        danger  = {
          if (is.null(self$danger_signals)) {
            stop("danger_signals must be provided for signal2_type='danger'")
          }
          matrix(self$danger_signals > self$threshold2, nrow = n, ncol = m)
        },
        entropy = if (is.null(y) || task != "classification") {
          # Fall back to density for non-classification tasks
          .gate_density_signal(affinity_matrix, self$threshold2)
        } else {
          .gate_entropy_signal(affinity_matrix, y, self$threshold2)
        })

      # Both signals must pass
      signal1 & signal2
    },

    #' @description Format a short summary, one element per line. Used by
    #'   \code{print()}.
    #' @param ... Not used.
    #' @return A character vector.
    format = function(...) {
      c(sprintf("<ActivationGate> signal2='%s'", self$signal2_type),
        sprintf("  Threshold1 (affinity): %.3f", self$threshold1),
        sprintf("  Threshold2 (context):  %.3f", self$threshold2))
    }
  )
)


#' Signal 2 by data density
#'
#' A data point passes (for every antibody) when its summed affinity to all
#' antibodies is above the q-th quantile.
#' @keywords internal
#' @noRd
.gate_density_signal <- function(affinity_matrix, q) {
  density_per_point <- rowSums(affinity_matrix)
  dense <- density_per_point > quantile(density_per_point, q)
  matrix(dense, nrow = nrow(affinity_matrix), ncol = ncol(affinity_matrix))
}

#' Signal 2 by label purity
#'
#' Each data point is assigned to its best-binding antibody. The pair passes
#' when the label entropy of that antibody's points is below the threshold
#' (lower entropy means a purer, more trustworthy context).
#' @keywords internal
#' @noRd
.gate_entropy_signal <- function(affinity_matrix, y, threshold) {
  signal2 <- matrix(FALSE, nrow = nrow(affinity_matrix),
                    ncol = ncol(affinity_matrix))
  assignments <- apply(affinity_matrix, 1, which.max)
  for (j in seq_len(ncol(affinity_matrix))) {
    assigned <- which(assignments == j)
    if (length(assigned) > 1) {
      freqs <- table(y[assigned]) / length(assigned)
      entropy <- -sum(freqs * log(freqs + 1e-15))
      if (entropy < threshold) signal2[assigned, j] <- TRUE
    }
  }
  signal2
}
