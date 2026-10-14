#' Null-coalesce operator
#'
#' Returns \code{x} if it is not \code{NULL}, otherwise returns \code{y}.
#' Used internally by bHIVE R6 classes for parameter defaults.
#'
#' @param x Value to test.
#' @param y Default value if \code{x} is NULL.
#' @return \code{x} if not NULL, otherwise \code{y}.
#' @name null-coalesce
#' @keywords internal
`%||%` <- function(x, y) if (is.null(x)) y else x

.validate_bHIVE_input <- function(X, 
                                  y = NULL) {
  if (methods::is(X, "SummarizedExperiment")) {
    stop("X is a SummarizedExperiment. Pass t(assay(X)) or reducedDim(X) ",
         "here, or use bHIVE() or bHIVEParam(), which accept it directly.",
         call. = FALSE)
  }
  if (!is.matrix(X) && !is.data.frame(X)) {
    stop("Input X must be a matrix or data frame.")
  }
  if (!is.null(y)) {
    if (!is.factor(y)) {
      stop("y must be a factor (classification).")
    }
    if (nrow(X) != length(y)) {
      stop("X and y must have the same number of rows.")
    }
  }
  if (anyNA(X)) {
    stop("Input X contains missing values. Please handle missing values before running bHIVE.")
  }
  invisible(TRUE)
}

#' Gradient direction for one prototype and one sample
#'
#' Clustering always pulls the prototype toward the sample. Classification
#' pulls when the sample shares the prototype's dominant label and, with
#' push_away, pushes it away otherwise. The "mae" loss uses only the sign of
#' that direction.
#' @return Numeric vector of length(ab_vec).
#' @keywords internal
#' @noRd
.update_prototype <- function(ab_vec,
                              x_i,
                              same_label = NA,
                              task = c("clustering","classification"),
                              loss = c("categorical_crossentropy", "binary_crossentropy",
                                       "kullback_leibler", "cosine", "mae"),
                              push_away = TRUE) {
  task <- match.arg(task)
  loss <- match.arg(loss)
  pull <- x_i - ab_vec
  if (task == "clustering") return(pull)
  # No label to compare against, or a mismatch without push_away: no move.
  if (is.na(same_label)) return(numeric(length(ab_vec)))
  if (isTRUE(same_label)) {
    dir <- pull
  } else if (push_away) {
    dir <- -pull
  } else {
    return(numeric(length(ab_vec)))
  }
  if (loss == "mae") sign(dir) else dir
}
