#' bHIVE clustering through bluster
#'
#' Run \code{\link{bHIVE}} clustering through \code{\link[bluster]{clusterRows}},
#' the \pkg{bluster} interface that \code{scran::clusterCells()} and
#' \code{bluster::clusterSweep()} build on. This makes bHIVE a drop-in
#' alternative to the k-means, SOM and graph-based clustering that
#' \pkg{bluster} already wraps.
#'
#' @param nAntibodies,beta,epsilon,maxIter,k,affinityFunc,distFunc,initMethod,scale,targetK,epsilonQuantile
#'   Scalar arguments passed to \code{\link{bHIVE}}. \code{NULL} (the default)
#'   uses the \code{bHIVE()} default. These can be varied with
#'   \code{bluster::clusterSweep()}.
#' @param affinityParams Named list passed to \code{bHIVE()}, for example
#'   \code{list(alpha = "auto")}. An empty list uses the \code{bHIVE()}
#'   default.
#' @param ... Further named arguments passed to \code{bHIVE()} and
#'   \code{\link{AINet}}, such as \code{coverageBoost}, \code{consolidate},
#'   \code{verbose} or immunology modules (\code{shm}, \code{idiotypic}, ...).
#'   Note that R6 modules carry state, so a \code{MemoryPool} passed here is
#'   shared across every run that uses this parameter object.
#' @param x A numeric matrix-like object with observations in rows, or a
#'   \code{bHIVEParam} for \code{show()}.
#' @param object A \code{bHIVEParam} object.
#' @param BLUSPARAM A \code{bHIVEParam} object.
#' @param full Logical. Return the full \code{bHIVE()} result as well as the
#'   clusters.
#'
#' @details
#' bHIVE learns a set of prototypes (antibodies) by clonal selection, mutation
#' and network suppression, then assigns each observation to its nearest
#' prototype. By default the number of clusters emerges from suppression. Set
#' \code{targetK} to force exactly that many clusters. For mass cytometry,
#' \code{scale = "arcsinh"} and \code{affinityParams = list(alpha = "auto")}
#' are good starting points. Progress output is off unless \code{verbose =
#' TRUE} is passed through \code{...}.
#'
#' @return \code{bHIVEParam()} returns a \code{bHIVEParam} object, a subclass
#'   of \linkS4class{BlusterParam}.
#'
#'   \code{clusterRows()} returns a factor of cluster labels, one per row of
#'   \code{x}. With \code{full = TRUE} it returns a list with
#'   \code{clusters} and \code{objects}, where \code{objects$bHIVE} is the
#'   \code{bHIVE()} result, including the fitted \code{\link{AINet}} model
#'   in \code{objects$bHIVE$model} for labelling new cells.
#'
#' @seealso \code{\link{bHIVE}}, \code{\link[bluster]{clusterRows}},
#'   \code{\link[bluster]{clusterSweep}}.
#'
#' @examples
#' library(bluster)
#' X <- as.matrix(iris[, 1:4])
#'
#' p <- bHIVEParam(targetK = 3, maxIter = 10)
#' p
#' set.seed(1)
#' clusters <- clusterRows(X, p)
#' table(clusters, iris$Species)
#'
#' # Full output keeps the fitted model for labelling new cells
#' out <- clusterRows(X, p, full = TRUE)
#' out$objects$bHIVE$model$predict(X[1:5, ])
#'
#' # Sweep over parameters
#' sweep <- clusterSweep(X, bHIVEParam(maxIter = 5), nAntibodies = c(10, 30))
#' sweep$parameters
#'
#' # Through scran on a SingleCellExperiment
#' if (requireNamespace("scran", quietly = TRUE)) {
#'     library(SingleCellExperiment)
#'     sce <- SingleCellExperiment(assays = list(logcounts = t(X)))
#'     reducedDim(sce, "PCA") <- prcomp(X)$x[, 1:3]
#'     colLabels(sce) <- scran::clusterCells(sce, use.dimred = "PCA",
#'                                           BLUSPARAM = p)
#'     table(colLabels(sce), iris$Species)
#' }
#'
#' @name bHIVEParam
#' @aliases bHIVEParam-class
NULL

setClassUnion("bHIVE_numeric_OR_NULL", c("numeric", "NULL"))
setClassUnion("bHIVE_character_OR_NULL", c("character", "NULL"))

#' @export
#' @importClassesFrom bluster BlusterParam
#' @importFrom methods setClass setClassUnion setValidity new
setClass("bHIVEParam", contains = "BlusterParam",
         slots = c(nAntibodies     = "bHIVE_numeric_OR_NULL",
                   beta            = "bHIVE_numeric_OR_NULL",
                   epsilon         = "bHIVE_numeric_OR_NULL",
                   maxIter         = "bHIVE_numeric_OR_NULL",
                   k               = "bHIVE_numeric_OR_NULL",
                   affinityFunc    = "bHIVE_character_OR_NULL",
                   distFunc        = "bHIVE_character_OR_NULL",
                   initMethod      = "bHIVE_character_OR_NULL",
                   scale           = "bHIVE_character_OR_NULL",
                   targetK         = "bHIVE_numeric_OR_NULL",
                   epsilonQuantile = "bHIVE_numeric_OR_NULL",
                   affinityParams  = "list",
                   extra           = "list"))

#' @export
#' @rdname bHIVEParam
bHIVEParam <- function(nAntibodies = NULL, beta = NULL, epsilon = NULL,
                       maxIter = NULL, k = NULL, affinityFunc = NULL,
                       distFunc = NULL, initMethod = NULL, scale = NULL,
                       targetK = NULL, epsilonQuantile = NULL,
                       affinityParams = list(), ...) {
  new("bHIVEParam", nAntibodies = nAntibodies, beta = beta,
      epsilon = epsilon, maxIter = maxIter, k = k,
      affinityFunc = affinityFunc, distFunc = distFunc,
      initMethod = initMethod, scale = scale, targetK = targetK,
      epsilonQuantile = epsilonQuantile, affinityParams = affinityParams,
      extra = list(...))
}

# Allowed values for the character slots.
.BHIVE_PARAM_CHOICES <- list(
  affinityFunc = c("gaussian", "laplace", "polynomial", "cosine", "hamming"),
  distFunc     = c("euclidean", "manhattan", "minkowski", "cosine",
                   "mahalanobis", "hamming"),
  initMethod   = c("sample", "random", "random_uniform", "kmeans++"),
  scale        = c("none", "zscore", "robust", "arcsinh"))

# Lower bounds for the numeric slots ("gt" = strictly greater).
.BHIVE_PARAM_BOUNDS <- list(
  nAntibodies = c(ge = 1), beta = c(gt = 0), epsilon = c(ge = 0),
  maxIter = c(ge = 1), k = c(ge = 1), targetK = c(ge = 1))

#' Validity checks for bHIVEParam
#' @keywords internal
#' @noRd
.bhive_param_validity <- function(object) {
  msg <- character(0)
  for (s in names(.BHIVE_PARAM_BOUNDS)) {
    v <- methods::slot(object, s)
    b <- .BHIVE_PARAM_BOUNDS[[s]]
    if (is.null(v)) next
    ok <- length(v) == 1L && is.finite(v) &&
      if (names(b) == "gt") v > b else v >= b
    if (!ok) {
      msg <- c(msg, sprintf("'%s' must be a single number %s %s", s,
                            if (names(b) == "gt") ">" else ">=", b))
    }
  }
  q <- object@epsilonQuantile
  if (!is.null(q) && !(length(q) == 1L && is.finite(q) && q > 0 && q < 1)) {
    msg <- c(msg, "'epsilonQuantile' must be a single number in (0, 1)")
  }
  for (s in names(.BHIVE_PARAM_CHOICES)) {
    v <- methods::slot(object, s)
    if (!is.null(v) && !(length(v) == 1L && v %in% .BHIVE_PARAM_CHOICES[[s]])) {
      msg <- c(msg, sprintf("'%s' must be one of: %s", s,
                            paste(.BHIVE_PARAM_CHOICES[[s]], collapse = ", ")))
    }
  }
  msg <- c(msg, .bhive_param_check_lists(object))
  if (length(msg)) msg else TRUE
}

#' Check that affinityParams and extra are uniquely named and do not clash
#' @keywords internal
#' @noRd
.bhive_param_check_lists <- function(object) {
  msg <- character(0)
  named <- function(l) {
    length(l) == 0L ||
      (!is.null(names(l)) && all(nzchar(names(l))) && !anyDuplicated(names(l)))
  }
  if (!named(object@affinityParams)) {
    msg <- c(msg, "'affinityParams' must be a uniquely named list")
  }
  extra <- object@extra
  if (!named(extra)) {
    msg <- c(msg, "arguments in '...' must be uniquely named")
  }
  clash <- intersect(names(extra),
                     c(methods::slotNames("bHIVEParam"), "X", "y", "task"))
  if (length(clash)) {
    msg <- c(msg, sprintf("'...' cannot set: %s",
                          paste(clash, collapse = ", ")))
  }
  msg
}

setValidity("bHIVEParam", .bhive_param_validity)

#' @export
#' @rdname bHIVEParam
#' @importFrom bluster .defaultScalarArguments
#' @importFrom methods callNextMethod
setMethod(".defaultScalarArguments", "bHIVEParam", function(x) {
  c(callNextMethod(),
    nAntibodies = "numeric", beta = "numeric", epsilon = "numeric",
    maxIter = "numeric", k = "numeric", affinityFunc = "character",
    distFunc = "character", initMethod = "character", scale = "character",
    targetK = "numeric", epsilonQuantile = "numeric")
})

#' @export
#' @rdname bHIVEParam
#' @importFrom bluster .showScalarArguments
#' @importFrom methods show
setMethod("show", "bHIVEParam", function(object) {
  callNextMethod()
  .showScalarArguments(object)
  show_list <- function(label, l) {
    cat(sprintf("%s(%i): %s\n", label, length(l),
                paste(names(l), collapse = " ")))
  }
  show_list("affinityParams", object@affinityParams)
  show_list("extra", object@extra)
})

#' @export
#' @rdname bHIVEParam
#' @importFrom bluster clusterRows .extractScalarArguments
setMethod("clusterRows", c("ANY", "bHIVEParam"),
          function(x, BLUSPARAM, full = FALSE) {
  args <- c(list(X = as.matrix(x), task = "clustering", verbose = FALSE),
            .extractScalarArguments(BLUSPARAM))
  if (length(BLUSPARAM@affinityParams)) {
    args$affinityParams <- BLUSPARAM@affinityParams
  }
  # Extra arguments can override the defaults above (e.g. verbose = TRUE).
  args[names(BLUSPARAM@extra)] <- BLUSPARAM@extra
  res <- do.call(.bHIVE_matrix, args)
  clusters <- factor(res$assignments)
  if (full) {
    list(clusters = clusters, objects = list(bHIVE = res))
  } else {
    clusters
  }
})
