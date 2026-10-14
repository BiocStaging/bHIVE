#' Visualize bHIVE/honeycombHIVE Results
#'
#' Plots the data together with the prototypes learned by \code{\link{bHIVE}}
#' or \code{\link{honeycombHIVE}}. One or more layers can be shown, faceted by
#' layer:
#'
#' \itemize{
#'   \item \code{"scatter"}: The observations, colored by cluster or class,
#'   with the prototypes drawn on top in the same coordinates. When \code{X}
#'   has more than two columns and \code{transform} is \code{TRUE}, a
#'   dimensionality reduction is fitted to \code{X} and the prototypes are
#'   projected into it.
#'   \item \code{"boxplot"}, \code{"violin"}: A selected feature by group, with
#'   the prototypes' values overlaid.
#'   \item \code{"density"}: Density of a selected feature by group, with the
#'   prototypes' values marked.
#' }
#'
#' For scatterplots the reduction is \code{"PCA"}, \code{"UMAP"},
#' \code{"tSNE"} or \code{"none"}. PCA and UMAP are fitted to \code{X} only and
#' the prototypes are placed with the fitted model's \code{predict()} method,
#' so both share one coordinate system. UMAP's projection of new points is an
#' approximation. tSNE has no such projection, so the data and prototypes are
#' embedded together and the prototypes take a small part in shaping the
#' layout.
#'
#' The prototypes drawn are, for a \code{bHIVE()} result, the antibodies in
#' original feature units (\code{antibodies_unscaled}), and for a
#' \code{honeycombHIVE()} layer, the collapsed \code{prototypes} that summarize
#' each cluster.
#'
#' @param result A list produced by \code{bHIVE} or \code{honeycombHIVE}.
#' For multilayer results, each element represents one layer.
#' @param X A numeric matrix or data frame of the original input features,
#' with the observations in the same order as the input to \code{bHIVE} or
#' \code{honeycombHIVE}. Required.
#' @param plot_type Character string specifying the type of plot to generate.
#' Options are:
#'   \itemize{
#'     \item \code{"scatter"}: A scatterplot of data points and prototypes.
#'     \item \code{"boxplot"}: A boxplot of a selected feature by group with
#'     prototypes overlaid.
#'     \item \code{"violin"}: A violin plot of a selected feature by group
#'     with prototypes overlaid.
#'     \item \code{"density"}: Density plots of a selected feature by group
#'     with prototype markers.
#'   }
#' @param feature Optional. For \code{"boxplot"}, \code{"violin"}, or
#' \code{"density"} plots, the name or index of the feature
#' in \code{X} to display. If \code{NULL}, the first column is used.
#' @param transform Logical. If \code{TRUE} and \code{X} has more than two
#' columns, the dimensionality reduction in \code{transformation_method} is
#' applied for scatterplots. Otherwise the first two columns are plotted.
#' @param transformation_method Character. The method used for dimensionality
#' reduction. Options are \code{"PCA"}, \code{"UMAP"} (needs the \pkg{umap}
#' package), \code{"tSNE"} (needs the \pkg{Rtsne} package), or \code{"none"}.
#' @param title Character. Title for the plot.
#' @param layer Integer vector indicating which layer(s) of the result to visualize.
#'   For bHIVE outputs, default is \code{1}. For multilayer honeycombHIVE
#'   outputs, specify one or more layer indices.
#' @param task Character, \code{"clustering"} or \code{"classification"}, or
#'   \code{NULL} (default) to use the task stored in \code{result}. Decides
#'   whether observations are grouped by cluster or by predicted class.
#' @param ... Additional arguments passed to the dimensionality reduction:
#'   \code{stats::prcomp()}, \code{umap::umap()} or \code{Rtsne::Rtsne()}.
#'
#' @return A \code{ggplot} object representing the visualization.
#'
#' @examples
#' data(iris)
#' X <- as.matrix(iris[, 1:4])
#'
#' # A single bHIVE clustering, with the prototypes projected onto the PCA
#' set.seed(1)
#' res <- bHIVE(X, targetK = 3, maxIter = 10, verbose = FALSE)
#' visualizeHIVE(res, X = X, plot_type = "scatter", title = "bHIVE prototypes")
#'
#' # Layers 1 and 2 of a honeycombHIVE hierarchy
#' res_h <- honeycombHIVE(X = X,
#'                        task = "clustering",
#'                        epsilon = 0.05,
#'                        layers = 3,
#'                        nAntibodies = 30,
#'                        beta = 5,
#'                        maxIter = 10,
#'                        verbose = FALSE)
#' visualizeHIVE(res_h, X = X, plot_type = "scatter", layer = 1:2,
#'               title = "honeycombHIVE layers")
#'
#' # A feature by class for a classifier, with each antibody's value overlaid
#' res_c <- bHIVE(X, y = iris$Species, maxIter = 10, verbose = FALSE)
#' visualizeHIVE(res_c, X = X, plot_type = "violin",
#'               feature = "Petal.Length", title = "Classifier")
#'
#' @importFrom stats prcomp predict
#' @import ggplot2
#' @export
visualizeHIVE <- function(result, X = NULL,
                          plot_type = c("scatter", "boxplot", "violin",
                                        "density"),
                          feature = NULL, transform = TRUE,
                          transformation_method = c("PCA", "UMAP", "tSNE",
                                                    "none"),
                          title = "HIVE Results", layer = 1, task = NULL,
                          ...) {
  plot_type <- match.arg(plot_type)
  method    <- match.arg(transformation_method)
  sel <- .viz_select_layers(result, layer)
  if (is.null(X)) {
    stop("X must be provided: visualizeHIVE() draws the data together ",
         "with the prototypes.", call. = FALSE)
  }
  X <- .viz_as_matrix(X)
  if (is.null(task)) task <- sel$layers[[1]]$task %||% "clustering"
  task <- match.arg(task, c("clustering", "classification"))

  groups <- lapply(sel$layers, function(l) .viz_point_groups(l, nrow(X), task))
  protos <- lapply(sel$layers, function(l) .viz_prototypes(l, task, X))
  if (plot_type == "scatter") {
    if (!transform || ncol(X) <= 2) method <- "none"
    .viz_scatter(X, groups, protos, sel$labels, method, title, list(...))
  } else {
    .viz_feature_plot(X, groups, protos, sel$labels, feature, plot_type,
                      title)
  }
}


# ============================================================================
# Internal helpers
# ============================================================================

#' Pick the layers to plot from a bHIVE or honeycombHIVE result
#' @return list(layers, labels)
#' @keywords internal
#' @noRd
.viz_select_layers <- function(result, layer) {
  if (!is.numeric(layer)) stop("'layer' must be numeric (layer indices).")
  if (!is.null(result$antibodies)) {
    return(list(layers = list(result), labels = "Layer 1"))
  }
  is_layer <- function(x) is.list(x) && !is.null(x$antibodies)
  if (!is.list(result) || length(result) == 0 ||
      !all(vapply(result, is_layer, logical(1)))) {
    stop("Invalid result object format.")
  }
  if (any(layer < 1 | layer > length(result))) {
    stop(sprintf("'layer' must be between 1 and %d.", length(result)))
  }
  list(layers = result[layer], labels = paste("Layer", layer))
}

#' Coerce X to a numeric matrix
#' @keywords internal
#' @noRd
.viz_as_matrix <- function(X) {
  X <- as.matrix(X)
  if (!is.numeric(X)) stop("X must be numeric.", call. = FALSE)
  if (ncol(X) < 2) {
    stop("X must have at least two columns for visualization.", call. = FALSE)
  }
  X
}

#' Per-observation groups for one layer
#'
#' Clustering: honeycomb membership, else the bHIVE cluster assignments.
#' Classification: honeycomb predictions, else the bHIVE class assignments.
#' Numeric IDs are ordered numerically.
#' @keywords internal
#' @noRd
.viz_point_groups <- function(layer, n, task) {
  candidates <- if (task == "classification") {
    list(layer$predictions, layer$assignments)
  } else {
    list(layer$membership, layer$assignments)
  }
  for (g in candidates) {
    if (!is.null(g) && length(g) == n) {
      if (is.factor(g)) return(g)
      return(factor(g, levels = sort(unique(g))))
    }
  }
  stop("The result has no per-observation labels matching nrow(X). ",
       "Is X the data the result was fitted on?", call. = FALSE)
}

#' Prototype matrix (in X units) and the group of each prototype
#' @return list(P, group)
#' @keywords internal
#' @noRd
.viz_prototypes <- function(layer, task, X) {
  if (!is.null(layer$prototypes)) {
    # honeycomb: row k summarizes the observations with membership == k
    P <- as.matrix(layer$prototypes)
    k <- seq_len(nrow(P))
    group <- if (task == "classification" && !is.null(layer$predictions)) {
      as.character(layer$predictions[match(k, layer$membership)])
    } else {
      as.character(k)
    }
  } else {
    P <- as.matrix(layer$antibodies_unscaled %||% layer$antibodies)
    group <- if (task == "classification" && !is.null(layer$antibody_classes)) {
      as.character(layer$antibody_classes)
    } else {
      as.character(seq_len(nrow(P)))
    }
  }
  if (ncol(P) != ncol(X)) {
    stop("The prototypes and X have different numbers of features.",
         call. = FALSE)
  }
  colnames(P) <- colnames(X)
  list(P = P, group = group)
}

#' Embed X and project the prototypes P into the same 2-D coordinates
#' @param args Extra arguments for the reduction (from `...`).
#' @return list(X = n x 2, P = m x 2)
#' @keywords internal
#' @noRd
.viz_embed <- function(X, P, method, args) {
  two <- function(m) m[, 1:2, drop = FALSE]
  if (method == "none") return(list(X = two(X), P = two(P)))
  if (method == "PCA") {
    if (is.null(args$scale.)) args$scale. <- TRUE
    pr <- do.call(stats::prcomp, c(list(X), args))
    return(list(X = two(pr$x), P = two(predict(pr, P))))
  }
  if (method == "UMAP") {
    .require_suggested("umap", method)
    fit <- do.call(umap::umap, c(list(X), args))
    return(list(X = two(fit$layout), P = two(predict(fit, P))))
  }
  # tSNE has no out-of-sample projection: embed data and prototypes jointly.
  .require_suggested("Rtsne", method)
  joint <- rbind(X, P)
  if (is.null(args$perplexity)) {
    args$perplexity <- max(1, min(30, floor((nrow(joint) - 1) / 3)))
  }
  if (is.null(args$check_duplicates)) args$check_duplicates <- FALSE
  Y <- do.call(Rtsne::Rtsne, c(list(joint), args))$Y
  in_X <- seq_len(nrow(X))
  list(X = two(Y[in_X, , drop = FALSE]), P = two(Y[-in_X, , drop = FALSE]))
}

#' Stop with an install hint when a suggested package is missing
#' @keywords internal
#' @noRd
.require_suggested <- function(pkg, what) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop(sprintf("transformation_method = '%s' needs the '%s' package. ",
                 what, pkg),
         sprintf("Install it with install.packages('%s').", pkg),
         call. = FALSE)
  }
  invisible(TRUE)
}

#' Axis labels for each reduction
#' @keywords internal
#' @noRd
.viz_axis_labels <- function(method, X) {
  if (method == "none") {
    cn <- colnames(X)
    return(if (is.null(cn)) c("Feature 1", "Feature 2") else cn[1:2])
  }
  paste0(c(PCA = "PC", UMAP = "UMAP", tSNE = "tSNE")[[method]], 1:2)
}

#' Scatterplot of observations and prototypes in shared coordinates
#' @keywords internal
#' @noRd
.viz_scatter <- function(X, groups, protos, labels, method, title, args) {
  emb <- .viz_embed(X, do.call(rbind, lapply(protos, `[[`, "P")), method,
                    args)
  df_points <- do.call(rbind, lapply(seq_along(labels), function(i) {
    data.frame(Dim1 = emb$X[, 1], Dim2 = emb$X[, 2], Group = groups[[i]],
               Layer = labels[i])
  }))
  n_proto <- vapply(protos, function(p) nrow(p$P), integer(1))
  df_protos <- data.frame(Dim1 = emb$P[, 1], Dim2 = emb$P[, 2],
                          Layer = rep(labels, n_proto))
  axes <- .viz_axis_labels(method, X)
  p <- ggplot() +
    geom_point(data = df_points,
               aes(x = .data$Dim1, y = .data$Dim2, color = .data$Group),
               alpha = 0.6, size = 1.5) +
    geom_point(data = df_protos, aes(x = .data$Dim1, y = .data$Dim2),
               shape = 21, fill = "black", color = "white", size = 3.5,
               stroke = 0.8) +
    labs(title = title, x = axes[1], y = axes[2], color = "Group") +
    theme_minimal() +
    scale_color_viridis_d()
  # A legend with many near-identical shades is unreadable.
  if (nlevels(df_points$Group) > 12) p <- p + guides(color = "none")
  if (length(labels) > 1) p <- p + facet_wrap(~Layer)
  p
}

#' Values and display name of the selected feature
#' @keywords internal
#' @noRd
.viz_feature_values <- function(X, feature) {
  if (is.null(feature)) feature <- 1L
  if (is.character(feature)) {
    if (!feature %in% colnames(X)) {
      stop(sprintf("feature '%s' is not a column of X.", feature),
           call. = FALSE)
    }
    index <- match(feature, colnames(X))
  } else {
    index <- as.integer(feature)
    if (index < 1 || index > ncol(X)) {
      stop(sprintf("feature must be between 1 and %d.", ncol(X)),
           call. = FALSE)
    }
  }
  name <- colnames(X)[index]
  if (is.null(name)) name <- paste("Feature", index)
  list(index = index, name = name, values = X[, index])
}

#' Boxplot, violin or density of one feature by group, prototypes overlaid
#' @keywords internal
#' @noRd
.viz_feature_plot <- function(X, groups, protos, labels, feature, plot_type,
                              title) {
  f <- .viz_feature_values(X, feature)
  per_layer <- lapply(seq_along(labels), function(i) {
    lv <- unique(c(levels(groups[[i]]), protos[[i]]$group))
    list(data  = data.frame(Feature = f$values,
                            Group = factor(groups[[i]], levels = lv),
                            Layer = labels[i]),
         proto = data.frame(Prototype = protos[[i]]$P[, f$index],
                            Group = factor(protos[[i]]$group, levels = lv),
                            Layer = labels[i]))
  })
  df_feature <- do.call(rbind, lapply(per_layer, `[[`, "data"))
  df_proto   <- do.call(rbind, lapply(per_layer, `[[`, "proto"))
  if (plot_type == "density") {
    p <- ggplot(df_feature, aes(x = .data$Feature, fill = .data$Group)) +
      geom_density(alpha = 0.5) +
      geom_vline(data = df_proto, aes(xintercept = .data$Prototype),
                 linetype = "dashed", linewidth = 0.4) +
      labs(x = f$name, y = "Density", fill = "Group")
  } else {
    body <- if (plot_type == "boxplot") geom_boxplot(alpha = 0.7) else
      geom_violin(alpha = 0.7)
    p <- ggplot(df_feature, aes(x = .data$Group, y = .data$Feature,
                                fill = .data$Group)) +
      body +
      geom_point(data = df_proto, aes(x = .data$Group, y = .data$Prototype),
                 inherit.aes = FALSE, shape = 21, fill = "black",
                 color = "white", size = 3) +
      labs(x = "Group", y = f$name)
  }
  p <- p + labs(title = paste(title, ":", f$name)) + theme_minimal() +
    scale_fill_viridis_d()
  if (plot_type != "density") p <- p + theme(legend.position = "none")
  if (length(labels) > 1) {
    p <- p + facet_wrap(~Layer, scales = if (plot_type == "density") "fixed" else "free_x")
  }
  p
}
