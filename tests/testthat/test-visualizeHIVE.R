# Tests for visualizeHIVE()

X <- as.matrix(iris[, 1:4])

set.seed(1)
res_clust <- bHIVE(X, targetK = 3, maxIter = 5, verbose = FALSE)
set.seed(2)
res_class <- bHIVE(X, y = iris$Species, nAntibodies = 15, maxIter = 5,
                   verbose = FALSE)
set.seed(3)
res_honey <- honeycombHIVE(X, task = "clustering", layers = 2,
                           nAntibodies = 20, maxIter = 5, verbose = FALSE)

proto_layer <- function(p) p$layers[[2]]$data
point_layer <- function(p) p$layers[[1]]$data

# ==== Scatter: one coordinate system for data and prototypes ====

test_that("PCA projects the prototypes onto the PCA of X", {
  p <- visualizeHIVE(res_clust, X = X, plot_type = "scatter",
                     transformation_method = "PCA")
  expect_s3_class(p, "ggplot")
  pr <- prcomp(X, scale. = TRUE)
  expected <- predict(pr, res_clust$antibodies_unscaled)[, 1:2]
  got <- as.matrix(proto_layer(p)[, c("Dim1", "Dim2")])
  expect_equal(unname(got), unname(expected))
  expect_equal(unname(as.matrix(point_layer(p)[, c("Dim1", "Dim2")])),
               unname(pr$x[, 1:2]))
  expect_identical(p$labels$x, "PC1")
})

test_that("bHIVE points are colored by their cluster or class", {
  p <- visualizeHIVE(res_clust, X = X)
  expect_false(anyNA(point_layer(p)$Group))
  expect_equal(nlevels(point_layer(p)$Group), 3L)
  p <- visualizeHIVE(res_class, X = X)
  expect_setequal(levels(point_layer(p)$Group), levels(iris$Species))
})

test_that("honeycomb layers use their membership and collapsed prototypes", {
  p <- visualizeHIVE(res_honey, X = X, layer = 1:2)
  pts <- point_layer(p)
  expect_equal(nrow(pts), 2 * nrow(X))
  expect_equal(as.integer(as.character(pts$Group[pts$Layer == "Layer 2"])),
               res_honey[[2]]$membership)
  expect_equal(nrow(proto_layer(p)),
               nrow(res_honey[[1]]$prototypes) +
                 nrow(res_honey[[2]]$prototypes))
  expect_true(inherits(p$facet, "FacetWrap"))
})

test_that("transformation 'none' plots the first two features as numbers", {
  p <- visualizeHIVE(res_clust, X = as.data.frame(X),
                     transformation_method = "none")
  expect_type(proto_layer(p)$Dim1, "double")
  expect_equal(proto_layer(p)$Dim1, unname(res_clust$antibodies_unscaled[, 1]))
  expect_identical(p$labels$x, "Sepal.Length")
})

test_that("tSNE embeds data and prototypes jointly despite duplicates", {
  skip_if_not_installed("Rtsne")
  set.seed(4)
  p <- visualizeHIVE(res_clust, X = X, transformation_method = "tSNE")
  expect_equal(nrow(point_layer(p)), nrow(X))
  expect_equal(nrow(proto_layer(p)), nrow(res_clust$antibodies))
})

test_that("UMAP projects the prototypes with predict()", {
  skip_if_not_installed("umap")
  set.seed(5)
  p <- visualizeHIVE(res_clust, X = X, transformation_method = "UMAP")
  expect_equal(nrow(proto_layer(p)), nrow(res_clust$antibodies))
  expect_false(anyNA(proto_layer(p)$Dim1))
})

test_that("extra arguments reach the reduction", {
  p <- visualizeHIVE(res_clust, X = X, scale. = FALSE)
  pr <- prcomp(X, scale. = FALSE)
  expect_equal(unname(as.matrix(point_layer(p)[, c("Dim1", "Dim2")])),
               unname(pr$x[, 1:2]))
})

# ==== Feature plots ====

test_that("feature plots work for matrices with named or indexed features", {
  for (type in c("boxplot", "violin", "density")) {
    p <- visualizeHIVE(res_class, X = X, plot_type = type,
                       feature = "Petal.Length", title = "T")
    expect_s3_class(p, "ggplot")
    expect_identical(p$labels$title, "T : Petal.Length")
  }
  p <- visualizeHIVE(res_clust, X = X, plot_type = "boxplot", feature = 3)
  expect_equal(proto_layer(p)$Prototype,
               unname(res_clust$antibodies_unscaled[, 3]))
})

test_that("prototype groups line up with the data groups", {
  p <- visualizeHIVE(res_class, X = X, plot_type = "boxplot")
  expect_true(all(as.character(proto_layer(p)$Group) %in%
                    levels(iris$Species)))
  expect_identical(as.character(proto_layer(p)$Group),
                   res_class$antibody_classes)
})

test_that("multi-layer feature plots keep each layer's own groups", {
  p <- visualizeHIVE(res_honey, X = X, plot_type = "boxplot", layer = 1:2)
  d <- p$data
  expect_equal(nrow(d), 2 * nrow(X))
  expect_equal(as.integer(as.character(d$Group[d$Layer == "Layer 1"])),
               res_honey[[1]]$membership)
})

test_that("classification honeycomb layers are grouped by predicted class", {
  set.seed(6)
  res <- honeycombHIVE(X, iris$Species, task = "classification",
                       layers = 2, maxIter = 5, verbose = FALSE)
  p <- visualizeHIVE(res, X = X, layer = 2, plot_type = "violin")
  expect_setequal(as.character(unique(p$data$Group)),
                  unique(res[[2]]$predictions))
})

# ==== Errors ====

test_that("visualizeHIVE gives clear errors", {
  expect_error(visualizeHIVE(res_clust), "X must be provided")
  expect_error(visualizeHIVE(res_clust, X = X, layer = "one"),
               "'layer' must be numeric")
  expect_error(visualizeHIVE(res_honey, X = X, layer = 3),
               "between 1 and 2")
  expect_error(visualizeHIVE(list(1, 2), X = X), "Invalid result")
  expect_error(visualizeHIVE(res_clust, X = X[1:10, ]),
               "per-observation labels")
  expect_error(visualizeHIVE(res_clust, X = X[, 1:3]),
               "different numbers of features")
  expect_error(visualizeHIVE(res_clust, X = X, plot_type = "boxplot",
                             feature = "nope"), "not a column")
  expect_error(visualizeHIVE(res_clust, X = X, plot_type = "boxplot",
                             feature = 9), "between 1 and 4")
  expect_error(visualizeHIVE(res_clust, X = X[, 1, drop = FALSE]),
               "at least two columns")
  expect_error(.require_suggested("notARealPackage123", "UMAP"),
               "needs the 'notARealPackage123' package")
})

test_that("the color legend is dropped when there are many groups", {
  p <- visualizeHIVE(res_honey, X = X, layer = 1)
  expect_gt(nlevels(point_layer(p)$Group), 12)
  expect_identical(p$guides$guides$colour, "none")
})
