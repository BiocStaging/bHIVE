# Tests for bHIVEParam, the bluster interface

X <- as.matrix(iris[, 1:4])

test_that("bHIVEParam constructs, shows and extends BlusterParam", {
  p <- bHIVEParam(targetK = 3, affinityParams = list(alpha = "auto"),
                  coverageBoost = TRUE)
  expect_s4_class(p, "bHIVEParam")
  expect_s4_class(p, "BlusterParam")
  out <- capture.output(show(p))
  expect_true("class: bHIVEParam" %in% out)
  expect_true("targetK: 3" %in% out)
  expect_true("nAntibodies: [default]" %in% out)
  expect_true("affinityParams(1): alpha" %in% out)
  expect_true("extra(1): coverageBoost" %in% out)
})

test_that("clusterRows matches bHIVE() under the same seed", {
  p <- bHIVEParam(nAntibodies = 15, maxIter = 5)
  set.seed(1)
  cl <- bluster::clusterRows(X, p)
  set.seed(1)
  ref <- bHIVE(X, nAntibodies = 15, maxIter = 5, verbose = FALSE)
  expect_s3_class(cl, "factor")
  expect_length(cl, nrow(X))
  expect_equal(as.numeric(as.character(cl)), ref$assignments)
})

test_that("full = TRUE returns the bHIVE result and fitted model", {
  set.seed(2)
  out <- bluster::clusterRows(X, bHIVEParam(targetK = 3, maxIter = 5),
                              full = TRUE)
  expect_named(out, c("clusters", "objects"))
  expect_equal(nlevels(out$clusters), 3L)
  model <- out$objects$bHIVE$model
  expect_s3_class(model, "AINet")
  expect_equal(as.numeric(model$predict(X)),
               as.numeric(as.character(out$clusters)))
})

test_that("extra arguments and affinityParams reach bHIVE", {
  set.seed(3)
  out <- bluster::clusterRows(
    X, bHIVEParam(maxIter = 3, affinityParams = list(alpha = "auto"),
                  coverageBoost = TRUE, scale = "zscore"),
    full = TRUE)
  cfg <- out$objects$bHIVE$model$config
  expect_true(cfg$coverageBoost)
  expect_equal(cfg$scale, "zscore")
  expect_true(is.numeric(cfg$kernel$alpha))
})

test_that("progress is silent unless verbose is requested", {
  p <- bHIVEParam(nAntibodies = 10, maxIter = 2)
  expect_silent(bluster::clusterRows(X, p))
  msgs <- capture_messages(
    bluster::clusterRows(X, bHIVEParam(nAntibodies = 10, maxIter = 2,
                                       verbose = TRUE)))
  expect_true(any(grepl("Iteration 1", msgs)))
})

test_that("bluster accessors and clusterSweep work", {
  p <- bHIVEParam(maxIter = 3)
  p[["epsilon"]] <- 0.2
  expect_equal(p[["epsilon"]], 0.2)
  set.seed(4)
  sw <- bluster::clusterSweep(X, p, nAntibodies = c(10, 30))
  expect_equal(ncol(sw$clusters), 2L)
  expect_equal(sw$parameters$nAntibodies, c(10, 30))
})

test_that("bHIVEParam validates its arguments", {
  expect_error(bHIVEParam(nAntibodies = -1), "nAntibodies")
  expect_error(bHIVEParam(beta = 0), "beta")
  expect_error(bHIVEParam(epsilonQuantile = 2), "epsilonQuantile")
  expect_error(bHIVEParam(affinityFunc = "rbf"), "affinityFunc")
  expect_error(bHIVEParam(scale = c("zscore", "robust")), "scale")
  expect_error(bHIVEParam(affinityParams = list(1)), "affinityParams")
  expect_error(bHIVEParam(X = X), "cannot set: X")
  expect_error(methods::new("bHIVEParam", extra = list(TRUE)), "uniquely named")
})

test_that("clusterRows labels a SingleCellExperiment from its reducedDim", {
  sce <- SingleCellExperiment::SingleCellExperiment(
    assays = list(logcounts = t(X)))
  SingleCellExperiment::reducedDim(sce, "PCA") <- prcomp(X)$x[, 1:3]
  p <- bHIVEParam(targetK = 3, maxIter = 5)
  set.seed(5)
  SingleCellExperiment::colLabels(sce) <- bluster::clusterRows(
    SingleCellExperiment::reducedDim(sce, "PCA"), p)
  expect_length(SingleCellExperiment::colLabels(sce), ncol(sce))
  expect_equal(nlevels(SingleCellExperiment::colLabels(sce)), 3L)
})
