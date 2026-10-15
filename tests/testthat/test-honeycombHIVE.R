# test script for honeycombHIVE.R - testcases are NOT comprehensive!

set.seed(42)
test_that("honeycombHIVE runs successfully for clustering task", {
  data(iris)
  X <- as.matrix(iris[, 1:4])
  
  # Run honeycombHIVE for clustering
  res <- honeycombHIVE(X = X, 
                       task = "clustering", 
                       layers = 3, 
                       nAntibodies = 15, 
                       beta = 5, 
                       maxIter = 5,
                       verbose = FALSE)
  
  # Test structure and results
  expect_type(res, "list")
  expect_length(res, 3) # 3 layers
  for (layer in res) {
    expect_true(all(c("antibodies", "assignments", "task", "membership") %in%
                      names(layer)))
    expect_equal(layer$task, "clustering")
    expect_true(is.matrix(layer$antibodies))
    expect_true(is.numeric(layer$assignments))
    expect_equal(length(layer$membership), 150)
  }
})

test_that("honeycombHIVE runs successfully for classification task", {
  data(iris)
  X <- as.matrix(iris[, 1:4])
  y <- iris$Species
  
  # Run honeycombHIVE for classification
  res <- honeycombHIVE(X = X, 
                       y = as.factor(y), 
                       task = "classification", 
                       layers = 2, 
                       nAntibodies = 10, 
                       beta = 3, 
                       maxIter = 5, 
                       verbose = FALSE)
  
  # Test structure and results
  expect_type(res, "list")
  expect_length(res, 2) # 2 layers
  for (layer in res) {
    expect_true(all(c("antibodies", "assignments", "task", "predictions",
                      "membership") %in% names(layer)))
    expect_equal(layer$task, "classification")
    expect_true(is.matrix(layer$antibodies))
  }
})

test_that("honeycombHIVE handles different affinity functions", {
  data(iris)
  X <- as.matrix(iris[, 1:4])
  
  natural <- c(gaussian = "euclidean", laplace = "manhattan",
               polynomial = "euclidean", cosine = "cosine")
  for (aff in names(natural)) {
    expect_no_error(
      res <- honeycombHIVE(X = X, 
                           task = "clustering", 
                           layers = 2, 
                           nAntibodies = 10, 
                           affinityFunc = aff, 
                           distFunc = natural[[aff]], 
                           maxIter = 5,
                           verbose = FALSE)
    )
    expect_length(res, 2)
    for (layer in res) {
      expect_equal(layer$task, "clustering")
    }
  }
})

test_that("honeycombHIVE handles different distance functions", {
  data(iris)
  X <- as.matrix(iris[, 1:4])
  
  dist_funcs <- c("euclidean", "manhattan", "minkowski")
  for (dist in dist_funcs) {
    expect_no_error(
      res <- honeycombHIVE(X = X, 
                           task = "clustering", 
                           layers = 2, 
                           nAntibodies = 10, 
                           distFunc = dist, 
                           maxIter = 5,
                           verbose = FALSE)
    )
    expect_length(res, 2)
    for (layer in res) {
      expect_equal(layer$task, "clustering")
    }
  }
})

test_that("honeycombHIVE refine=TRUE for classification with cross-entropy", {
  data(iris)
  X <- as.matrix(iris[,1:4])
  y <- iris$Species
  
  # Compare refine=FALSE vs. refine=TRUE
  res_no_ref <- honeycombHIVE(
    X = X,
    y = y,
    task = "classification",
    layers = 2,
    nAntibodies = 10,
    maxIter = 5,
    verbose = FALSE,
    refine = FALSE
  )
  
  res_ref <- honeycombHIVE(
    X = X,
    y = y,
    task = "classification",
    layers = 2,
    nAntibodies = 10,
    maxIter = 5,
    verbose = FALSE,
    refine = TRUE,
    refineLoss = "categorical_crossentropy",
    refineSteps = 3,
    refineLR = 0.05,
    refinePushAway = TRUE
  )
  
  # Check structure
  expect_length(res_ref, 2)
  for (layer_obj in res_ref) {
    expect_true(all(c("antibodies","assignments","task",
                      "predictions","membership") %in% names(layer_obj)))
    expect_true(is.matrix(layer_obj$antibodies))
  }
  
  # Check that the final prototypes differ
  final_no_ref <- res_no_ref[[2]]$antibodies
  final_ref    <- res_ref[[2]]$antibodies
  expect_false(isTRUE(all.equal(final_no_ref, final_ref)))
  
  preds_no_ref <- res_no_ref[[2]]$predictions
  preds_ref    <- res_ref[[2]]$predictions
})

set.seed(42)
test_that("honeycombHIVE refineSteps=0 does not change prototypes", {
  data(iris)
  X <- as.matrix(iris[,1:4])
  
  # Clustering is stochastic (and now runs Lloyd consolidation, so the cluster
  # count itself varies run-to-run). Seed each call identically so all three
  # share the same underlying partition and prototype set; the ONLY difference
  # is the refinement, which is what this test isolates.

  # refine=TRUE but steps=0 => effectively no gradient updates
  set.seed(123)
  res_steps0 <- honeycombHIVE(
    X = X,
    task = "clustering",
    layers = 1,
    nAntibodies = 8,
    maxIter = 3,
    verbose = FALSE,
    refine = TRUE,
    refineLoss = "mae",
    refineSteps = 0,
    refineLR = 0.1
  )

  # normal refineSteps>0
  set.seed(123)
  res_steps5 <- honeycombHIVE(
    X = X,
    task = "clustering",
    layers = 1,
    nAntibodies = 8,
    maxIter = 3,
    verbose = FALSE,
    refine = TRUE,
    refineLoss = "mae",
    refineSteps = 5,
    refineLR = 0.1
  )

  # Compare final prototypes
  set.seed(123)
  no_refine <- honeycombHIVE(
    X = X,
    task = "clustering",
    layers = 1,
    nAntibodies = 8,
    maxIter = 3,
    verbose = FALSE,
    refine = FALSE
  )
  
  prot_steps0 <- res_steps0[[1]]$antibodies
  prot_steps5 <- res_steps5[[1]]$antibodies
  prot_no_ref <- no_refine[[1]]$antibodies
  
  diff0 <- sum(abs(as.matrix(prot_steps0) - as.matrix(prot_no_ref)))
  diff5 <- sum(abs(as.matrix(prot_steps5) - as.matrix(prot_no_ref)))
  expect_true(diff0 < diff5,
              info="Refinement steps=0 should yield prototypes closer to no_refine than steps=5.")
})

test_that("honeycombHIVE fails gracefully with invalid refineLoss for given task", {
  data(iris)
  X <- as.matrix(iris[,1:4])
  y <- iris$Species

  expect_error({
    honeycombHIVE(
      X = X,
      y = y,
      task = "classification",
      layers = 1,
      refine = TRUE,
      refineLoss = "bogus_loss",
      refineSteps = 3,
      refineLR = 0.01,
      verbose = FALSE
    )
  }, regexp="Invalid refineLoss 'bogus_loss' for task 'classification'")
})

test_that("honeycombHIVE clustering works with refine=TRUE, default refineLoss=mae", {
  data(iris)
  X <- as.matrix(iris[, 1:4])
  
  # refine=FALSE => normal honeycombHIVE
  res_no_ref <- honeycombHIVE(
    X = X,
    task = "clustering",
    layers = 2,
    nAntibodies = 8,
    maxIter = 5,
    refine = FALSE,
    verbose = FALSE
  )
  
  # refine=TRUE => add refine step
  res_ref <- honeycombHIVE(
    X = X,
    task = "clustering",
    layers = 2,
    nAntibodies = 8,
    maxIter = 5,
    refine = TRUE,            
    refineLoss = "mae",       
    refineSteps = 3,
    refineLR = 0.01,
    verbose = FALSE
  )
  
  # Check structure
  expect_length(res_ref, 2)  # 2 layers
  for (layer_obj in res_ref) {
    expect_true("antibodies" %in% names(layer_obj))
    expect_true("assignments" %in% names(layer_obj))
    expect_true("membership"  %in% names(layer_obj))
    expect_equal(layer_obj$task, "clustering")
  }
  
  # Compare final prototypes 
  final_no_ref <- res_no_ref[[2]]$antibodies
  final_ref    <- res_ref[[2]]$antibodies
  # Usually not identical
  expect_false(isTRUE(all.equal(final_no_ref, final_ref)))
})

test_that("honeycombHIVE classification: refine=TRUE with cross_entropy", {
  data(iris)
  X <- as.matrix(iris[, 1:4])
  y <- iris$Species
  
  # Non-refined
  res_no_ref <- honeycombHIVE(
    X = X,
    y = y,
    task = "classification",
    layers = 2,
    nAntibodies = 10,
    maxIter = 5,
    refine = FALSE,
    verbose = FALSE
  )
  
  # Refined with cross-entropy
  res_ref <- honeycombHIVE(
    X = X,
    y = y,
    task = "classification",
    layers = 2,
    nAntibodies = 10,
    maxIter = 5,
    refine = TRUE,
    refineLoss = "categorical_crossentropy",
    refineSteps = 4,
    refineLR = 0.01,
    refinePushAway = TRUE,
    verbose = FALSE
  )
  
  # Structure checks
  expect_length(res_ref, 2)
  for (ly in res_ref) {
    expect_equal(ly$task, "classification")
    expect_true(is.matrix(ly$antibodies))
    # We also expect "assignments", "predictions", "membership"
    expect_true("assignments" %in% names(ly))
    expect_true("membership"  %in% names(ly))
  }
  
  # Compare final layer prototypes
  final_no_ref <- res_no_ref[[2]]$antibodies
  final_ref    <- res_ref[[2]]$antibodies
  expect_false(isTRUE(all.equal(final_no_ref, final_ref)))
})

test_that("honeycombHIVE refineSteps=0 yields same result as refine=FALSE", {
  data(iris)
  X <- as.matrix(iris[,1:4])

  set.seed(42)
  # refine=FALSE
  res_no_ref <- honeycombHIVE(
    X = X,
    task = "clustering",
    layers = 1,
    nAntibodies = 6,
    maxIter = 3,
    refine = FALSE,
    verbose = FALSE
  )

  set.seed(42)
  res_steps0 <- honeycombHIVE(
    X = X,
    task = "clustering",
    layers = 1,
    nAntibodies = 6,
    maxIter = 3,
    refine = TRUE,
    refineLoss = "mae",
    refineSteps = 0,
    refineLR = 0.01,
    verbose = FALSE
  )

  # Compare final prototypes
  final_no_ref <- res_no_ref[[1]]$antibodies
  final_steps0 <- res_steps0[[1]]$antibodies

  expect_equal(final_no_ref, final_steps0, tolerance=1e-14)
})

# ==== Classification hierarchy, refinement space, final-layer model ====

test_that("classification keeps several prototypes per class past layer 1", {
  # Collapsing by class label used to leave 3 prototypes per layer, and the
  # layer-2 predictions lost their labels (all one class).
  X <- iris[, 1:4]
  y <- iris$Species
  set.seed(1)
  res <- honeycombHIVE(X, y, task = "classification", layers = 2,
                       verbose = FALSE)
  expect_gt(nrow(res[[1]]$prototypes), nlevels(y))
  for (layer in res) {
    expect_gt(mean(layer$predictions == y), 0.7)
    expect_equal(nrow(layer$prototypes), max(layer$membership))
  }
})

test_that("predictions are the majority label of each cluster", {
  X <- iris[, 1:4]
  y <- iris$Species
  set.seed(2)
  res <- honeycombHIVE(X, y, task = "classification", layers = 2,
                       verbose = FALSE)
  layer <- res[[2]]
  for (k in unique(layer$membership)) {
    idx <- which(layer$membership == k)
    expect_equal(unique(layer$predictions[idx]), .majority_label(y[idx]))
  }
})

test_that("only the final layer keeps the fitted model", {
  set.seed(3)
  res <- honeycombHIVE(iris[, 1:4], task = "clustering", layers = 3,
                       nAntibodies = 15, maxIter = 5, verbose = FALSE)
  expect_null(res[[1]]$model)
  expect_null(res[[2]]$model)
  expect_s3_class(res[[3]]$model, "AINet")
})

test_that("refinement runs in the model's scaled space", {
  X <- iris[, 1:4]
  set.seed(4)
  res <- honeycombHIVE(X, iris$Species, task = "classification", layers = 2,
                       maxIter = 5, refine = TRUE, refineSteps = 3,
                       scale = "zscore", verbose = FALSE)
  final <- res[[2]]
  # Unscaled prototypes stay in the data range, and the stored model uses
  # the refined antibodies.
  rng <- range(as.matrix(X))
  expect_true(all(final$antibodies_unscaled >= rng[1] - 1 &
                    final$antibodies_unscaled <= rng[2] + 1))
  expect_identical(final$model$result$antibodies, final$antibodies)
})

test_that("medoid prototypes keep the feature names", {
  set.seed(5)
  res <- honeycombHIVE(iris[, 1:4], task = "clustering", layers = 2,
                       nAntibodies = 10, maxIter = 3,
                       collapseMethod = "medoid", verbose = FALSE)
  expect_identical(colnames(res[[2]]$antibodies), colnames(iris)[1:4])
})

# ==== Small clusters and collapse methods ====

test_that("small clusters are merged, kept or dropped as requested", {
  X <- iris[, 1:4]
  fit <- function(action) {
    set.seed(30)
    honeycombHIVE(X, task = "clustering", layers = 1, nAntibodies = 30,
                  maxIter = 5, minClusterSize = 8,
                  smallClusterAction = action, verbose = FALSE)[[1]]
  }
  merged <- fit("merge")
  kept   <- fit("keep")
  dropped <- fit("drop")
  # merge keeps every observation and only clusters of at least 8
  expect_false(anyNA(merged$membership))
  expect_true(all(table(merged$membership) >= 8))
  # keep retains every cluster, including the small ones
  expect_gt(nrow(kept$prototypes), nrow(merged$prototypes))
  expect_false(anyNA(kept$membership))
  # drop discards small clusters and their members
  expect_true(anyNA(dropped$membership))
  expect_equal(nrow(dropped$prototypes), nrow(merged$prototypes))
})

test_that("small-cluster handling is reported with message()", {
  X <- iris[, 1:4]
  for (action in c("merge", "keep", "drop")) {
    set.seed(31)
    msgs <- capture_messages(
      honeycombHIVE(X, task = "clustering", layers = 1, nAntibodies = 30,
                    maxIter = 3, minClusterSize = 8,
                    smallClusterAction = action, verbose = TRUE))
    expect_true(any(grepl("small cluster", msgs)), label = action)
    expect_true(any(grepl("Layer 1 completed", msgs)), label = action)
  }
})

test_that("median and mode collapse methods summarize each cluster", {
  X <- round(iris[, 1:4])
  for (cm in c("median", "mode")) {
    set.seed(32)
    res <- honeycombHIVE(X, task = "clustering", layers = 1,
                         nAntibodies = 10, maxIter = 3, collapseMethod = cm,
                         verbose = FALSE)[[1]]
    k <- res$membership[1]
    members <- X[res$membership == k, , drop = FALSE]
    expected <- if (cm == "median") {
      apply(members, 2, median)
    } else {
      apply(members, 2, function(v) as.numeric(names(which.max(table(v)))))
    }
    expect_equal(unname(res$prototypes[k, ]), unname(expected), label = cm)
  }
})

test_that("honeycombHIVE rejects an unknown refinement loss", {
  expect_error(honeycombHIVE(iris[, 1:4], refine = TRUE, refineLoss = "mse",
                             verbose = FALSE), "Invalid refineLoss")
})
