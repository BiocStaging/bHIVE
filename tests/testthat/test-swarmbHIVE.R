# test script for swarmbHIVE.R - testcases are NOT comprehensive!
set.seed(42)

test_that("swarmbHIVE works for clustering (silhouette)", {
  data(iris)
  X <- as.matrix(iris[, 1:4])
  grid <- expand.grid(
    nAntibodies = c(10, 20),
    beta        = c(3, 5),
    epsilon     = c(0.01, 0.05)
  )
  
  # Tune hyperparameters
  res <- swarmbHIVE(
    X     = X,
    task  = "clustering",
    grid  = grid,
    metric = "silhouette",
    maxIter = 10,
    verbose = FALSE
  )
  
  # Check results structure
  expect_named(res, c("best_params", "results"))
  expect_true(is.data.frame(res$results))
  expect_equal(ncol(res$results), ncol(grid) + 1) # grid cols + metric_value
  expect_true(is.data.frame(res$best_params))
  expect_equal(nrow(res$best_params), 1)
})

test_that("swarmbHIVE works for classification (accuracy)", {
  data(iris)
  X <- as.matrix(iris[, 1:4])
  y <- iris$Species
  
  grid <- expand.grid(
    nAntibodies = c(10, 20),
    beta        = c(3, 5),
    epsilon     = c(0.01, 0.05)
  )
  
  # Tune hyperparameters
  res <- swarmbHIVE(
    X     = X, 
    y     = y, 
    task  = "classification",
    grid  = grid, 
    metric = "accuracy",
    maxIter = 10,
    verbose = FALSE
  )
  
  # Check results structure
  expect_named(res, c("best_params", "results"))
  expect_true(is.data.frame(res$results))
  expect_equal(ncol(res$results), ncol(grid) + 1)
  expect_true(is.data.frame(res$best_params))
  expect_equal(nrow(res$best_params), 1)
})

test_that("swarmbHIVE throws an error for invalid metric", {
  data(iris)
  X <- as.matrix(iris[, 1:4])
  grid <- expand.grid(
    nAntibodies = c(10, 20),
    beta        = c(3, 5),
    epsilon     = c(0.01, 0.05)
  )
  
  # Invalid metric
  expect_error(
    swarmbHIVE(
      X     = X,
      task  = "clustering",
      grid  = grid,
      metric = "invalid_metric"
    )
  )
})

test_that("swarmbHIVE works for classification (balanced_accuracy)", {
  data(iris)
  X <- as.matrix(iris[, 1:4])
  y <- iris$Species
  
  grid <- expand.grid(
    nAntibodies = 10,
    beta        = 3,
    epsilon     = 0.01
  )
  
  res <- swarmbHIVE(
    X         = X,
    y         = y,
    task      = "classification",
    grid      = grid,
    metric    = "balanced_accuracy",
    maxIter   = 5, 
    verbose   = FALSE
  )
  
  expect_named(res, c("best_params", "results"))
  expect_equal(nrow(res$best_params), 1)
  expect_true(all(res$results$metric_value <= 1))
})

test_that("swarmbHIVE works for classification (f1)", {
  data(iris)
  X <- as.matrix(iris[, 1:4])
  y <- iris$Species
  
  grid <- expand.grid(
    nAntibodies = 10,
    beta        = 3,
    epsilon     = 0.01
  )
  
  res <- swarmbHIVE(
    X      = X,
    y      = y,
    task   = "classification",
    grid   = grid,
    metric = "f1",
    maxIter= 5,
    verbose= FALSE
  )
  
  expect_named(res, c("best_params", "results"))
  expect_equal(nrow(res$best_params), 1)
  # F1 generally in [0,1]
  expect_true(all(res$results$metric_value >= 0 & res$results$metric_value <= 1))
})

test_that("swarmbHIVE works for classification (kappa)", {
  data(iris)
  X <- as.matrix(iris[, 1:4])
  y <- iris$Species
  
  grid <- expand.grid(
    nAntibodies = 10,
    beta        = 3,
    epsilon     = 0.01
  )
  
  res <- swarmbHIVE(
    X      = X,
    y      = y,
    task   = "classification",
    grid   = grid,
    metric = "kappa",
    maxIter= 5,
    verbose= FALSE
  )
  
  expect_named(res, c("best_params", "results"))
  expect_equal(nrow(res$best_params), 1)
  # Kappa can be negative in worst cases, but let's just check it doesn't break
  expect_true(all(!is.na(res$results$metric_value)))
})

test_that("swarmbHIVE works for clustering (davies_bouldin)", {
  data(iris)
  X <- as.matrix(iris[, 1:4])
  grid <- expand.grid(
    nAntibodies = 5,
    beta        = 3,
    epsilon     = 0.01
  )
  
  res <- swarmbHIVE(
    X       = X,
    task    = "clustering",
    grid    = grid,
    metric  = "davies_bouldin",
    maxIter = 5,
    verbose = FALSE
  )
  
  # DB index is smaller=better
  # We won't check absolute range, just that we got a real number
  expect_named(res, c("best_params", "results"))
  expect_equal(nrow(res$best_params), 1)
  expect_true(all(!is.na(res$results$metric_value)))
})

test_that("swarmbHIVE works for clustering (calinski_harabasz)", {
  data(iris)
  X <- as.matrix(iris[, 1:4])
  grid <- expand.grid(
    nAntibodies = 5,
    beta        = 3,
    epsilon     = 0.01
  )
  
  res <- swarmbHIVE(
    X       = X,
    task    = "clustering",
    grid    = grid,
    metric  = "calinski_harabasz",
    maxIter = 5,
    verbose = FALSE
  )
  
  expect_named(res, c("best_params", "results"))
  expect_equal(nrow(res$best_params), 1)
  expect_true(all(!is.na(res$results$metric_value)))
})

# ==== Grid forwarding, progress output and metric helpers ====

test_that("every tunable grid column reaches bHIVE, factors as character", {
  seen <- list()
  local_mocked_bindings(.bHIVE_matrix = function(...) {
    seen[[length(seen) + 1L]] <<- list(...)
    list(assignments = rep(1:2, length.out = 150))
  })
  X <- as.matrix(iris[, 1:4])
  grid <- expand.grid(affinityFunc = c("gaussian", "laplace"),
                      coverageBoost = TRUE, maxIter = 7, junk = 1)
  swarmbHIVE(X, task = "clustering", grid = grid, metric = "silhouette",
             maxIter = 50, verbose = FALSE)
  expect_length(seen, 2L)
  expect_identical(seen[[2]]$affinityFunc, "laplace")
  expect_true(seen[[1]]$coverageBoost)
  # A maxIter grid column wins over the function argument.
  expect_equal(seen[[1]]$maxIter, 7)
  expect_null(seen[[1]]$junk)
  expect_false(seen[[1]]$verbose)
})

test_that("a single-column grid keeps its column name", {
  X <- as.matrix(iris[, 1:4])
  set.seed(3)
  res <- swarmbHIVE(X, task = "clustering",
                    grid = data.frame(nAntibodies = c(5, 15)),
                    metric = "silhouette", maxIter = 3, verbose = FALSE)
  expect_equal(nrow(res$results), 2L)
  expect_equal(res$results$nAntibodies, c(5, 15))
})

test_that("progress uses message() and handles non-integer grid values", {
  X <- as.matrix(iris[, 1:4])
  grid <- data.frame(nAntibodies = 8, beta = 2.5, epsilon = 0.05)
  set.seed(4)
  msgs <- capture_messages(
    swarmbHIVE(X, y = iris$Species, task = "classification", grid = grid,
               maxIter = 3, verbose = TRUE))
  expect_true(any(grepl("Evaluating combo 1/1: nAntibodies=8, beta=2.5",
                        msgs, fixed = TRUE)))
  expect_true(any(grepl("Best parameters:", msgs, fixed = TRUE)))
  expect_silent(suppressMessages(
    swarmbHIVE(X, y = iris$Species, task = "classification", grid = grid,
               maxIter = 3, verbose = TRUE)))
})

test_that("classification metric helpers match hand calculations", {
  lv <- c("a", "b", "c")
  tb <- table(factor(c("a", "a", "b", "b", "c"), levels = lv),
              factor(c("a", "b", "b", "b", "b"), levels = lv))
  expect_equal(.balanced_accuracy(tb), mean(c(0.5, 1, 0)))
  # F1: a = 2/3, b = 2*(2/4)*1/(2/4+1) = 2/3, c never predicted = 0
  expect_equal(.macro_f1(tb), mean(c(2 / 3, 2 / 3, 0)))
  perfect <- table(factor(lv, levels = lv), factor(lv, levels = lv))
  expect_equal(.cohen_kappa(perfect), 1)
  expect_equal(.swarm_class_metric(c("a", "b"), factor(c("a", "c")),
                                   "accuracy"), 0.5)
})

test_that("cluster metrics are NA for a single cluster", {
  X <- as.matrix(iris[, 1:4])
  expect_true(is.na(.swarm_cluster_metric(rep(1, 150), X, "silhouette",
                                          dist(X))))
  expect_true(is.na(.swarm_cluster_metric(rep(1, 150), X, "davies_bouldin",
                                          NULL)))
})
