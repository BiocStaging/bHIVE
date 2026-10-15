# Tests for the caret wrappers bHIVEmodel and honeycombHIVEmodel

X <- as.matrix(iris[, 1:4])
y <- iris$Species
# 80 antibodies, not 40. bHIVE is stochastic and only its own RNG stream is
# reproducible: the number of draws per iteration depends on how many
# antibodies survive the C++ epsilon-suppression step, so a one-ULP difference
# in the distance kernels on another architecture moves the whole trajectory.
# With 40 antibodies too few prototypes survive to split versicolor from
# virginica on an unlucky trajectory and accuracy falls to 0.68 (seen on the
# Bioconductor linux-arm64 builder). Over 200 seeds x beta in {3, 5}, 80
# antibodies stayed at or above 0.85 while 40 fell below 0.80 on 12% of runs.
one_row <- data.frame(nAntibodies = 80, beta = 5, epsilon = 0.05)
# caret signals a failed resample with these warnings. Other warnings (for
# example from attaching caret's dependencies) still surface as usual.
fold_failure <- "fit failed|missing values in resampled"

test_that("bHIVEmodel predicts with the fitted model (reviewer example)", {
  # The reviewer's snippet: direct fit/predict calls, no caret needed.
  set.seed(42)
  fit <- bHIVEmodel$fit(x = X, y = y, wts = NULL, lev = levels(y),
                        last = TRUE, classProbs = FALSE, param = one_row,
                        verbose = FALSE)
  pred <- bHIVEmodel$predict(fit, X)
  expect_s3_class(pred, "factor")
  expect_identical(levels(pred), levels(y))
  expect_identical(as.character(pred), fit$assignments)
  expect_gt(mean(pred == y), 0.75)
  # Every class is predicted, not just the first one.
  expect_setequal(unique(as.character(pred)), levels(y))
})

test_that("bHIVEmodel declares classification only", {
  expect_identical(bHIVEmodel$type, "Classification")
  expect_null(bHIVEmodel$prob)
  expect_error(bHIVEmodel$fit(x = X, y = as.numeric(y), param = one_row),
               "classification only")
})

test_that("bHIVEmodel grids support grid and random search", {
  g <- bHIVEmodel$grid(X, y, len = 2)
  expect_equal(nrow(g), 8L)
  expect_true(all(g$nAntibodies == round(g$nAntibodies)))
  set.seed(1)
  r <- bHIVEmodel$grid(X, y, len = 4, search = "random")
  expect_equal(nrow(r), 4L)
})

test_that("bHIVEmodel trains and predicts through caret", {
  skip_if_not_installed("caret")
  set.seed(42)
  expect_no_warning(
    fit <- caret::train(x = X, y = y, method = bHIVEmodel,
                        trControl = caret::trainControl(method = "cv",
                                                        number = 3),
                        tuneGrid = rbind(one_row, transform(one_row, beta = 3)),
                        maxIter = 10),
    message = fold_failure)
  expect_s3_class(fit, "train")
  expect_false(anyNA(fit$results$Accuracy))
  expect_gt(mean(predict(fit, X) == y), 0.75)
})

test_that("bHIVEmodel works with tuneLength and warns on classProbs", {
  skip_if_not_installed("caret")
  set.seed(43)
  fit <- caret::train(x = X, y = y, method = bHIVEmodel, tuneLength = 2,
                      trControl = caret::trainControl(method = "cv",
                                                      number = 2),
                      maxIter = 3)
  expect_equal(nrow(fit$results), 8L)
  expect_warning(
    caret::train(x = X, y = y, method = bHIVEmodel, tuneGrid = one_row,
                 trControl = caret::trainControl(method = "none",
                                                 classProbs = TRUE),
                 maxIter = 3),
    "do(es)? not implement")
})

test_that("honeycombHIVEmodel trains on every fold and predicts", {
  skip_if_not_installed("caret")
  grid <- data.frame(nAntibodies = 20, beta = 5, epsilon = 0.05, layers = 2,
                     refineOptimizer = "adam", refineSteps = 3,
                     refineLR = 0.01)
  set.seed(44)
  # expand.grid-style factor columns used to fail every fold in match.arg().
  grid$refineOptimizer <- factor(grid$refineOptimizer)
  expect_no_warning(
    fit <- caret::train(x = X, y = y, method = honeycombHIVEmodel,
                        trControl = caret::trainControl(method = "cv",
                                                        number = 3),
                        tuneGrid = grid, maxIter = 10),
    message = fold_failure)
  expect_false(anyNA(fit$results$Accuracy))
  pred <- predict(fit, X)
  expect_identical(levels(pred), levels(y))
  expect_gt(mean(pred == y), 0.6)
})

test_that("honeycombHIVEmodel declares classification only", {
  expect_identical(honeycombHIVEmodel$type, "Classification")
  expect_null(honeycombHIVEmodel$prob)
  g <- honeycombHIVEmodel$grid(X, y, len = 2)
  expect_true(all(c("refineOptimizer", "refineSteps", "refineLR") %in%
                    names(g)))
})
