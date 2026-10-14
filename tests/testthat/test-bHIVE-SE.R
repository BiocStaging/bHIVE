# Tests for the SummarizedExperiment / SingleCellExperiment methods of bHIVE()

suppressPackageStartupMessages({
  library(SummarizedExperiment)
  library(SingleCellExperiment)
})

X <- as.matrix(iris[, 1:4])
make_se <- function() {
  SummarizedExperiment(assays = list(logcounts = t(X)),
                       colData = DataFrame(species = iris$Species))
}
make_sce <- function() {
  sce <- SingleCellExperiment(assays = list(logcounts = t(X), counts = t(X)),
                              colData = DataFrame(species = iris$Species))
  reducedDim(sce, "PCA") <- prcomp(X)$x[, 1:2]
  sce
}

test_that("bHIVE is an S4 generic that still accepts matrices", {
  expect_true(isGeneric("bHIVE"))
  set.seed(1)
  res <- bHIVE(X, nAntibodies = 15, maxIter = 5, verbose = FALSE)
  expect_type(res, "list")
  expect_equal(length(res$assignments), nrow(X))
  set.seed(1)
  res_df <- bHIVE(as.data.frame(X), nAntibodies = 15, maxIter = 5,
                  verbose = FALSE)
  expect_identical(res_df$assignments, res$assignments)
})

test_that("SummarizedExperiment clustering matches the matrix result", {
  set.seed(2)
  ref <- bHIVE(X, nAntibodies = 15, maxIter = 5, verbose = FALSE)
  set.seed(2)
  se <- bHIVE(make_se(), nAntibodies = 15, maxIter = 5, verbose = FALSE)
  expect_s4_class(se, "SummarizedExperiment")
  expect_s3_class(se$bHIVE, "factor")
  expect_equal(as.numeric(as.character(se$bHIVE)), ref$assignments)
  expect_equal(metadata(se)$bHIVE$antibodies, ref$antibodies)
  expect_s3_class(metadata(se)$bHIVE$model, "AINet")
})

test_that("y can name a colData column for classification", {
  set.seed(3)
  se <- bHIVE(make_se(), y = "species", nAntibodies = 20, maxIter = 5,
              name = "pred", verbose = FALSE)
  expect_identical(levels(se$pred), levels(iris$Species))
  expect_gt(mean(se$pred == se$species), 0.6)
  expect_equal(metadata(se)$pred$task, "classification")
})

test_that("y can be a vector with one label per cell", {
  set.seed(4)
  se <- bHIVE(make_se(), y = as.character(iris$Species), nAntibodies = 20,
              maxIter = 5, verbose = FALSE)
  expect_setequal(levels(se$bHIVE), levels(iris$Species))
})

test_that("SummarizedExperiment input errors are readable", {
  expect_error(bHIVE(make_se(), assay.type = "exprs", verbose = FALSE),
               "assay 'exprs' not found")
  expect_error(bHIVE(make_se(), y = "nope", verbose = FALSE),
               "colData\\(\\) column name")
  se <- make_se()
  se$species[1] <- NA
  expect_error(bHIVE(se, y = "species", verbose = FALSE), "missing labels")
})

test_that("SingleCellExperiment can cluster an assay or a reducedDim", {
  set.seed(5)
  ref <- bHIVE(prcomp(X)$x[, 1:2], nAntibodies = 15, maxIter = 5,
               verbose = FALSE)
  set.seed(5)
  sce <- bHIVE(make_sce(), use.dimred = "PCA", nAntibodies = 15,
               maxIter = 5, verbose = FALSE)
  expect_s4_class(sce, "SingleCellExperiment")
  expect_equal(as.numeric(as.character(sce$bHIVE)), ref$assignments)

  set.seed(6)
  a <- bHIVE(make_sce(), assay.type = "counts", nAntibodies = 15,
             maxIter = 5, verbose = FALSE)
  set.seed(6)
  b <- bHIVE(X, nAntibodies = 15, maxIter = 5, verbose = FALSE)
  expect_equal(as.numeric(as.character(a$bHIVE)), b$assignments)

  expect_error(bHIVE(make_sce(), use.dimred = "UMAP", verbose = FALSE),
               "reducedDim 'UMAP' not found")
})

test_that("the stored model predicts new cells", {
  set.seed(7)
  sce <- bHIVE(make_sce(), y = "species", nAntibodies = 20, maxIter = 5,
               verbose = FALSE)
  pred <- metadata(sce)$bHIVE$model$predict(X[1:10, ])
  expect_identical(pred, as.character(sce$bHIVE[1:10]))
})

test_that("swarmbHIVE and honeycombHIVE point SE users to the right call", {
  expect_error(honeycombHIVE(make_se(), layers = 1, verbose = FALSE),
               "is a SummarizedExperiment")
  expect_error(swarmbHIVE(make_se(), grid = data.frame(nAntibodies = 5),
                          verbose = FALSE),
               "is a SummarizedExperiment")
})

test_that("tunable arguments cover bHIVE and AINet-only settings", {
  args <- .bhive_tunable_args()
  expect_true(all(c("nAntibodies", "epsilon", "targetK", "scale",
                    "coverageBoost", "consolidate", "scaleCofactor") %in% args))
  expect_false(any(c("X", "y", "task", "verbose", "...") %in% args))
})
