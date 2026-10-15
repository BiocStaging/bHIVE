# Golden-output harness for AINet$fit().
#
# Each config is a closure that builds fresh data, modules and model under a
# fixed seed, fits, and returns the fitted model. run_golden() reduces a model
# to the fields that must stay bit-identical across a pure refactor: the
# result, the repertoire, the per-iteration history, any memory-pool state,
# and .Random.seed after fit (catches RNG call-order drift even when the
# outputs happen to coincide).
#
# The fixture lives in tests/testthat/fixtures/ainet-golden.rds and is written
# by tests/testthat/fixtures/make-ainet-golden.R. Floating-point results depend
# on the BLAS, so the comparison only runs locally (skip_on_cran/skip_on_ci).

golden_configs <- function() {
  X_iris <- as.matrix(iris[, 1:4])
  y_iris <- iris$Species

  two_blobs <- function() {
    set.seed(7)
    rbind(matrix(rnorm(1800, 0), ncol = 3),
          matrix(rnorm(1800, 4), ncol = 3))
  }

  list(
    clust_default = function() {
      set.seed(101)
      m <- AINet$new(nAntibodies = 20, maxIter = 10, verbose = FALSE)
      m$fit(X_iris, task = "clustering")
      list(model = m)
    },
    clust_targetK_zscore = function() {
      set.seed(102)
      m <- AINet$new(nAntibodies = 20, maxIter = 10, targetK = 3,
                     scale = "zscore", verbose = FALSE)
      m$fit(X_iris, task = "clustering")
      list(model = m)
    },
    clust_auto_alpha_robust = function() {
      X <- two_blobs()
      set.seed(103)
      m <- AINet$new(nAntibodies = 15, maxIter = 5, scale = "robust",
                     affinityParams = list(alpha = "auto"),
                     epsilonQuantile = 0.1, verbose = FALSE)
      m$fit(X, task = "clustering")
      list(model = m)
    },
    class_default = function() {
      set.seed(104)
      m <- AINet$new(nAntibodies = 20, maxIter = 10, verbose = FALSE)
      m$fit(X_iris, y_iris, task = "classification")
      list(model = m)
    },
    class_all_modules = function() {
      set.seed(105)
      m <- AINet$new(
        nAntibodies = 40, maxIter = 6, verbose = FALSE,
        shm = SHMEngine$new(method = "adaptive"),
        init = VDJLibrary$new(),
        activation = ActivationGate$new(),
        idiotypic = IdiotypicNetwork$new(theta_low = 1e-4, theta_high = 0.999,
                                         survival_threshold = 0.01),
        germinalCenter = GerminalCenter$new(nTfh = 20),
        microenvironment = Microenvironment$new(),
        classSwitcher = ClassSwitcher$new(),
        memory = MemoryPool$new(archive_threshold = 0.01))
      m$fit(X_iris, y_iris, task = "classification")
      list(model = m, memory = m$modules$memory$memory_cells)
    },
    clust_all_modules = function() {
      set.seed(106)
      m <- AINet$new(
        nAntibodies = 40, maxIter = 6, verbose = FALSE,
        shm = SHMEngine$new(method = "hotspot"),
        activation = ActivationGate$new(),
        idiotypic = IdiotypicNetwork$new(theta_low = 1e-4, theta_high = 0.999,
                                         survival_threshold = 0.01),
        germinalCenter = GerminalCenter$new(nTfh = 20),
        microenvironment = Microenvironment$new(),
        classSwitcher = ClassSwitcher$new())
      m$fit(X_iris, task = "clustering")
      list(model = m)
    },
    clust_coverage_kmeanspp = function() {
      set.seed(107)
      m <- AINet$new(nAntibodies = 10, maxIter = 8, initMethod = "kmeans++",
                     coverageBoost = TRUE, verbose = FALSE)
      m$fit(X_iris, task = "clustering")
      list(model = m)
    },
    clust_memory_recall = function() {
      set.seed(108)
      mem <- MemoryPool$new(archive_threshold = 0.01)
      m1 <- AINet$new(nAntibodies = 15, maxIter = 5, memory = mem,
                      verbose = FALSE)
      m1$fit(X_iris, task = "clustering")
      m2 <- AINet$new(nAntibodies = 15, maxIter = 5, memory = mem,
                      verbose = FALSE)
      m2$fit(X_iris, task = "clustering")
      list(model = m2, first = run_golden_model(m1),
           memory = mem$memory_cells)
    },
    clust_cosine = function() {
      set.seed(109)
      m <- AINet$new(nAntibodies = 15, maxIter = 8, affinityFunc = "cosine",
                     distFunc = "cosine", epsilon = 5e-5, verbose = FALSE)
      m$fit(X_iris, task = "clustering")
      list(model = m)
    },
    clust_laplace_manhattan_early_stop = function() {
      set.seed(110)
      m <- AINet$new(nAntibodies = 15, maxIter = 30, affinityFunc = "laplace",
                     distFunc = "manhattan", consolidate = FALSE,
                     noImprovementLimit = 3, verbose = FALSE)
      m$fit(X_iris, task = "clustering")
      list(model = m)
    }
  )
}

run_golden_model <- function(m) {
  list(result     = m$result,
       history    = m$history,
       cells      = m$repertoire$cells,
       metadata   = m$repertoire$metadata)
}

run_golden <- function(cfg) {
  out <- cfg()
  res <- run_golden_model(out$model)
  res$seed <- .Random.seed
  res$extra <- out[setdiff(names(out), "model")]
  res
}
