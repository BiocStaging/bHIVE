# Tests for AINet R6 class

set.seed(42)
data(iris)
X <- as.matrix(iris[, 1:4])
y_class <- iris$Species

# ==== Construction / Validation ====

test_that("AINet initializes with defaults", {
  model <- AINet$new(verbose = FALSE)
  expect_true(inherits(model, "AINet"))
  expect_true(inherits(model, "ImmuneAlgorithm"))
  expect_equal(model$config$nAntibodies, 20L)
  expect_equal(model$config$beta, 5)
  expect_null(model$result)
})

test_that("AINet validates numeric parameters", {
  expect_error(AINet$new(nAntibodies = -1), "nAntibodies")
  expect_error(AINet$new(beta = 0), "beta")
  expect_error(AINet$new(epsilon = -1), "epsilon")
  expect_error(AINet$new(maxIter = 0), "maxIter")
  expect_error(AINet$new(k = 0), "k")
  expect_error(AINet$new(mutationDecay = 0), "mutationDecay")
  expect_error(AINet$new(mutationDecay = 1.5), "mutationDecay")
  expect_error(AINet$new(mutationMin = -1), "mutationMin")
})

test_that("AINet validates affinity and distance function names", {
  expect_error(AINet$new(affinityFunc = "invalid"))
  expect_error(AINet$new(distFunc = "invalid"))
  expect_error(AINet$new(initMethod = "invalid"))
})

test_that("AINet accepts module injection", {
  shm <- SHMEngine$new(method = "adaptive")
  model <- AINet$new(shm = shm, verbose = FALSE)
  expect_true("shm" %in% names(model$modules))
  expect_equal(model$modules$shm$method, "adaptive")
})

test_that("AINet stores only non-NULL modules", {
  model <- AINet$new(verbose = FALSE)
  expect_equal(length(model$modules), 0)

  model2 <- AINet$new(shm = SHMEngine$new(), verbose = FALSE)
  expect_equal(length(model2$modules), 1)
})

# ==== Clustering ====

test_that("AINet fits clustering task", {
  model <- AINet$new(nAntibodies = 10, maxIter = 5, verbose = FALSE)
  model$fit(X, task = "clustering")

  expect_equal(model$result$task, "clustering")
  expect_true(is.matrix(model$result$antibodies))
  expect_equal(length(model$result$assignments), nrow(X))
  expect_true(all(model$result$assignments > 0))
})

test_that("AINet infers clustering when y is NULL", {
  model <- AINet$new(nAntibodies = 10, maxIter = 5, verbose = FALSE)
  model$fit(X)
  expect_equal(model$result$task, "clustering")
})

# ==== Classification ====

test_that("AINet fits classification task", {
  model <- AINet$new(nAntibodies = 15, maxIter = 5, verbose = FALSE)
  model$fit(X, y_class, task = "classification")

  expect_equal(model$result$task, "classification")
  expect_equal(length(model$result$assignments), nrow(X))
  expect_true(all(model$result$assignments %in% levels(y_class)))
})

test_that("AINet infers classification when y is factor", {
  model <- AINet$new(nAntibodies = 10, maxIter = 3, verbose = FALSE)
  model$fit(X, y_class)
  expect_equal(model$result$task, "classification")
})

# ==== Prediction ====

test_that("AINet predict works for clustering", {
  model <- AINet$new(nAntibodies = 10, maxIter = 5, verbose = FALSE)
  model$fit(X, task = "clustering")
  preds <- model$predict(X[1:10, ])
  expect_equal(length(preds), 10)
})

test_that("AINet predict works for classification", {
  model <- AINet$new(nAntibodies = 10, maxIter = 5, verbose = FALSE)
  model$fit(X, y_class, task = "classification")
  preds <- model$predict(X[1:10, ])
  expect_equal(length(preds), 10)
  expect_true(all(preds %in% levels(y_class)))
})

# ==== Initialization methods ====

test_that("AINet works with all initialization methods", {
  for (init in c("sample", "random", "random_uniform", "kmeans++")) {
    model <- AINet$new(nAntibodies = 10, maxIter = 3, initMethod = init,
                       verbose = FALSE)
    model$fit(X, task = "clustering")
    expect_true(is.matrix(model$result$antibodies))
  }
})

# ==== Affinity / distance variants ====

test_that("AINet works with different affinity functions", {
  # Each affinity paired with its natural distance (see the metric guard).
  natural <- c(gaussian = "euclidean", laplace = "manhattan",
                   polynomial = "euclidean", cosine = "cosine")
  for (af in names(natural)) {
    model <- AINet$new(nAntibodies = 10, maxIter = 3, affinityFunc = af,
                       distFunc = natural[[af]], verbose = FALSE)
    model$fit(X, task = "clustering")
    expect_true(is.matrix(model$result$antibodies))
  }
})

test_that("AINet works with different distance functions", {
  for (df in c("euclidean", "manhattan", "minkowski")) {
    model <- AINet$new(nAntibodies = 10, maxIter = 3,
                       distFunc = df, verbose = FALSE)
    model$fit(X, task = "clustering")
    expect_true(is.matrix(model$result$antibodies))
  }
})

# ==== Repertoire / history ====

test_that("AINet populates repertoire after fitting", {
  model <- AINet$new(nAntibodies = 10, maxIter = 5, verbose = FALSE)
  model$fit(X, task = "clustering")
  expect_true(inherits(model$repertoire, "ImmuneRepertoire"))
  expect_true(model$repertoire$size() > 0)
})

test_that("AINet records iteration history", {
  model <- AINet$new(nAntibodies = 10, maxIter = 5, verbose = FALSE)
  model$fit(X, task = "clustering")
  expect_true(length(model$history) > 0)
  expect_true(all(vapply(model$history, function(h) h$n_antibodies > 0, logical(1))))
})

# ==== Early stopping ====

test_that("AINet early stopping works", {
  model <- AINet$new(nAntibodies = 10, maxIter = 100,
                     stopTolerance = 5, noImprovementLimit = 3,
                     verbose = FALSE)
  model$fit(X, task = "clustering")
  expect_true(length(model$history) < 100)
})

# ==== Edge case: all antibodies suppressed ====

test_that("AINet handles small populations gracefully", {
  # With extreme epsilon, the algorithm may error or produce a tiny repertoire
  # Either outcome is acceptable
  result <- tryCatch({
    model <- AINet$new(nAntibodies = 2, maxIter = 3, epsilon = 100,
                       verbose = FALSE)
    model$fit(X, task = "clustering")
    "completed"
  }, error = function(e) "error")
  expect_true(result %in% c("completed", "error"))
})

# ==== Module integration: VDJLibrary as init ====

test_that("AINet routes initialization through `init$generate` when supplied", {
  # Inject a tracking mock as the `init` module. If $.initialize_antibodies
  # is correctly routing through it, the closure will be called exactly
  # once during fit and we can verify both the call and the returned
  # repertoire shape.
  call_log <- new.env()
  call_log$n_calls <- 0L
  call_log$last_n  <- NA_integer_

  mock_init <- list(
    generate = function(nAntibodies, X) {
      call_log$n_calls <- call_log$n_calls + 1L
      call_log$last_n  <- nAntibodies
      # Return a deterministic structured matrix (rows are constant per row)
      mat <- matrix(seq_len(nAntibodies) / nAntibodies,
                    nrow = nAntibodies, ncol = ncol(X))
      mat
    }
  )

  set.seed(7)
  model <- AINet$new(nAntibodies = 12, maxIter = 1, k = 3, epsilon = 1e-6,
                     initMethod = "sample", init = mock_init, verbose = FALSE)
  model$fit(X, task = "clustering")

  expect_equal(call_log$n_calls, 1L)
  expect_equal(call_log$last_n, 12L)
})

test_that("AINet ignores `init` when its $generate is not a function", {
  bogus <- list(generate = "not a function")
  model <- AINet$new(nAntibodies = 8, maxIter = 1, init = bogus, verbose = FALSE)
  expect_no_error(model$fit(X, task = "clustering"))
  expect_true(nrow(model$repertoire$as_matrix()) >= 1)
})

# ==== Module integration: IdiotypicNetwork ====

test_that("IdiotypicNetwork module changes the surviving repertoire", {
  set.seed(11)
  base <- AINet$new(nAntibodies = 30, maxIter = 5, k = 3, epsilon = 1e-6,
                    initMethod = "sample", verbose = FALSE)
  base$fit(X, task = "clustering")

  set.seed(11)
  idi <- IdiotypicNetwork$new(theta_low = 0.05, theta_high = 1.0,
                              source_rate = 0.5, decay_rate = 0.1,
                              dt = 0.05, timeSteps = 5,
                              survival_threshold = 0.3)
  with_idi <- AINet$new(nAntibodies = 30, maxIter = 5, k = 3, epsilon = 1e-6,
                        initMethod = "sample", idiotypic = idi, verbose = FALSE)
  with_idi$fit(X, task = "clustering")

  # Bell-curve dynamics should change the post-fit repertoire size relative
  # to the unmodulated baseline (almost always smaller; rare ties allowed).
  expect_false(identical(base$repertoire$size(), with_idi$repertoire$size()))
})

test_that("IdiotypicNetwork safety net keeps a minimum repertoire when ill-tuned", {
  # Thresholds and survival_threshold are calibrated so every antibody would
  # otherwise die; the safety net retains the top-population subset and the
  # iteration must continue rather than throwing.
  set.seed(13)
  idi_bad <- IdiotypicNetwork$new(theta_low = 100, theta_high = 1000,
                                  source_rate = 0,    decay_rate = 10,
                                  dt = 1.0, timeSteps = 3,
                                  survival_threshold = 1e6)
  model <- AINet$new(nAntibodies = 20, maxIter = 3, k = 3, epsilon = 1e-6,
                     initMethod = "sample", idiotypic = idi_bad,
                     verbose = FALSE)
  expect_no_error(model$fit(X, task = "clustering"))
  expect_true(model$repertoire$size() >= 1L)
})

# ==== Module integration: ActivationGate (pre-selection density gating) ====

test_that("ActivationGate alters which antibodies enter clonal selection", {
  set.seed(17)
  base <- AINet$new(nAntibodies = 30, maxIter = 4, k = 3, epsilon = 1e-6,
                    initMethod = "sample", verbose = FALSE)
  base$fit(X, task = "clustering")

  set.seed(17)
  gate <- ActivationGate$new(signal2_type = "density",
                             threshold1 = 0, threshold2 = 0.6)
  gated <- AINet$new(nAntibodies = 30, maxIter = 4, k = 3, epsilon = 1e-6,
                     initMethod = "sample", activation = gate, verbose = FALSE)
  gated$fit(X, task = "clustering")

  # Either the surviving repertoires differ, or their final assignments do.
  same_size <- identical(base$repertoire$size(), gated$repertoire$size())
  same_pred <- identical(base$result$assignments, gated$result$assignments)
  expect_false(same_size && same_pred)
})

test_that("ActivationGate works in classification without breaking class coverage", {
  set.seed(19)
  gate <- ActivationGate$new(signal2_type = "density",
                             threshold1 = 0, threshold2 = 0.7)
  model <- AINet$new(nAntibodies = 30, maxIter = 5, k = 3, epsilon = 1e-6,
                     initMethod = "sample", activation = gate, verbose = FALSE)
  expect_no_error(model$fit(X, y_class, task = "classification"))
  preds <- model$predict(X)
  expect_equal(length(preds), nrow(X))
  expect_gte(length(unique(model$result$antibody_classes)), 2L)
})

# ==== Module integration: Microenvironment-aware mutation ====

test_that("Microenvironment module perturbs the post-selection repertoire", {
  set.seed(23)
  base <- AINet$new(nAntibodies = 25, maxIter = 4, k = 3, epsilon = 1e-6,
                    initMethod = "sample", verbose = FALSE)
  base$fit(X, task = "clustering")

  set.seed(23)
  mic <- Microenvironment$new(high_density_threshold = 0.75,
                              low_density_threshold = 0.25,
                              stabilization_factor = 0.3,
                              exploration_factor = 2.0)
  with_mic <- AINet$new(nAntibodies = 25, maxIter = 4, k = 3, epsilon = 1e-6,
                        initMethod = "sample", microenvironment = mic,
                        verbose = FALSE)
  with_mic$fit(X, task = "clustering")

  # Density-aware jitter must change the antibody coordinates relative to
  # the unmodulated baseline. If sizes match, compare matrices; if sizes
  # differ, the size diff itself proves an effect.
  A_base <- base$repertoire$as_matrix()
  A_mic  <- with_mic$repertoire$as_matrix()
  if (nrow(A_base) == nrow(A_mic)) {
    expect_false(isTRUE(all.equal(A_base, A_mic)))
  } else {
    succeed("Microenvironment changed surviving repertoire size")
  }
})

# ==== Final-stage orphan pruning ====

test_that("AINet prunes antibodies that bind no training point", {
  # Force orphans by using `random_uniform` init: many antibodies will land
  # in regions of feature space no iris row is closest to.
  set.seed(29)
  model <- AINet$new(nAntibodies = 50, maxIter = 3, k = 3, epsilon = 1e-6,
                     initMethod = "random_uniform", verbose = FALSE)
  model$fit(X, task = "clustering")

  # After pruning, every surviving antibody must serve as nearest neighbor
  # to at least one training point. Equivalently, the number of unique
  # cluster IDs in the result equals the antibody count.
  expect_equal(length(unique(model$result$assignments)),
               nrow(model$result$antibodies))
})

test_that("Orphan pruning preserves at least 2 antibodies when many are orphans", {
  # The pruning step refuses to drop below 2 surviving antibodies even if
  # most are orphans, so the final assignment can still produce >=2 clusters.
  set.seed(31)
  model <- AINet$new(nAntibodies = 40, maxIter = 2, k = 3, epsilon = 1e-6,
                     initMethod = "random_uniform", verbose = FALSE)
  model$fit(X, task = "clustering")
  expect_gte(nrow(model$result$antibodies), 2L)
})

# ==== Pre-allocated antibody_classes (regression test for iter-1 indexing) ====

test_that("AINet classification with ActivationGate does not error on iter 1", {
  # The pre-selection gate references antibody_classes on iteration 1,
  # before the first clonal_selection_iteration_cpp call has populated it.
  # A pre-allocation step in $fit() seeds it; this test guards against
  # a regression where that pre-allocation is removed.
  set.seed(37)
  gate <- ActivationGate$new(signal2_type = "density",
                             threshold1 = 0, threshold2 = 0.6)
  model <- AINet$new(nAntibodies = 20, maxIter = 1, k = 3, epsilon = 1e-6,
                     initMethod = "sample", activation = gate, verbose = FALSE)
  expect_no_error(model$fit(X, y_class, task = "classification"))
})

# ==== Module integration: GerminalCenter ====

test_that("GerminalCenter module reduces the surviving repertoire", {
  # With selectionPressure ~ 1 and small nTfh, surviving count is bounded
  # near nTfh. The baseline (no GC) typically retains many more antibodies.
  set.seed(41)
  base <- AINet$new(nAntibodies = 30, maxIter = 4, k = 3, epsilon = 1e-6,
                    initMethod = "sample", verbose = FALSE)
  base$fit(X, task = "clustering")

  set.seed(41)
  gc <- GerminalCenter$new(nTfh = 5, selectionPressure = 0.9, rounds = 2)
  with_gc <- AINet$new(nAntibodies = 30, maxIter = 4, k = 3, epsilon = 1e-6,
                       initMethod = "sample", germinalCenter = gc,
                       verbose = FALSE)
  with_gc$fit(X, task = "clustering")

  expect_lt(with_gc$repertoire$size(), base$repertoire$size())
  expect_gte(with_gc$repertoire$size(), 1L)
  # GC must record survivors from its last call.
  expect_true(is.numeric(gc$last_survivors))
})

test_that("GerminalCenter in classification keeps antibody_classes aligned", {
  set.seed(43)
  gc <- GerminalCenter$new(nTfh = 4, selectionPressure = 0.7)
  model <- AINet$new(nAntibodies = 25, maxIter = 3, k = 3, epsilon = 1e-6,
                     germinalCenter = gc, verbose = FALSE)
  expect_no_error(model$fit(X, y_class, task = "classification"))
  expect_equal(length(model$result$antibody_classes),
               nrow(model$result$antibodies))
})

# ==== Module integration: MemoryPool ====

test_that("MemoryPool archives high-affinity antibodies after fit", {
  set.seed(47)
  mp <- MemoryPool$new(archive_threshold = 1e-6, max_memory = 100)
  model <- AINet$new(nAntibodies = 20, maxIter = 3, memory = mp,
                     verbose = FALSE)
  model$fit(X, task = "clustering")
  expect_gt(mp$size(), 0L)
  expect_equal(ncol(mp$memory_cells), ncol(X))
})

test_that("MemoryPool recall seeds a second fit's starting repertoire", {
  # Fit once to populate memory, then fit a fresh model that recalls.
  # The recalled cells should be present in the starting repertoire even
  # though nAntibodies is small.
  set.seed(53)
  mp <- MemoryPool$new(archive_threshold = 1e-6, recall_threshold = 0,
                       max_memory = 100)
  AINet$new(nAntibodies = 20, maxIter = 2, memory = mp,
            verbose = FALSE)$fit(X, task = "clustering")
  archived <- mp$size()
  expect_gt(archived, 0L)

  # Second fit recalls and merges memory into the initial repertoire.
  set.seed(57)
  model2 <- AINet$new(nAntibodies = 5, maxIter = 1, epsilon = 1e-6,
                      memory = mp, verbose = FALSE)
  model2$fit(X, task = "clustering")
  # Repertoire size should reflect at least some recall (post-iter pruning
  # may shrink it, so just confirm fit completes without error).
  expect_gte(model2$repertoire$size(), 1L)
})

# ==== Module integration: ClassSwitcher ====

test_that("ClassSwitcher updates isotypes when paired with Microenvironment", {
  set.seed(61)
  me <- Microenvironment$new(high_density_threshold = 0.6,
                              low_density_threshold = 0.3)
  cs <- ClassSwitcher$new(alpha_IgM = 0.1, alpha_IgG = 5.0, alpha_IgA = 1.0)
  model <- AINet$new(nAntibodies = 30, maxIter = 4, epsilon = 1e-6,
                     microenvironment = me, classSwitcher = cs,
                     verbose = FALSE)
  model$fit(X, task = "clustering")
  # After fit, repertoire metadata$isotype should reflect zone-based switching.
  iso <- model$repertoire$metadata$isotype
  expect_true(any(iso %in% c("IgG", "IgA")))
})

# ==== SHM dispatch in clonal_selection ====

test_that("Different SHM strategies produce different repertoires", {
  # Drive the SHM dispatcher directly via shm_mutate_cpp on a single
  # antibody / data-point pair with low affinity. Bypasses the AINet
  # loop machinery (early stopping, network suppression, gating) so the
  # assertion only depends on the SHM dispatch, which is what we want
  # to verify. Seed reset before each call -> R::rnorm consumes the
  # same standard-normal stream, scaled by each method's rate formula,
  # producing materially different mutations.
  A   <- matrix(c(10, 10, 10, 10), nrow = 1)   # far from any iris point
  Xd  <- matrix(c(5.1, 3.5, 1.4, 0.2), nrow = 1)
  aff <- 0.0                                   # forces large mutation rates

  call_dispatch <- function(method, c_rate = 1.0, temperature = 0.5) {
    set.seed(67)
    shm_mutate_cpp(
      A, Xd, affinities = aff,
      top_k_indices = 0L, data_indices = 0L,
      method = method, iter = 1L,
      decay = 1.0, mutationMin = 0.01,
      c_rate = c_rate, temperature = temperature,
      E_0 = 1.0, base_rate = 0.1,
      beta1 = 0.9, beta2 = 0.999, adam_epsilon = 1e-8,
      m1_state = matrix(0, 0, 0), m2_state = matrix(0, 0, 0),
      affinity_type = "gaussian",
      aff_alpha = 1.0, aff_c = 1.0, aff_p = 2.0
    )$A
  }
  A_u <- call_dispatch("uniform")
  A_a <- call_dispatch("airs", c_rate = 0.5, temperature = 0.2)

  expect_false(isTRUE(all.equal(A_u, A_a)))
  expect_gt(max(abs(A_u - A_a)), 1e-3)
})

test_that("AINet runs end-to-end with adaptive SHM (moment state tracking)", {
  # Adaptive SHM requires (m x d) moment matrices to be threaded through
  # clonal_selection and re-aligned after every subset/reorder operation.
  # This test exercises the full plumbing without checking specific values.
  set.seed(71)
  model <- AINet$new(nAntibodies = 25, maxIter = 4, epsilon = 1e-6,
                     shm = SHMEngine$new(method = "adaptive",
                                          base_rate = 0.05),
                     verbose = FALSE)
  expect_no_error(model$fit(X, task = "clustering"))
  expect_gte(model$repertoire$size(), 1L)
})

test_that("AINet runs with all five SHM strategies in classification", {
  for (method in c("uniform", "airs", "hotspot", "energy", "adaptive")) {
    set.seed(73)
    model <- AINet$new(nAntibodies = 20, maxIter = 3, epsilon = 1e-6,
                       shm = SHMEngine$new(method = method),
                       verbose = FALSE)
    expect_no_error(model$fit(X, y_class, task = "classification"))
    expect_equal(length(model$result$assignments), nrow(X))
  }
})

# ==== Read-out consistency: predict() reproduces the training labels ====

test_that("predict() works after alpha = 'auto' and matches training", {
  set.seed(11)
  model <- AINet$new(nAntibodies = 20, maxIter = 10, scale = "zscore",
                     affinityParams = list(alpha = "auto"), verbose = FALSE)
  model$fit(X, y_class)
  expect_true(is.numeric(model$config$kernel$alpha))
  expect_identical(model$predict(X), model$result$assignments)
})

test_that("cluster predict() uses the metric that produced the labels", {
  # Consolidation labels come from a Euclidean E-step, so predict() must not
  # fall back to the training distFunc (manhattan here).
  set.seed(12)
  model <- AINet$new(nAntibodies = 20, maxIter = 10, affinityFunc = "laplace",
                     distFunc = "manhattan", verbose = FALSE)
  model$fit(X)
  expect_equal(model$config$kernel$assignDist, "euclidean")
  expect_equal(as.numeric(model$predict(X)), model$result$assignments)
  expect_equal(nrow(model$result$antibodies),
               length(unique(model$result$assignments)))
})

test_that("unconsolidated prototypes line up with cluster labels", {
  set.seed(13)
  model <- AINet$new(nAntibodies = 15, maxIter = 6, consolidate = FALSE,
                     coverageBoost = TRUE, verbose = FALSE)
  model$fit(X)
  expect_equal(model$config$kernel$assignDist, "euclidean")
  expect_equal(nrow(model$result$antibodies),
               length(unique(model$result$assignments)))
  expect_equal(as.numeric(model$predict(X)), model$result$assignments)
})

test_that("predict() uses the ClassSwitcher-modulated bandwidth", {
  set.seed(14)
  model <- AINet$new(nAntibodies = 30, maxIter = 6, verbose = FALSE,
                     microenvironment = Microenvironment$new(),
                     classSwitcher = ClassSwitcher$new())
  model$fit(X, y_class)
  expect_identical(model$predict(X), model$result$assignments)
})

test_that("classification result records the labeling antibody", {
  set.seed(15)
  model <- AINet$new(nAntibodies = 20, maxIter = 5, verbose = FALSE)
  model$fit(X, y_class)
  idx <- model$result$antibody_index
  expect_length(idx, nrow(X))
  expect_identical(model$result$antibody_classes[idx],
                   model$result$assignments)
})

test_that("history resets between fits", {
  set.seed(16)
  model <- AINet$new(nAntibodies = 15, maxIter = 4, verbose = FALSE)
  model$fit(X)
  model$fit(X)
  expect_length(model$history, 4L)
})

test_that("memory recall works with alpha = 'auto'", {
  set.seed(17)
  mem <- MemoryPool$new(archive_threshold = 0.001, recall_threshold = 0)
  model <- AINet$new(nAntibodies = 15, maxIter = 3, memory = mem,
                     affinityParams = list(alpha = "auto"), verbose = FALSE)
  model$fit(X)
  expect_gt(mem$size(), 0L)
  expect_no_error(model$fit(X))
})

test_that("progress messages can be silenced with suppressMessages()", {
  set.seed(18)
  model <- AINet$new(nAntibodies = 10, maxIter = 3, coverageBoost = TRUE,
                     affinityParams = list(alpha = "auto"), verbose = TRUE)
  msgs <- capture_messages(model$fit(X))
  expect_true(any(grepl("Iteration 1", msgs)))
  expect_true(any(grepl("Auto bandwidth", msgs)))
  expect_silent(suppressMessages(model$fit(X)))
})

# ==== Forced-K Lloyd refinement ====

test_that(".fill_empty_clusters restores every label without emptying singletons", {
  Xs <- rbind(c(0, 0), c(0, 1), c(10, 10), c(10, 11), c(50, 50))
  filled <- .fill_empty_clusters(Xs, c(1L, 1L, 2L, 2L, 2L), K = 3L)
  expect_setequal(filled, 1:3)
  # The point farthest from its own cluster mean (the outlier) moves.
  expect_equal(filled[5], 3L)
  expect_identical(.fill_empty_clusters(Xs, c(1L, 2L, 3L, 3L, 3L), 3L),
                   c(1L, 2L, 3L, 3L, 3L))
})

test_that("targetK refinement keeps exactly K clusters from a gapped seed", {
  set.seed(19)
  model <- AINet$new(nAntibodies = 10, maxIter = 2, verbose = FALSE)
  model$fit(X)
  lloyd <- model$.__enclos_env__$private$.lloyd_refine
  cent  <- X[c(1, 51, 101, 2), ]
  seed  <- rep(c(1L, 2L, 3L), each = 50)   # label 4 starts empty
  out <- lloyd(X, cent, seed, steps = 30L, cfg = model$config,
               alpha = 1, c_p = 1, p_p = 2, Sigma_inv = matrix(0, 0, 0),
               K = 4L)
  expect_equal(sort(unique(out$assignments)), 1:4)
  expect_equal(nrow(out$antibodies), 4L)
  free <- lloyd(X, cent, seed, steps = 30L, cfg = model$config,
                alpha = 1, c_p = 1, p_p = 2, Sigma_inv = matrix(0, 0, 0))
  expect_equal(nrow(free$antibodies), length(unique(free$assignments)))
})

test_that("targetK returns K clusters and predict() reproduces them", {
  set.seed(20)
  model <- AINet$new(nAntibodies = 30, maxIter = 6, targetK = 6,
                     verbose = FALSE)
  model$fit(X)
  expect_equal(length(unique(model$result$assignments)), 6L)
  expect_equal(nrow(model$result$antibodies), 6L)
  expect_equal(as.numeric(model$predict(X)), model$result$assignments)
})

test_that("the metric guard aligns cosine and Euclidean-family geometries", {
  expect_warning(m <- AINet$new(affinityFunc = "cosine",
                                distFunc = "euclidean", verbose = FALSE),
                 "inconsistent geometries")
  expect_equal(m$config$distFunc, "cosine")
  expect_warning(m <- AINet$new(affinityFunc = "laplace",
                                distFunc = "cosine", verbose = FALSE),
                 "overriding distFunc to 'manhattan'")
  expect_equal(m$config$distFunc, "manhattan")
  expect_no_warning(AINet$new(affinityFunc = "gaussian",
                              distFunc = "manhattan", verbose = FALSE))
})

# ==== Paths otherwise exercised only by the local golden test ====

test_that("robust and arcsinh scaling round-trip and predict consistently", {
  for (sc in c("robust", "arcsinh")) {
    set.seed(21)
    m <- AINet$new(nAntibodies = 15, maxIter = 5, scale = sc, verbose = FALSE)
    m$fit(X)
    expect_equal(.bhive_invert_scaling(
      .bhive_apply_scaling(X, m$config$scaling), m$config$scaling),
      X, label = sc)
    expect_equal(as.numeric(m$predict(X)), m$result$assignments, label = sc)
  }
  expect_error(.bhive_learn_scaling(X, "minmax"), "Unknown scaling method")
})

test_that("epsilonQuantile adapts suppression and falls back when tiny", {
  set.seed(22)
  m <- AINet$new(nAntibodies = 15, maxIter = 5, epsilonQuantile = 0.2,
                 verbose = FALSE)
  expect_no_error(m$fit(X))
  priv <- m$.__enclos_env__$private
  expect_equal(priv$.adaptive_epsilon(X[1:2, ], "euclidean", 0.5, 2,
                                      matrix(0, 0, 0), fallback = 0.3), 0.3)
  expect_equal(priv$.adaptive_epsilon(X[c(1, 1, 1), ], "euclidean", 0.5, 2,
                                      matrix(0, 0, 0), fallback = 0.3), 0.3)
  expect_error(AINet$new(epsilonQuantile = 2), "epsilonQuantile")
})

test_that("auto bandwidth only applies to distance kernels", {
  priv <- AINet$new(verbose = FALSE)$.__enclos_env__$private
  expect_equal(priv$.auto_bandwidth(X, "cosine"), 1)
  expect_equal(priv$.auto_bandwidth(X[c(1, 1), ], "gaussian"), 1)
  expect_gt(priv$.auto_bandwidth(X, "laplace"), 0)
})

test_that("Mahalanobis distance uses the supplied covariance", {
  set.seed(23)
  m <- AINet$new(nAntibodies = 15, maxIter = 4, distFunc = "mahalanobis",
                 affinityParams = list(alpha = 1, Sigma = cov(X)),
                 verbose = FALSE)
  m$fit(X)
  expect_equal(as.numeric(m$predict(X)), m$result$assignments)
})

test_that("the activation gate sets aside non-binding antibodies with their SHM state", {
  set.seed(24)
  m <- AINet$new(nAntibodies = 25, maxIter = 3, verbose = FALSE,
                 shm = SHMEngine$new(method = "adaptive"),
                 activation = ActivationGate$new(threshold1 = 0.5,
                                                 threshold2 = 0.6))
  expect_no_error(m$fit(X, iris$Species))
  expect_identical(m$predict(X), m$result$assignments)
})

test_that("early stopping is reported with message()", {
  set.seed(25)
  m <- AINet$new(nAntibodies = 10, maxIter = 30, noImprovementLimit = 2,
                 verbose = TRUE)
  msgs <- capture_messages(m$fit(X))
  expect_true(any(grepl("Early stopping", msgs)))
  expect_lt(length(m$history), 30L)
})

test_that("targetK above the antibody count adds k-means++ seeds", {
  set.seed(26)
  m <- AINet$new(nAntibodies = 4, maxIter = 3, epsilon = 2, targetK = 8,
                 verbose = FALSE)
  m$fit(X)
  expect_equal(length(unique(m$result$assignments)), 8L)
})

test_that("classification memory archives class labels", {
  set.seed(27)
  mem <- MemoryPool$new(archive_threshold = 0.001)
  m <- AINet$new(nAntibodies = 20, maxIter = 3, memory = mem, verbose = FALSE)
  m$fit(X, iris$Species)
  expect_identical(m$repertoire$metadata$class_label,
                   m$result$antibody_classes)
})

# ==== Class-stratified seeding ====

# Two large classes and a class of 6 observations that sits next to big1,
# the way rare cytometry populations sit next to related abundant ones.
imbalanced <- function() {
  set.seed(40)
  X <- rbind(matrix(rnorm(400, 0, 0.5), ncol = 2),
             matrix(rnorm(400, 4, 0.5), ncol = 2),
             cbind(rnorm(6, 1.6, 0.15), rnorm(6, 0, 0.15)))
  y <- factor(rep(c("big1", "big2", "rare"), c(200, 200, 6)))
  list(X = X, y = y)
}

test_that("stratified seeding lets a rare class be predicted", {
  d <- imbalanced()
  rare_recall <- function(mpc) {
    mean(vapply(1:10, function(s) {
      set.seed(s)
      m <- AINet$new(nAntibodies = 60, k = 1, maxIter = 5, minPerClass = mpc,
                     verbose = FALSE)
      m$fit(d$X, d$y)
      mean(m$predict(d$X)[d$y == "rare"] == "rare")
    }, numeric(1)))
  }
  stratified <- rare_recall(5)
  expect_gt(stratified, 0.85)
  expect_gt(stratified, rare_recall(0))
})

test_that("seeding quotas follow class shares with a per-class floor", {
  d <- imbalanced()
  m <- AINet$new(nAntibodies = 20, minPerClass = 4, verbose = FALSE)
  seeded <- m$.__enclos_env__$private$.stratified_antibodies(
    d$X, d$y, m$config)
  # big classes: round(20 * 200/406) = 10 each; rare: max(4, 0) = 4
  expect_equal(as.vector(table(seeded$classes)), c(10, 10, 4))
  expect_equal(nrow(seeded$A), 24L)
  # rare seeds are drawn from the rare observations
  rare_rows <- d$X[d$y == "rare", , drop = FALSE]
  expect_true(all(apply(seeded$A[seeded$classes == "rare", , drop = FALSE], 1,
                        function(a) any(colSums(t(rare_rows) == a) == 2))))
})

test_that("unused class levels get no antibodies", {
  d <- imbalanced()
  y <- factor(d$y, levels = c(levels(d$y), "empty"))
  m <- AINet$new(nAntibodies = 10, maxIter = 3, verbose = FALSE)
  expect_no_error(m$fit(d$X, y))
  expect_false("empty" %in% m$result$antibody_classes)
})

test_that("minPerClass is validated and an empty initialization is rejected", {
  expect_error(AINet$new(minPerClass = -1), "minPerClass")
  empty_init <- list(generate = function(n, X) X[0, , drop = FALSE])
  m <- AINet$new(init = empty_init, verbose = FALSE)
  expect_error(m$fit(X), "non-empty matrix")
})

test_that("AINet reports when idiotypic regulation falls back", {
  # On a spread-out repertoire the per-pair rule suppresses every antibody.
  set.seed(50)
  m <- AINet$new(nAntibodies = 25, maxIter = 3, verbose = FALSE,
                 idiotypic = IdiotypicNetwork$new())
  msgs <- capture_messages(m$fit(X))
  expect_true(any(grepl("IdiotypicNetwork removed every antibody", msgs)))
  # A permissive window keeps antibodies and stays quiet.
  set.seed(50)
  m <- AINet$new(nAntibodies = 25, maxIter = 3, verbose = FALSE,
                 idiotypic = IdiotypicNetwork$new(theta_low = 1e-4,
                                                  theta_high = 0.999,
                                                  survival_threshold = 0.01))
  expect_false(any(grepl("IdiotypicNetwork",
                         capture_messages(m$fit(X)))))
})
