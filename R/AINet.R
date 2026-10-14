#' @title AINet
#' @description R6 implementation of the Artificial Immune Network algorithm.
#' This is the core bHIVE algorithm using C++ backends for performance-critical
#' operations. Supports composable modules for somatic hypermutation, idiotypic
#' network regulation, germinal center selection, and more.
#'
#' @examples
#' # Clustering with Iris data
#' data(iris)
#' X <- as.matrix(iris[, 1:4])
#' model <- AINet$new(nAntibodies = 15, maxIter = 10, verbose = FALSE)
#' model$fit(X, task = "clustering")
#' table(model$result$assignments)
#'
#' # Classification
#' model2 <- AINet$new(nAntibodies = 20, maxIter = 10, verbose = FALSE)
#' model2$fit(X, iris$Species, task = "classification")
#' mean(model2$result$assignments == as.character(iris$Species))
#'
#' # Predict on new data
#' preds <- model2$predict(X[1:10, ])
#'
#' @return An \code{\link[R6]{R6Class}} generator. \code{AINet$new()} returns
#'   an \code{AINet} object. Call \code{$fit()} on it to train, then
#'   \code{$predict()} to label new data.
#'
#' @importFrom R6 R6Class
#' @importFrom stats rnorm runif sd quantile
#' @export
AINet <- R6::R6Class(
  "AINet",
  inherit = ImmuneAlgorithm,
  public = list(

    #' @description Create a new AINet algorithm instance.
    #' @param nAntibodies Integer. Initial antibody population size.
    #' @param beta Numeric. Clone multiplier.
    #' @param epsilon Numeric. Suppression distance threshold.
    #' @param maxIter Integer. Maximum iterations.
    #' @param k Integer. Top-k antibodies to clone per data point.
    #' @param affinityFunc Character. Affinity function name.
    #' @param distFunc Character. Distance function name.
    #' @param affinityParams List. Parameters for affinity/distance functions.
    #' @param mutationDecay Numeric. Per-iteration mutation rate decay.
    #' @param mutationMin Numeric. Minimum mutation rate.
    #' @param maxClones Numeric. Maximum clones per antibody.
    #' @param stopTolerance Numeric. Early stopping tolerance.
    #' @param noImprovementLimit Integer. Early stopping patience.
    #' @param initMethod Character. Initialization method.
    #' @param shm An SHMEngine instance or NULL for default uniform mutation.
    #' @param init A VDJLibrary instance or NULL for default initialization.
    #' @param activation An ActivationGate instance or NULL.
    #' @param idiotypic An IdiotypicNetwork instance or NULL.
    #' @param germinalCenter A GerminalCenter instance or NULL.
    #' @param microenvironment A Microenvironment instance or NULL.
    #' @param memory A MemoryPool instance or NULL.
    #' @param classSwitcher A ClassSwitcher instance or NULL.
    #' @param consolidate Logical. For clustering, run Lloyd-style
    #'   consolidation (an M-step) after affinity maturation so antibodies are
    #'   pulled onto the data manifold and become true data-space prototypes.
    #'   Has no effect on classification. Default TRUE.
    #' @param consolidationSteps Integer. Maximum consolidation iterations.
    #' @param scale Character. Per-feature input scaling applied at \code{fit()}
    #'   and re-applied to new data at \code{predict()}. One of \code{"none"}
    #'   (default, no transform), \code{"zscore"} (center/SD), \code{"robust"}
    #'   (median / IQR, outlier-tolerant), or \code{"arcsinh"} (inverse
    #'   hyperbolic sine with cofactor \code{scaleCofactor}, the standard
    #'   mass-cytometry transform). Because \code{epsilon}, mutation scale, and
    #'   all distances live in feature units, scaling makes the same defaults
    #'   behave consistently across datasets of different magnitude.
    #' @param scaleCofactor Numeric. Cofactor for \code{scale = "arcsinh"}
    #'   (\code{asinh(x / cofactor)}). Default 5 (CyTOF convention; use ~150 for
    #'   fluorescence flow).
    #' @param targetK Integer or NULL. If set, force the clustering solution to
    #'   exactly \code{targetK} clusters. Affinity maturation still discovers
    #'   where prototypes belong, but the final consolidation seeds a K-means
    #'   (Lloyd) refinement at exactly \code{targetK} centroids: surviving
    #'   antibodies are agglomerated (if more than K) or split with k-means++ (if
    #'   fewer than K) before refinement. This decouples the reported cluster
    #'   count from the emergent suppression dynamics. NULL (default) keeps the
    #'   emergent, self-selected K. Ignored for classification.
    #' @param epsilonQuantile Numeric in (0, 1) or NULL. If set, the suppression
    #'   threshold is recomputed each iteration as this quantile of the pairwise
    #'   distances among the current antibodies, making suppression scale-free
    #'   and adaptive instead of using the fixed \code{epsilon}. NULL (default)
    #'   uses the fixed \code{epsilon}.
    #' @param coverageBoost Logical. Clustering only. After maturation, find
    #'   data points that no surviving antibody covers well (max affinity in the
    #'   bottom \code{coverageQuantile} tail) and seed fresh antibodies there
    #'   with k-means++. Counters the clonal-selection bias toward dense regions,
    #'   which otherwise leaves rare populations unrepresented. Pairs naturally
    #'   with \code{targetK}: the extra seeds give the forced-K refinement
    #'   candidate prototypes for sparse populations. Default FALSE.
    #' @param coverageQuantile Numeric in (0, 1). Affinity-coverage tail that
    #'   defines "poorly covered" points for \code{coverageBoost}. Default 0.05.
    #' @param minPerClass Integer. Classification only. Antibodies are seeded
    #'   per class: each class gets \code{max(minPerClass, round(nAntibodies *
    #'   class share))} antibodies drawn from its own observations with
    #'   \code{initMethod}, and starts with that class label. Without this, a
    #'   rare class (for example a population of a few cells in a cytometry
    #'   sample) often receives no antibody at all and can never be
    #'   predicted. The total can slightly exceed \code{nAntibodies}. Set to 0
    #'   to sample all observations without regard to class. Ignored when an
    #'   \code{init} module is supplied. Default 5.
    #' @param verbose Logical. Report progress with \code{message()}, so it can
    #'   be silenced with \code{suppressMessages()}.
    initialize = function(nAntibodies = 20, beta = 5, epsilon = 0.01,
                          maxIter = 50, k = 3,
                          affinityFunc = "gaussian", distFunc = "euclidean",
                          affinityParams = list(alpha = 1, c = 1, p = 2,
                                                Sigma = NULL),
                          mutationDecay = 1.0, mutationMin = 0.01,
                          maxClones = Inf, stopTolerance = 0.0,
                          noImprovementLimit = Inf, initMethod = "sample",
                          consolidate = TRUE, consolidationSteps = 10L,
                          scale = c("none", "zscore", "robust", "arcsinh"),
                          scaleCofactor = 5, targetK = NULL,
                          epsilonQuantile = NULL, coverageBoost = FALSE,
                          coverageQuantile = 0.05, minPerClass = 5,
                          shm = NULL, init = NULL, activation = NULL,
                          idiotypic = NULL, germinalCenter = NULL,
                          microenvironment = NULL, memory = NULL,
                          classSwitcher = NULL, verbose = TRUE) {
      args <- as.list(environment())
      .ainet_validate_args(args)
      metric <- .ainet_align_metric(affinityFunc, distFunc)
      config <- .ainet_build_config(args, metric)
      # Keep only the modules the caller supplied, in a fixed order.
      modules <- Filter(Negate(is.null), args[.AINET_MODULES])
      super$initialize(config = config, modules = modules)
    },

    #' @description Fit the AINet algorithm to data.
    #' @param X Numeric matrix or data frame (n x d).
    #' @param y Optional factor target for classification.
    #' @param task Character: "clustering" or "classification".
    #'   Inferred from y if NULL.
    #' @param ... Additional arguments (currently unused).
    #' @return Invisible self, with \code{result} populated.
    fit = function(X, y = NULL, task = NULL, ...) {
      # Each stage is a private method. They share one per-fit environment
      # (ctx) and run in a fixed order, because several stages draw random
      # numbers and results are reproducible only under this order.
      ctx <- private$.fit_setup(X, y, task)
      private$.fit_seed_antibodies(ctx)
      private$.fit_kernel_setup(ctx)
      private$.fit_recall_memory(ctx)
      private$.fit_task_setup(ctx)
      private$.fit_shm_setup(ctx)
      private$.fit_loop_setup(ctx)
      for (iter in seq_len(ctx$cfg$maxIter)) {
        if (private$.fit_iteration(ctx, iter)) break
      }
      private$.report_idiotypic_fallbacks(ctx)
      private$.fit_prune_orphans(ctx)
      private$.fit_coverage_boost(ctx)
      private$.fit_final_assignment(ctx)
      private$.fit_archive_memory(ctx)
      invisible(self)
    }
  ),

  private = list(

    # ======================================================================
    # fit() stages
    # ======================================================================
    # ctx is an environment created by .fit_setup(). It holds the scaled data
    # (X), the task bookkeeping, the kernel parameters (iter_alpha, c_p, p_p,
    # Sigma_inv), the SHM moment state (m1, m2) and the per-antibody class
    # labels. Stages mutate it in place so large matrices are not copied.

    # --- 0. Validate, infer task, learn and apply input scaling -------------
    .fit_setup = function(X, y, task) {
      .validate_bHIVE_input(X, y)
      if (is.null(task)) {
        task <- if (is.null(y)) "clustering" else "classification"
      }
      task <- match.arg(task, c("clustering", "classification"))
      X <- as.matrix(X)

      # Learn the transform on this training matrix and stash the stats so
      # predict() applies the identical map to new data. epsilon, mutation
      # scale and every distance are in feature units, so a fixed default
      # only behaves sensibly once features share a comparable scale.
      scaling <- private$.learn_scaling(X, self$config$scale,
                                        self$config$scaleCofactor)
      self$config$scaling <- scaling
      self$history <- list()

      ctx <- new.env(parent = emptyenv())
      ctx$X    <- private$.apply_scaling(X, scaling)
      ctx$y    <- y
      ctx$task <- task
      ctx$n    <- nrow(ctx$X)
      ctx$d    <- ncol(ctx$X)
      ctx$cfg  <- self$config
      ctx
    },

    # --- 1. Initial antibody population -------------------------------------
    # Classification without an init module seeds antibodies per class so
    # rare classes are represented from the start (see minPerClass).
    .fit_seed_antibodies = function(ctx) {
      cfg <- ctx$cfg
      if (ctx$task == "classification" && cfg$minPerClass > 0 &&
          is.null(self$modules$init)) {
        seeded <- private$.stratified_antibodies(ctx$X, ctx$y, cfg)
        ctx$A0 <- seeded$A
        ctx$seed_classes <- seeded$classes
      } else {
        ctx$A0 <- private$.initialize_antibodies(ctx$X, cfg$nAntibodies,
                                                 cfg$initMethod,
                                                 init_lib = self$modules$init)
      }
      if (!is.matrix(ctx$A0) || nrow(ctx$A0) == 0L || ncol(ctx$A0) != ctx$d) {
        stop("Antibody initialization must return a non-empty matrix with ",
             "one column per feature.", call. = FALSE)
      }
      invisible(NULL)
    },

    # Memory recall (clustering only): merge relevant prior memory cells into
    # the starting repertoire, then build it. Runs after the kernel setup so
    # recall uses the resolved bandwidth when alpha = "auto". Disabled for
    # classification because MemoryPool$recall returns cells without their
    # original class labels.
    .fit_recall_memory = function(ctx) {
      A   <- ctx$A0
      mem <- self$modules$memory
      if (ctx$task == "clustering" && !is.null(mem) && mem$size() > 0L) {
        recalled <- mem$recall(ctx$X,
                               affinityFunc = ctx$cfg$affinityFunc,
                               affinityParams = list(alpha = ctx$iter_alpha,
                                                     c = ctx$c_p,
                                                     p = ctx$p_p))
        if (nrow(recalled) > 0 && ncol(recalled) == ctx$d) {
          A <- rbind(A, recalled)
        }
      }
      self$repertoire <- ImmuneRepertoire$new(A)
      ctx$A0 <- NULL
      invisible(NULL)
    },

    # --- 2. Task-specific setup ---------------------------------------------
    .fit_task_setup = function(ctx) {
      ctx$task_int <- switch(ctx$task, clustering = 0L, classification = 1L)
      ctx$nClasses <- 0L
      ctx$classes  <- NULL
      if (ctx$task == "classification") {
        ctx$classes  <- levels(ctx$y)
        ctx$nClasses <- length(ctx$classes)
        ctx$y_num    <- as.numeric(ctx$y) - 1  # 0-indexed class
      } else {
        ctx$y_num <- rep(0, ctx$n)
      }
      invisible(NULL)
    },

    # --- Kernel parameters --------------------------------------------------
    # iter_alpha starts at the base alpha and may be modulated per iteration
    # by ClassSwitcher. alpha = "auto" sets the RBF/Laplace bandwidth by the
    # median heuristic, alpha = 1 / median(||x_i - x_j||^2) over a subsample
    # of the (scaled) data, so the kernel resolves structure at the data's
    # own length scale instead of the fixed alpha = 1 that is arbitrary once
    # features are scaled. Only meaningful for gaussian/laplace kernels.
    .fit_kernel_setup = function(ctx) {
      cfg <- ctx$cfg
      base_alpha <- cfg$affinityParams$alpha %||% 1
      if (is.character(base_alpha) && identical(base_alpha, "auto")) {
        base_alpha <- private$.auto_bandwidth(ctx$X, cfg$affinityFunc)
        if (cfg$verbose) {
          message(sprintf("Auto bandwidth: alpha = %.5g", base_alpha))
        }
      }
      ctx$c_p        <- cfg$affinityParams$c %||% 1
      ctx$p_p        <- cfg$affinityParams$p %||% 2
      ctx$iter_alpha <- base_alpha
      ctx$Sigma_inv  <- if (!is.null(cfg$affinityParams$Sigma)) {
        solve(cfg$affinityParams$Sigma)
      } else {
        matrix(0, 0, 0)
      }
      invisible(NULL)
    },

    # --- SHM dispatch -------------------------------------------------------
    # Parameters are forwarded to clonal_selection_iteration_cpp. Without an
    # SHM module, "uniform" reproduces the legacy decay-based mutation. The
    # adaptive moment matrices are only allocated when used.
    .fit_shm_setup = function(ctx) {
      shm_engine <- self$modules$shm
      ctx$shm_method <- if (is.null(shm_engine)) "uniform" else shm_engine$method
      ctx$shm_p <- if (is.null(shm_engine)) {
        list(c_rate = 1, temperature = 0.5, E_0 = 1, base_rate = 0.1,
             beta1 = 0.9, beta2 = 0.999, adam_epsilon = 1e-8)
      } else {
        shm_engine$params
      }
      ctx$use_adaptive <- identical(ctx$shm_method, "adaptive")
      m <- self$repertoire$size()
      ctx$m1 <- if (ctx$use_adaptive) matrix(0, m, ctx$d) else matrix(0, 0, 0)
      ctx$m2 <- if (ctx$use_adaptive) matrix(0, m, ctx$d) else matrix(0, 0, 0)
      if (ctx$use_adaptive && is.function(shm_engine$init_state)) {
        shm_engine$init_state(m, ctx$d)
      }
      invisible(NULL)
    },

    # --- Early-stopping state and initial class labels ----------------------
    # Classes are pre-allocated so the gated-aside bookkeeping has a value to
    # index into on the first iteration.
    .fit_loop_setup = function(ctx) {
      m <- self$repertoire$size()
      ctx$noImproveCount <- 0
      ctx$prevCount      <- m
      ctx$antibody_classes <- if (ctx$task == "classification") {
        if (length(ctx$seed_classes) == m) {
          ctx$seed_classes
        } else {
          sample(ctx$classes, size = m, replace = TRUE)
        }
      } else {
        rep(NA_character_, m)
      }
      invisible(NULL)
    },

    # --- 3. One maturation iteration. Returns TRUE to stop early -------------
    .fit_iteration = function(ctx, iter) {
      aside <- private$.stage_activation_gate(ctx)
      sel   <- private$.stage_clonal_selection(ctx, aside$A, iter)
      private$.stage_rejoin(ctx, aside, sel)
      private$.stage_germinal_center(ctx)
      env   <- private$.stage_microenvironment(ctx, iter)
      private$.stage_class_switch(ctx, env)
      private$.stage_idiotypic(ctx)
      m_new <- private$.stage_suppression(ctx)
      private$.stage_record(ctx, iter, m_new)
    },

    # Keep the per-antibody bookkeeping (class labels, SHM moments) aligned
    # with the repertoire after it is subset. idx may be negative.
    .subset_aux = function(ctx, idx) {
      if (ctx$task == "classification") {
        ctx$antibody_classes <- ctx$antibody_classes[idx]
      }
      if (ctx$use_adaptive) {
        ctx$m1 <- ctx$m1[idx, , drop = FALSE]
        ctx$m2 <- ctx$m2[idx, , drop = FALSE]
      }
      invisible(NULL)
    },

    # (a0) ActivationGate: antibodies in over-dense neighborhoods (or binding
    # nothing) sit out this round of clonal selection unchanged, then rejoin.
    # Prevents runaway cloning into already-crowded regions of feature space.
    # Returns the selectable matrix A plus the set-aside rows.
    .stage_activation_gate = function(ctx) {
      A_current <- self$repertoire$as_matrix()
      m   <- nrow(A_current)
      out <- list(A = A_current, A_aside = NULL)
      gate <- self$modules$activation
      if (is.null(gate) || m <= 4L) return(out)

      cfg <- ctx$cfg
      Ab_Ab_aff <- compute_affinity_matrix(A_current, A_current,
                                           cfg$affinityFunc,
                                           ctx$iter_alpha, ctx$c_p, ctx$p_p)
      diag(Ab_Ab_aff) <- 0
      density <- rowSums(Ab_Ab_aff)
      # threshold2 in [0, 1] is the density quantile above which an antibody
      # is gated (e.g. 0.75 = top quartile sits out), clamped to [0.5, 0.95].
      q_cut <- max(0.5, min(0.95, gate$threshold2 %||% 0.75))
      gated_idx <- which(density > quantile(density, q_cut, na.rm = TRUE))

      # Signal 1: an antibody whose max affinity to any data point is below
      # threshold1 is also gated out (it isn't binding anything).
      if (!is.null(gate$threshold1) && gate$threshold1 > 0) {
        Ab_X_aff <- compute_affinity_matrix(ctx$X, A_current,
                                            cfg$affinityFunc,
                                            ctx$iter_alpha, ctx$c_p, ctx$p_p)
        low_aff_idx <- which(apply(Ab_X_aff, 2, max) < gate$threshold1)
        gated_idx   <- union(gated_idx, low_aff_idx)
      }

      # At least 2 antibodies must still enter clonal selection.
      if (length(gated_idx) == 0 || (m - length(gated_idx)) < 2L) return(out)
      out$A_aside <- A_current[gated_idx, , drop = FALSE]
      if (ctx$use_adaptive) {
        out$m1 <- ctx$m1[gated_idx, , drop = FALSE]
        out$m2 <- ctx$m2[gated_idx, , drop = FALSE]
      }
      private$.subset_aux(ctx, -gated_idx)
      out$A <- A_current[-gated_idx, , drop = FALSE]
      out
    },

    # (a) Clonal selection + SHM-dispatched mutation [C++]. Returns the
    # selected antibodies and, for classification, their majority labels.
    .stage_clonal_selection = function(ctx, A, iter) {
      cfg   <- ctx$cfg
      shm_p <- ctx$shm_p
      cs <- clonal_selection_iteration_cpp(
        A, ctx$X, ctx$y_num, ctx$task_int, cfg$k, cfg$beta,
        cfg$maxClones, cfg$mutationDecay, cfg$mutationMin,
        iter, cfg$affinityFunc, ctx$iter_alpha, ctx$c_p, ctx$p_p,
        ctx$nClasses, ctx$shm_method,
        shm_p$c_rate, shm_p$temperature, shm_p$E_0, shm_p$base_rate,
        shm_p$beta1, shm_p$beta2, shm_p$adam_epsilon,
        ctx$m1, ctx$m2
      )
      if (ctx$use_adaptive) {
        ctx$m1 <- cs$m1_state
        ctx$m2 <- cs$m2_state
      }
      new_classes <- NULL
      if (ctx$task == "classification") {
        classes  <- ctx$classes
        nClasses <- ctx$nClasses
        new_classes <- apply(cs$class_counts, 1, function(row) {
          if (all(row == 0)) classes[sample(nClasses, 1)]
          else classes[which.max(row)]
        })
      }
      list(A = cs$A, classes = new_classes)
    },

    # Rejoin gated-aside antibodies to the post-selection repertoire. In
    # classification, refresh their labels by majority vote of their nearest
    # data points so stale random labels from pre-allocation don't poison
    # final predictions.
    .stage_rejoin = function(ctx, aside, sel) {
      if (is.null(aside$A_aside)) {
        self$repertoire$cells <- sel$A
        if (ctx$task == "classification") ctx$antibody_classes <- sel$classes
        return(invisible(NULL))
      }
      self$repertoire$cells <- rbind(sel$A, aside$A_aside)
      if (ctx$task == "classification") {
        cfg <- ctx$cfg
        y   <- ctx$y
        ga_aff <- compute_affinity_matrix(ctx$X, aside$A_aside,
                                          cfg$affinityFunc,
                                          ctx$iter_alpha, ctx$c_p, ctx$p_p)
        refreshed <- vapply(seq_len(nrow(aside$A_aside)), function(j) {
          top <- order(ga_aff[, j], decreasing = TRUE)[
            seq_len(min(cfg$k, nrow(ga_aff)))]
          tab <- table(y[top])
          names(tab)[which.max(tab)]
        }, character(1))
        ctx$antibody_classes <- c(sel$classes, refreshed)
      }
      if (ctx$use_adaptive) {
        ctx$m1 <- rbind(ctx$m1, aside$m1)
        ctx$m2 <- rbind(ctx$m2, aside$m2)
      }
      invisible(NULL)
    },

    # (a1) GerminalCenter: Tfh-mediated quality selection. Probabilistic
    # survival weighted by a task-aware quality score (clustering: average
    # affinity to assigned points; classification: majority-class purity).
    .stage_germinal_center = function(ctx) {
      gc_mod <- self$modules$germinalCenter
      if (is.null(gc_mod) || self$repertoire$size() <= gc_mod$nTfh) {
        return(invisible(NULL))
      }
      gc_surv <- gc_mod$select(
        self$repertoire, ctx$X, ctx$y, ctx$task,
        affinityFunc   = ctx$cfg$affinityFunc,
        affinityParams = list(alpha = ctx$iter_alpha, c = ctx$c_p,
                              p = ctx$p_p))
      private$.subset_aux(ctx, gc_surv)
      invisible(NULL)
    },

    # (a2) Microenvironment-aware mutation jitter. Antibodies in over-dense
    # regions get small jitter (stabilize, memory-like). Antibodies in sparse
    # regions get large jitter (explore), countering the clonal-selection
    # drift toward the data centroid. Class labels are preserved. Returns the
    # zone assessment for ClassSwitcher, or NULL when the module is off.
    .stage_microenvironment = function(ctx, iter) {
      microenv <- self$modules$microenvironment
      if (is.null(microenv) || self$repertoire$size() <= 4L) return(NULL)
      cfg <- ctx$cfg
      env <- microenv$assess(self$repertoire, ctx$X,
                             affinityFunc = cfg$affinityFunc,
                             affinityParams = list(alpha = ctx$iter_alpha,
                                                   c = ctx$c_p, p = ctx$p_p))
      A_post   <- self$repertoire$as_matrix()
      x_sd     <- apply(ctx$X, 2, sd)
      decay    <- cfg$mutationDecay ^ max(iter - 1, 0)
      base_amp <- 0.005 * decay  # ~0.5% of feature SD on iter 1, decaying
      mods     <- env$mutation_modifiers
      for (j in seq_len(nrow(A_post))) {
        if (mods[j] <= 0) next
        sigma <- base_amp * mods[j] * x_sd
        A_post[j, ] <- A_post[j, ] + rnorm(ncol(A_post), 0, sigma)
      }
      self$repertoire$cells <- A_post
      env
    },

    # (a3) ClassSwitcher: bind isotype to microenvironment zone and use the
    # population-mean per-isotype alpha for the NEXT iteration's affinity
    # calls. Needs this iteration's zones. The C++ kernels take a scalar
    # alpha, so per-antibody alphas are averaged.
    .stage_class_switch = function(ctx, env) {
      cs_mod <- self$modules$classSwitcher
      if (!is.null(cs_mod) && !is.null(env)) {
        alphas <- cs_mod$switch_isotypes(self$repertoire, env$zones)
        ctx$iter_alpha <- mean(alphas)
      }
      invisible(NULL)
    },

    # (b) Idiotypic regulation [C++]. Bell-curve Ab-Ab dynamics cull
    # antibodies in over-crowded niches (over-stimulation) and isolated ones
    # (under-stimulation), leaving a diversity-preserving repertoire. Runs
    # BEFORE epsilon-ball suppression so the two operators ablate
    # independently. See IdiotypicNetwork for parameter semantics.
    .stage_idiotypic = function(ctx) {
      idi <- self$modules$idiotypic
      if (is.null(idi)) return(invisible(NULL))
      idi_out <- idiotypic_dynamics_cpp(
        self$repertoire$as_matrix(),
        ctx$cfg$affinityFunc, ctx$iter_alpha, ctx$c_p, ctx$p_p,
        idi$theta_low, idi$theta_high, idi$source_rate, idi$decay_rate,
        idi$dt, as.integer(idi$timeSteps), idi$survival_threshold
      )
      surv_idx <- which(as.logical(idi_out$keep))
      # Safety net: if the dynamics would kill every antibody (ill-tuned
      # thresholds for this data scale), keep the top-population 10% so the
      # run continues and a downstream sweep can still score it low. Counted
      # so fit() can report it once rather than acting silently.
      if (length(surv_idx) == 0L) {
        pop    <- as.numeric(idi_out$population)
        keep_n <- max(1L, floor(0.1 * length(pop)))
        surv_idx <- order(pop, decreasing = TRUE)[seq_len(keep_n)]
        ctx$idi_fallbacks <- (ctx$idi_fallbacks %||% 0L) + 1L
      }
      self$repertoire$subset(surv_idx)
      private$.subset_aux(ctx, surv_idx)
      invisible(NULL)
    },

    # (c) Network suppression [C++]. Removes near-duplicate antibodies within
    # an epsilon-ball in distFunc. With epsilonQuantile set, the threshold
    # tracks the current antibody spread so suppression is scale-free.
    # Returns the surviving repertoire size.
    .stage_suppression = function(ctx) {
      cfg <- ctx$cfg
      eps_iter <- cfg$epsilon
      if (!is.null(cfg$epsilonQuantile)) {
        eps_iter <- private$.adaptive_epsilon(
          self$repertoire$cells, cfg$distFunc, cfg$epsilonQuantile,
          ctx$p_p, ctx$Sigma_inv, fallback = cfg$epsilon)
      }
      keep <- network_suppression_cpp(self$repertoire$cells, cfg$distFunc,
                                      eps_iter, ctx$p_p, ctx$Sigma_inv)
      kept_idx <- which(keep)
      self$repertoire$subset(kept_idx)
      private$.subset_aux(ctx, kept_idx)
      m_new <- self$repertoire$size()
      if (m_new == 0) {
        stop("All antibodies were suppressed. Increase nAntibodies or decrease epsilon.")
      }
      m_new
    },

    # Record history, update early stopping, report progress. Returns TRUE
    # when the run should stop.
    .stage_record = function(ctx, iter, m_new) {
      cfg <- ctx$cfg
      self$history[[iter]] <- list(n_antibodies = m_new)
      changeCount <- abs(m_new - ctx$prevCount)
      ctx$noImproveCount <- if (changeCount <= cfg$stopTolerance) {
        ctx$noImproveCount + 1
      } else {
        0
      }
      ctx$prevCount <- m_new
      if (ctx$noImproveCount >= cfg$noImprovementLimit) {
        if (cfg$verbose) {
          message("Early stopping: no improvement for ", ctx$noImproveCount,
                  " iterations.")
        }
        return(TRUE)
      }
      if (cfg$verbose) {
        message(sprintf("Iteration %d | #Antibodies: %d | noImproveCount: %d",
                        iter, m_new, ctx$noImproveCount))
      }
      FALSE
    },

    # The idiotypic rule applies the bell activation to each antibody pair,
    # so in a large or spread-out repertoire every distant antibody counts
    # as a suppressor and the dynamics remove everything. Say so once per
    # fit (see ?IdiotypicNetwork) instead of falling back silently.
    .report_idiotypic_fallbacks = function(ctx) {
      n_fb <- ctx$idi_fallbacks %||% 0L
      if (n_fb > 0L) {
        message(sprintf(paste0(
          "IdiotypicNetwork removed every antibody in %d of %d iterations; ",
          "the top 10%% by population were kept each time. See ",
          "?IdiotypicNetwork for when the module applies."),
          n_fb, length(self$history)))
      }
      invisible(NULL)
    },

    # --- 4a. Orphan-antibody pruning ----------------------------------------
    # Drop antibodies that are not the nearest neighbor of any training
    # point. These "ghost" antibodies survived suppression but bind nothing.
    # Pruning stops test-time mis-assignments and tightens the link between
    # cluster count and effective antibody count.
    .fit_prune_orphans = function(ctx) {
      cfg <- ctx$cfg
      A_final <- self$repertoire$as_matrix()
      fa_pre <- final_assignment_cpp(ctx$X, A_final, cfg$affinityFunc,
                                     cfg$distFunc, ctx$task_int,
                                     ctx$iter_alpha, ctx$c_p, ctx$p_p,
                                     ctx$Sigma_inv)
      assigned_to <- if (ctx$task == "clustering") {
        as.integer(fa_pre$assignments)
      } else {
        as.integer(fa_pre$best_antibody_idx)
      }
      orphan <- tabulate(assigned_to, nbins = nrow(A_final)) == 0L
      if (any(orphan) && sum(!orphan) >= 2L) {
        keep_idx <- which(!orphan)
        A_final  <- A_final[keep_idx, , drop = FALSE]
        private$.subset_aux(ctx, keep_idx)
        self$repertoire$subset(keep_idx)
      }
      ctx$A_final <- A_final
      invisible(NULL)
    },

    # --- 4a'. Coverage boost (clustering only) ------------------------------
    # Clonal expansion follows density, so rare populations can end the run
    # with no nearby prototype. Seed k-means++ antibodies among the
    # worst-covered points (lowest max affinity to any survivor) so the final
    # assignment or forced-K refinement has candidates for those regions.
    .fit_coverage_boost = function(ctx) {
      cfg <- ctx$cfg
      A_final <- ctx$A_final
      if (ctx$task != "clustering" || !isTRUE(cfg$coverageBoost) ||
          nrow(A_final) < 1L) {
        return(invisible(NULL))
      }
      cov_aff <- compute_affinity_matrix(ctx$X, A_final, cfg$affinityFunc,
                                         ctx$iter_alpha, ctx$c_p, ctx$p_p)
      coverage <- apply(cov_aff, 1, max)
      thr <- stats::quantile(coverage, cfg$coverageQuantile, names = FALSE,
                             na.rm = TRUE)
      under <- which(coverage <= thr)
      if (length(under) < 2L) return(invisible(NULL))
      n_new <- min(length(under), max(2L, round(0.5 * cfg$nAntibodies)))
      new_seeds <- init_kmeanspp_cpp(ctx$X[under, , drop = FALSE], n_new)
      ctx$A_final <- rbind(A_final, new_seeds)
      self$repertoire$cells <- ctx$A_final
      if (cfg$verbose) {
        message(sprintf("Coverage boost: +%d antibodies for %d under-covered points",
                        n_new, length(under)))
      }
      invisible(NULL)
    },

    # --- 4b. Final assignment -----------------------------------------------
    # The kernel actually used for the read-out is stored on the config so
    # predict() reproduces the training assignments exactly: the resolved
    # bandwidth (alpha = "auto", ClassSwitcher modulation) and the metric
    # that produced the labels (Euclidean after Lloyd refinement).
    .fit_final_assignment = function(ctx) {
      ctx$assign_dist <- ctx$cfg$distFunc
      if (ctx$task == "clustering") {
        private$.final_clustering(ctx)
      } else {
        private$.final_classification(ctx)
      }
      self$config$kernel <- list(alpha = ctx$iter_alpha, c = ctx$c_p,
                                 p = ctx$p_p, assignDist = ctx$assign_dist)
      invisible(NULL)
    },

    # Clustering read-out. targetK forces the cluster count through a seeded
    # Lloyd refinement. Otherwise consolidation pulls antibodies onto the
    # data manifold so they are genuine data-space prototypes. The repertoire
    # itself stays the matured biological population that memory archives,
    # since consolidated centroids have no 1:1 identity with its rows.
    .final_clustering = function(ctx) {
      cfg <- ctx$cfg
      X <- ctx$X
      A_final <- ctx$A_final
      if (!is.null(cfg$targetK)) {
        out <- private$.force_k(X, A_final, cfg$targetK, cfg, ctx$iter_alpha,
                                ctx$c_p, ctx$p_p, ctx$Sigma_inv)
        ctx$assign_dist <- "euclidean"
      } else if (cfg$consolidate && cfg$consolidationSteps > 0L &&
                 nrow(A_final) >= 2L) {
        out <- private$.consolidate_clusters(X, A_final, cfg, ctx$iter_alpha,
                                             ctx$c_p, ctx$p_p, ctx$Sigma_inv)
        ctx$assign_dist <- "euclidean"
      } else {
        fa <- final_assignment_cpp(X, A_final, cfg$affinityFunc, cfg$distFunc,
                                   0L, ctx$iter_alpha, ctx$c_p, ctx$p_p,
                                   ctx$Sigma_inv)
        # Keep only antibodies that win points so row i is cluster i.
        used <- sort(unique(fa$assignments))
        out  <- list(antibodies  = A_final[used, , drop = FALSE],
                     assignments = as.numeric(match(fa$assignments, used)))
      }
      self$result <- list(
        antibodies  = out$antibodies,
        assignments = out$assignments,
        # Prototypes back-transformed to the original feature space (for
        # interpretation); identical to `antibodies` when scale = "none".
        antibodies_unscaled =
          private$.invert_scaling(out$antibodies, cfg$scaling),
        task        = ctx$task
      )
      invisible(NULL)
    },

    .final_classification = function(ctx) {
      cfg <- ctx$cfg
      A_final <- ctx$A_final
      fa <- final_assignment_cpp(ctx$X, A_final, cfg$affinityFunc,
                                 cfg$distFunc, 1L, ctx$iter_alpha, ctx$c_p,
                                 ctx$p_p, ctx$Sigma_inv)
      self$result <- list(
        antibodies       = A_final,
        assignments      = ctx$antibody_classes[fa$best_antibody_idx],
        antibody_classes = ctx$antibody_classes,
        antibodies_unscaled =
          private$.invert_scaling(A_final, cfg$scaling),
        task             = ctx$task,
        # Row of `antibodies` that labelled each training point.
        antibody_index   = as.integer(fa$best_antibody_idx)
      )
      invisible(NULL)
    },

    # --- 5. Memory archive --------------------------------------------------
    # High-affinity antibodies become long-lived memory cells that persist on
    # the MemoryPool across fit() calls. Classification labels ride along in
    # the repertoire metadata so recall consumers can use them later.
    .fit_archive_memory = function(ctx) {
      mem <- self$modules$memory
      if (is.null(mem)) return(invisible(NULL))
      if (ctx$task == "classification") {
        self$repertoire$metadata$class_label <- ctx$antibody_classes
      }
      mem$archive(self$repertoire, ctx$X,
                  affinityFunc   = ctx$cfg$affinityFunc,
                  affinityParams = list(alpha = ctx$iter_alpha,
                                        c     = ctx$c_p,
                                        p     = ctx$p_p))
      invisible(NULL)
    },

    # ======================================================================
    # Helpers
    # ======================================================================

    # --- Input scaling (C) --------------------------------------------------
    # Thin method wrappers around the file-level scaling helpers so the fit()
    # body reads cleanly; predict() (in the base class) calls the helpers
    # directly since it has no access to these private methods.
    .learn_scaling = function(X, method, cofactor) {
      .bhive_learn_scaling(X, method, cofactor)
    },
    .apply_scaling = function(X, scaling) {
      .bhive_apply_scaling(X, scaling)
    },
    .invert_scaling = function(A, scaling) {
      .bhive_invert_scaling(A, scaling)
    },

    # --- Median-heuristic RBF bandwidth (C) ---------------------------------
    # alpha = 1 / median(||x_i - x_j||^2). Distance-kernel only; identity for
    # cosine/polynomial/hamming where a Euclidean length scale is meaningless.
    .auto_bandwidth = function(X, affinityFunc) {
      if (!affinityFunc %in% c("gaussian", "laplace")) return(1)
      n <- nrow(X)
      idx <- if (n > 1000L) sample.int(n, 1000L) else seq_len(n)
      Xs <- X[idx, , drop = FALSE]
      d2 <- as.vector(stats::dist(Xs))^2
      d2 <- d2[d2 > 0]
      if (length(d2) == 0) return(1)
      med <- stats::median(d2)
      if (!is.finite(med) || med <= 0) return(1)
      1 / med
    },

    # --- Adaptive suppression threshold (C) ---------------------------------
    # Return a low quantile of the pairwise antibody distances so the epsilon
    # ball tracks the repertoire's own spread. Falls back to the fixed epsilon
    # when there are too few antibodies to form a distribution.
    .adaptive_epsilon = function(A, distFunc, q, p, Sigma_inv, fallback) {
      m <- nrow(A)
      if (is.null(m) || m < 3L) return(fallback)
      D <- compute_distance_matrix(A, A, distFunc, p, Sigma_inv)
      du <- D[upper.tri(D)]
      du <- du[is.finite(du) & du > 0]
      if (length(du) == 0) return(fallback)
      as.numeric(stats::quantile(du, probs = q, names = FALSE))
    },

    # --- Force K clusters (A) -----------------------------------------------
    # Use the matured antibodies as informed seeds, coerce to K centroids,
    # then run Euclidean Lloyd. If the network kept more than K prototypes,
    # agglomerate the closest (Ward on the antibodies) down to K. If fewer,
    # add k-means++ seeds drawn from the data so under-clustering can be
    # corrected. The seed assignment is by affinity argmax to the seeds.
    .force_k = function(X, A, K, cfg, alpha, c_p, p_p, Sigma_inv) {
      d <- ncol(X)
      K <- min(K, nrow(X))  # cannot ask for more clusters than points
      m <- nrow(A)
      if (m > K) {
        hc   <- stats::hclust(stats::dist(A), method = "ward.D2")
        grp  <- stats::cutree(hc, k = K)
        cent <- t(vapply(sort(unique(grp)), function(g)
          colMeans(A[grp == g, , drop = FALSE]), numeric(d)))
      } else if (m < K) {
        extra <- init_kmeanspp_cpp(X, K)
        cent  <- rbind(A, extra)[seq_len(K), , drop = FALSE]
      } else {
        cent <- A
      }
      assign <- as.integer(final_assignment_cpp(
        X, cent, cfg$affinityFunc, cfg$distFunc, 1L,
        alpha, c_p, p_p, Sigma_inv)$best_antibody_idx)
      private$.lloyd_refine(X, cent, assign, max(cfg$consolidationSteps, 30L),
                            cfg, alpha, c_p, p_p, Sigma_inv, K = K)
    },

    # Lloyd-style consolidation of matured antibodies into data-space
    # prototypes. Affinity maturation finds where prototypes should point,
    # but under magnitude-blind affinities (cosine) it leaves them off the
    # data manifold, and even Euclidean affinities do not guarantee a
    # prototype equals the centroid of the points it wins.
    #
    # The seed assignment is by AFFINITY (argmax), not distance: a Euclidean
    # seed over off-manifold antibodies can collapse every point onto one
    # far-flung antibody, whereas the affinity seed preserves the partition
    # the network actually learned.
    .consolidate_clusters = function(X, A, cfg, alpha, c_p, p_p, Sigma_inv) {
      seed <- final_assignment_cpp(X, A, cfg$affinityFunc, cfg$distFunc, 1L,
                                   alpha, c_p, p_p, Sigma_inv)$best_antibody_idx
      private$.lloyd_refine(X, A, as.integer(seed), cfg$consolidationSteps,
                            cfg, alpha, c_p, p_p, Sigma_inv)
    },

    # Euclidean Lloyd refinement from a seed partition. The arithmetic mean
    # is the centroid that minimizes within-cluster L2, so reassignment is
    # also Euclidean, which keeps the steps consistent and monotone. This is
    # independent of the training affinity: the affinity drove the search,
    # refinement commits the prototypes to the data manifold.
    #
    # With K set (targetK), a cluster that empties is re-seeded at the point
    # farthest from its own centroid before each M-step, so the partition
    # keeps exactly K clusters instead of silently losing some.
    #
    # The returned labels are the last E-step against the returned
    # prototypes, and only prototypes that win points are kept, so row i of
    # `antibodies` is cluster i and predict() on the training data
    # reproduces `assignments`.
    .lloyd_refine = function(X, cent, assign, steps, cfg, alpha, c_p, p_p,
                             Sigma_inv, K = NULL) {
      prev <- NULL
      for (s in seq_len(steps)) {
        if (!is.null(K)) assign <- .fill_empty_clusters(X, assign, K)
        # M-step: each prototype is the mean of the points assigned to it.
        cent <- .cluster_means(X, assign)
        # E-step: reassign points to the nearest prototype (L2).
        new_assign <- as.integer(final_assignment_cpp(
          X, cent, cfg$affinityFunc, "euclidean", 0L,
          alpha, c_p, p_p, Sigma_inv)$assignments)
        converged <- !is.null(prev) && identical(new_assign, prev)
        prev   <- new_assign
        assign <- new_assign
        if (converged) break
      }
      used <- sort(unique(assign))
      if (!is.null(K) && length(used) < min(K, nrow(X))) {
        # Ran out of steps with a cluster still empty: re-seed it and take a
        # final M-step so the partition has exactly K clusters.
        assign <- .fill_empty_clusters(X, assign, K)
        return(list(antibodies  = .cluster_means(X, assign),
                    assignments = as.numeric(assign)))
      }
      list(antibodies  = cent[used, , drop = FALSE],
           assignments = as.numeric(match(assign, used)))
    },

    # Per-class antibody seeds: quota_c = max(minPerClass, round(nAntibodies *
    # n_c / n)), drawn from class c's own observations with initMethod.
    # Classes with no observations (unused factor levels) get none.
    .stratified_antibodies = function(X, y, cfg) {
      counts  <- table(y)
      present <- names(counts)[counts > 0]
      share   <- as.numeric(counts[present]) / sum(counts)
      quota   <- pmax(as.integer(cfg$minPerClass),
                      as.integer(round(cfg$nAntibodies * share)))
      parts <- lapply(seq_along(present), function(i) {
        Xc <- X[y == present[i], , drop = FALSE]
        private$.initialize_antibodies(Xc, quota[i], cfg$initMethod)
      })
      list(A = do.call(rbind, parts), classes = rep(present, quota))
    },

    .initialize_antibodies = function(X, nAntibodies, method, init_lib = NULL) {
      # When a VDJLibrary (or any object exposing $generate(n, X)) is supplied
      # via the `init` module, route initialization through V(D)J combinatorial
      # assembly. This produces a structured, diverse starting repertoire that
      # spans the data manifold rather than clumping near the centroid.
      if (!is.null(init_lib) && is.function(init_lib$generate)) {
        return(init_lib$generate(nAntibodies, X))
      }

      n <- nrow(X)
      d <- ncol(X)
      switch(
        method,
        "sample" = X[sample.int(n, size = nAntibodies, replace = TRUE), , drop = FALSE],
        "random" = {
          xMean <- colMeans(X)
          xSd   <- apply(X, 2, sd) + 1e-8
          mat   <- matrix(rnorm(nAntibodies * d), nrow = nAntibodies)
          mat   <- sweep(mat, 2, xSd, `*`)
          sweep(mat, 2, xMean, `+`)
        },
        "random_uniform" = {
          xMin <- apply(X, 2, min)
          xMax <- apply(X, 2, max)
          mat  <- matrix(runif(nAntibodies * d), nrow = nAntibodies)
          sweep(sweep(mat, 2, xMax - xMin, `*`), 2, xMin, `+`)
        },
        "kmeans++" = init_kmeanspp_cpp(X, nAntibodies)
      )
    }
  )
)


# ============================================================================
# Internal Lloyd helpers
# ============================================================================

#' Per-cluster means, one row per label in sort(unique(assign))
#' @keywords internal
#' @noRd
.cluster_means <- function(X, assign) {
  ks <- sort(unique(assign))
  t(vapply(ks, function(g) colMeans(X[assign == g, , drop = FALSE]),
           numeric(ncol(X))))
}

#' Re-seed empty clusters so labels 1..K are all present
#'
#' Each missing label takes the point farthest (squared L2) from its own
#' cluster mean, never emptying a singleton. This is the standard k-means fix
#' for empty clusters and keeps a forced-K partition at exactly K.
#' @param assign Integer labels in 1..K, possibly with gaps.
#' @param K Number of clusters wanted. Must not exceed nrow(X).
#' @keywords internal
#' @noRd
.fill_empty_clusters <- function(X, assign, K) {
  missing <- setdiff(seq_len(K), assign)
  for (g in missing) {
    ks    <- sort(unique(assign))
    cent  <- .cluster_means(X, assign)
    own   <- cent[match(assign, ks), , drop = FALSE]
    d2    <- rowSums((X - own)^2)
    sizes <- tabulate(assign, nbins = K)
    d2[sizes[assign] <= 1L] <- -Inf
    assign[which.max(d2)] <- g
  }
  assign
}


# ============================================================================
# Internal constructor helpers for AINet$initialize()
# ============================================================================

# Module slots, in the order they are stored on the model.
.AINET_MODULES <- c("shm", "init", "activation", "idiotypic",
                    "germinalCenter", "microenvironment", "memory",
                    "classSwitcher")

#' Validate AINet numeric arguments
#' @param args Named list of AINet$initialize() arguments.
#' @keywords internal
#' @noRd
.ainet_validate_args <- function(args) {
  with(args, stopifnot(
    "nAntibodies must be a positive integer" = is.numeric(nAntibodies) && nAntibodies >= 1,
    "beta must be positive" = is.numeric(beta) && beta > 0,
    "epsilon must be non-negative" = is.numeric(epsilon) && epsilon >= 0,
    "maxIter must be a positive integer" = is.numeric(maxIter) && maxIter >= 1,
    "k must be a positive integer" = is.numeric(k) && k >= 1,
    "mutationDecay must be in (0, 1]" = is.numeric(mutationDecay) && mutationDecay > 0 && mutationDecay <= 1,
    "mutationMin must be non-negative" = is.numeric(mutationMin) && mutationMin >= 0
  ))
  if (!is.null(args$targetK)) {
    stopifnot("targetK must be a positive integer" =
                is.numeric(args$targetK) && length(args$targetK) == 1L &&
                args$targetK >= 1)
  }
  stopifnot("minPerClass must be a non-negative integer" =
              is.numeric(args$minPerClass) && length(args$minPerClass) == 1L &&
              args$minPerClass >= 0)
  q <- args$epsilonQuantile
  if (!is.null(q)) {
    stopifnot("epsilonQuantile must be in (0, 1)" =
                is.numeric(q) && length(q) == 1L && q > 0 && q < 1)
  }
  invisible(TRUE)
}

#' Resolve and align the affinity and distance metrics
#'
#' Clonal selection matures antibodies to maximize `affinityFunc`, but cluster
#' assignment uses `distFunc`. If the two disagree on geometry the model
#' optimizes one space and is read out in another. The worst case is cosine
#' affinity (scores angle only, leaves antibody magnitude a free random walk)
#' paired with a Euclidean-family distance: antibodies drift far off the data
#' manifold and every point collapses onto one of them. Align the distance to
#' the affinity's natural metric and warn rather than silently returning a
#' degenerate clustering.
#' @return list(affinityFunc, distFunc)
#' @keywords internal
#' @noRd
.ainet_align_metric <- function(affinityFunc, distFunc) {
  affinityFunc <- match.arg(affinityFunc, c("gaussian", "laplace",
                                            "polynomial", "cosine", "hamming"))
  distFunc     <- match.arg(distFunc, c("euclidean", "manhattan",
                                        "minkowski", "cosine",
                                        "mahalanobis", "hamming"))
  natural_dist <- c(gaussian = "euclidean", laplace = "manhattan",
                    polynomial = "euclidean", cosine = "cosine",
                    hamming = "hamming")[[affinityFunc]]
  if ((affinityFunc == "cosine") != (distFunc == "cosine")) {
    warning(sprintf(
      "affinityFunc='%s' and distFunc='%s' use inconsistent geometries; %s. Set distFunc='%s' to silence this.",
      affinityFunc, distFunc,
      sprintf("overriding distFunc to '%s'", natural_dist),
      natural_dist), call. = FALSE)
    distFunc <- natural_dist
  }
  list(affinityFunc = affinityFunc, distFunc = distFunc)
}

#' Build the AINet config list
#' @param args Named list of AINet$initialize() arguments.
#' @param metric Output of .ainet_align_metric().
#' @keywords internal
#' @noRd
.ainet_build_config <- function(args, metric) {
  list(
    nAntibodies        = as.integer(args$nAntibodies),
    beta               = args$beta,
    epsilon            = args$epsilon,
    maxIter            = args$maxIter,
    k                  = args$k,
    affinityFunc       = metric$affinityFunc,
    distFunc           = metric$distFunc,
    affinityParams     = args$affinityParams,
    mutationDecay      = args$mutationDecay,
    mutationMin        = args$mutationMin,
    maxClones          = args$maxClones,
    stopTolerance      = args$stopTolerance,
    noImprovementLimit = args$noImprovementLimit,
    initMethod         = match.arg(args$initMethod,
                                   c("sample", "random", "random_uniform",
                                     "kmeans++")),
    consolidate        = isTRUE(args$consolidate),
    consolidationSteps = as.integer(args$consolidationSteps),
    scale              = match.arg(args$scale,
                                   c("none", "zscore", "robust", "arcsinh")),
    scaleCofactor      = args$scaleCofactor,
    targetK            = if (is.null(args$targetK)) NULL else as.integer(args$targetK),
    epsilonQuantile    = args$epsilonQuantile,
    coverageBoost      = isTRUE(args$coverageBoost),
    coverageQuantile   = args$coverageQuantile,
    minPerClass        = as.integer(args$minPerClass),
    scaling            = NULL,  # populated at fit() with learned stats
    verbose            = args$verbose
  )
}



# ============================================================================
# Internal scaling helpers (shared by AINet$fit and ImmuneAlgorithm$predict)
# ============================================================================
# Kept as free functions, not R6 methods, so the base-class predict() can apply
# the identical transform to new data without reaching into AINet internals.

#' Learn a per-feature scaling from a training matrix
#' @param X numeric matrix
#' @param method one of "none","zscore","robust","arcsinh"
#' @param cofactor arcsinh cofactor
#' @return list(method, center, scale, cofactor) with stats in original units
#' @keywords internal
#' @noRd
.bhive_learn_scaling <- function(X, method = "none", cofactor = 5) {
  if (is.null(method) || method == "none") {
    return(list(method = "none"))
  }
  if (method == "arcsinh") {
    return(list(method = "arcsinh", cofactor = cofactor))
  }
  if (method == "zscore") {
    ctr <- colMeans(X, na.rm = TRUE)
    scl <- apply(X, 2, stats::sd, na.rm = TRUE)
  } else if (method == "robust") {
    ctr <- apply(X, 2, stats::median, na.rm = TRUE)
    scl <- apply(X, 2, stats::IQR, na.rm = TRUE)
  } else {
    stop("Unknown scaling method: ", method)
  }
  # Guard against zero-variance columns (constant features) -> divide by 1.
  scl[!is.finite(scl) | scl <= 0] <- 1
  list(method = method, center = ctr, scale = scl)
}

#' Apply a learned scaling to a matrix
#' @keywords internal
#' @noRd
.bhive_apply_scaling <- function(X, scaling) {
  if (is.null(scaling) || scaling$method == "none") return(X)
  if (scaling$method == "arcsinh") return(asinh(X / scaling$cofactor))
  sweep(sweep(X, 2, scaling$center, `-`), 2, scaling$scale, `/`)
}

#' Invert a learned scaling (map prototypes back to original units)
#' @keywords internal
#' @noRd
.bhive_invert_scaling <- function(A, scaling) {
  if (is.null(scaling) || scaling$method == "none") return(A)
  if (scaling$method == "arcsinh") return(sinh(A) * scaling$cofactor)
  sweep(sweep(A, 2, scaling$scale, `*`), 2, scaling$center, `+`)
}
