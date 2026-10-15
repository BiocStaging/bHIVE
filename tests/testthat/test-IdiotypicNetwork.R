# Tests for IdiotypicNetwork R6 class

set.seed(42)

# ==== Construction ====

test_that("IdiotypicNetwork initializes with valid defaults", {
  idi <- IdiotypicNetwork$new()
  expect_equal(idi$theta_low, 0.01)
  expect_equal(idi$theta_high, 0.5)
  expect_equal(idi$source_rate, 0.5)
  expect_equal(idi$decay_rate, 0.1)
  expect_equal(idi$dt, 0.1)
  expect_equal(idi$timeSteps, 20L)
  expect_equal(idi$survival_threshold, 0.5)
  expect_null(idi$last_dynamics)
})

test_that("IdiotypicNetwork validates parameters", {
  expect_error(IdiotypicNetwork$new(theta_low = -1), "theta_low")
  expect_error(IdiotypicNetwork$new(theta_low = 0.5, theta_high = 0.3), "theta_high")
  expect_error(IdiotypicNetwork$new(source_rate = -1), "source_rate")
  expect_error(IdiotypicNetwork$new(decay_rate = -1), "decay_rate")
  expect_error(IdiotypicNetwork$new(dt = 0), "dt")
  expect_error(IdiotypicNetwork$new(timeSteps = 0), "timeSteps")
  expect_error(IdiotypicNetwork$new(survival_threshold = -1), "survival_threshold")
})

# ==== Regulation ====

# A deterministic repertoire with known dynamics under the default window
# [theta_low, theta_high] = [0.01, 0.5] and gaussian alpha = 1. Antibodies
# 1-4 sit on a regular simplex with pairwise affinity exp(-2a^2) = 0.25, in
# the middle of the window, so they stimulate each other and survive.
# Antibody 5 is far from all of them (affinity ~0, below theta_low), so it
# is under-stimulated and dies. Random fixtures put almost every pairwise
# affinity outside the window, so the module removed every antibody.
idiotypic_fixture <- function() {
  a <- sqrt(log(4) / 2)
  rbind(cbind(diag(4) * a, 0), c(0, 0, 0, 0, 10))
}

test_that("IdiotypicNetwork keeps the stimulated cluster and culls the isolate", {
  rep <- ImmuneRepertoire$new(idiotypic_fixture())
  idi <- IdiotypicNetwork$new()
  expect_no_warning(idi$regulate(rep, "gaussian"))
  expect_equal(which(idi$last_dynamics$keep), 1:4)
  pop <- as.numeric(idi$get_population())
  expect_true(all(pop[1:4] > idi$survival_threshold))
  expect_lt(pop[5], idi$survival_threshold)
  expect_equal(rep$size(), 4L)
})

test_that("IdiotypicNetwork regulate keeps metadata in sync", {
  rep <- ImmuneRepertoire$new(idiotypic_fixture())
  idi <- IdiotypicNetwork$new()
  expect_no_warning(idi$regulate(rep, "gaussian", list(alpha = 1)))
  expect_equal(rep$size(), nrow(rep$metadata))
  expect_equal(rep$as_matrix(), idiotypic_fixture()[1:4, ])
})

test_that("IdiotypicNetwork get_network returns the Ab-Ab affinities", {
  rep <- ImmuneRepertoire$new(idiotypic_fixture())
  idi <- IdiotypicNetwork$new()
  expect_null(idi$get_network())
  expect_no_warning(idi$regulate(rep, "gaussian"))
  net <- idi$get_network()
  expect_true(is.matrix(net))
  expect_equal(dim(net), c(5L, 5L))
  expect_equal(net[1, 2], 0.25)
  expect_equal(unname(diag(net)), rep(0, 5))
})

test_that("IdiotypicNetwork get_population returns one level per antibody", {
  rep <- ImmuneRepertoire$new(idiotypic_fixture())
  idi <- IdiotypicNetwork$new()
  expect_null(idi$get_population())
  expect_no_warning(idi$regulate(rep, "gaussian"))
  pop <- idi$get_population()
  expect_true(is.numeric(pop))
  expect_length(pop, 5)
})

test_that("IdiotypicNetwork warns when all antibodies removed", {
  # Deliberate wipe-out: a survival threshold no antibody can reach. The
  # repertoire is left unchanged rather than emptied.
  A <- matrix(0.5, nrow = 5, ncol = 3) + matrix(rnorm(15, sd = 0.01), 5, 3)
  rep <- ImmuneRepertoire$new(A)
  idi <- IdiotypicNetwork$new(theta_low = 0.01, theta_high = 0.5,
                               survival_threshold = 10)
  expect_warning(idi$regulate(rep, "gaussian", list(alpha = 1)),
                 "removed all")
  expect_equal(rep$size(), 5L)
})

# ==== Print ====

test_that("IdiotypicNetwork print works before regulation", {
  idi <- IdiotypicNetwork$new()
  expect_output(print(idi), "IdiotypicNetwork")
})

test_that("IdiotypicNetwork print works after regulation", {
  rep <- ImmuneRepertoire$new(idiotypic_fixture())
  idi <- IdiotypicNetwork$new()
  expect_no_warning(idi$regulate(rep, "gaussian"))
  expect_output(print(idi), "4/5 antibodies survived")
})
