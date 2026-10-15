# Bit-identical regression guard for AINet$fit(). See helper-golden.R.

test_that("AINet$fit() reproduces the golden outputs", {
  skip_on_cran()
  skip_on_ci()
  golden <- readRDS(test_path("fixtures", "ainet-golden.rds"))
  cfgs <- golden_configs()
  expect_setequal(names(cfgs), names(golden))
  for (nm in names(golden)) {
    now <- run_golden(cfgs[[nm]])
    ref <- golden[[nm]]
    # Compare only the result fields the fixture recorded, so additive fields
    # in later versions do not break the guard.
    now$result <- now$result[names(ref$result)]
    expect_identical(now, ref, label = nm)
  }
})
