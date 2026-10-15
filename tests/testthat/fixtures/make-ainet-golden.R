# Regenerate the AINet golden-output fixture.
#
# Run from the package root on a build you trust:
#   Rscript tests/testthat/fixtures/make-ainet-golden.R
# Only regenerate when a change is meant to alter fit() output, and say which
# configs changed and why in the commit message.
pkgload::load_all(".", quiet = TRUE)
source("tests/testthat/helper-golden.R")
golden <- lapply(golden_configs(), run_golden)
saveRDS(golden, "tests/testthat/fixtures/ainet-golden.rds")
message("Wrote ", length(golden), " golden configs")
