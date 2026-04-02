# tests/testthat/test-downstream-eval.R
# Unit tests for downstream evaluation helpers

test_that("plot_ko_delta rejects invalid result object", {
  expect_error(plot_ko_delta(list()), "delta_scores|non-numeric")
})

test_that("select_advanced_features handles small matrices", {
  set.seed(123)
  mat <- matrix(rnorm(100), nrow = 10, ncol = 10,
                dimnames = list(paste0("S", 1:10), paste0("G", 1:10)))

  # Request more features than available
  result <- select_advanced_features(mat, n_features = 100)

  # Should return at most available columns
  expect_lte(ncol(result), 10)
})
