test_that("scale_axis works", {
  X <- rnorm(n = 96, mean = 2, sd = 1.5)
  X_scaled <- scale_axis(X)

  expect_equal(min(X_scaled), 0)
  expect_equal(max(X_scaled), 1)
})
