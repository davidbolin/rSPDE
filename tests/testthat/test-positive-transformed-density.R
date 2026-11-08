test_that("transformed positive-parameter density keeps positive support", {
  set.seed(101)
  samples <- exp(rnorm(5000, -1, 0.8))
  dens <- .rspde_positive_density(samples, n = 512)
  expect_true(all(dens$x > 0))
  expect_true(all(dens$y >= 0))
  area <- sum(diff(dens$x) *
                (head(dens$y, -1) + tail(dens$y, -1)) / 2)
  expect_equal(area, 1, tolerance = 0.05)
})
