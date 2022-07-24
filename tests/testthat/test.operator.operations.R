context("Operator.operations")

test_that("Operator algebra", {
  x <- seq(from = 0, to = 1, length.out = 31)
  fem <- rSPDE.fem1d(x)

  nu <- 0.8
  tau <- seq(from = 1, to = 2, length.out = length(x))
  kappa <- seq(from = 2, to = 20, length.out = length(x))

  op <- spde.matern.operators(
    kappa = kappa, tau = tau, nu = nu,
    G = fem$G, C = fem$C, d = 1
  )
  v <- x


  # Pr multiplication
  expect_equal(op$Pr %*% v, Pr.mult(op, v), tolerance = 1e-10)
  expect_equal(t(op$Pr) %*% v, Pr.mult(op, v, transpose = TRUE),
  tolerance = 1e-10)

  # Pl multiplication
  expect_equal(op$Pl %*% v, Pl.mult(op, v), tolerance = 1e-10)
  expect_equal(t(op$Pl) %*% v, Pl.mult(op, v, transpose = TRUE),
  tolerance = 1e-10)

  # Pr solve
  expect_equal(solve(op$Pr, v), Pr.solve(op, v), tolerance = 1e-10)
  expect_equal(solve(t(op$Pr), v), Pr.solve(op, v, transpose = TRUE),
  tolerance = 1e-10)

  # Pl solve
  expect_equal(solve(op$Pl, v), Pl.solve(op, v), tolerance = 1e-10)
  expect_equal(solve(t(op$Pl), v), Pl.solve(op, v, transpose = TRUE),
  tolerance = 1e-10)

  # Q mult
  expect_equal(op$Q %*% v, Q.mult(op, v), tolerance = 1e-10)
  expect_equal(solve(op$Q, v), Q.solve(op, v), tolerance = 1e-10)

  # Qr mult
  expect_equal(sqrt(op$Ci) %*% op$Pl %*% v, Qsqrt.mult(op, v),
  tolerance = 1e-10)
  expect_equal(t(sqrt(op$Ci) %*% op$Pl) %*% v, Qsqrt.mult(op, v,
  transpose = TRUE), tolerance = 1e-10)

  # Qr solve
  expect_equal(solve(sqrt(op$Ci) %*% op$Pl, v), Qsqrt.solve(op, v),
  tolerance = 1e-10)
  expect_equal(solve(t(sqrt(op$Ci) %*% op$Pl), v), Qsqrt.solve(op, v,
  transpose = TRUE), tolerance = 1e-10)
})
