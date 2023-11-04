context("factional.operators")

test_that("Operator construction for fractional stationary Matern", {
  library(fmesher)
  x <- seq(from = 0, to = 1, length.out = 51)
  mesh_1d <- fm_mesh_1d(x)
  # fem <- rSPDE.fem1d(x)
  fem <- fm_fem(mesh_1d)

  d <- 1
  nu <- 0.8
  sigma <- 0.5
  kappa <- 20
  alpha <- nu + d/2
  range <- sqrt(8*nu)/kappa

  op1 <- matern.operators(
    range = range, sigma = sigma, nu = nu,,
    loc_mesh = x, d = 1,
    type = "operator",
    parameterization = "matern"
  )

  tau <- sqrt(gamma(nu) / (sigma^2 * kappa^(2 * nu) *
  (4 * pi)^(d / 2) * gamma(nu + d / 2)))
  beta <- (nu + d / 2) / 2

  op2 <- spde.matern.operators(
    kappa = kappa, tau = tau, alpha = alpha,
    loc_mesh = x, d = d, type = "operator",
    parameterization = "spde"
  )

  L <- fem$g1 + kappa^2 * fem$c0
  op3 <- fractional.operators(
    L = L, scale.factor = kappa^2, tau = tau,
    beta = beta, C = fem$c0
  )
  # v <- t(rSPDE.A1d(x, 0.5))
  v <- t(fm_basis(mesh_1d, 0.5))
  c1 <- as.vector(Sigma.mult(op1, v))
  c2 <- as.vector(Sigma.mult(op2, v))
  c3 <- as.vector(Sigma.mult(op3, v))
  c0 <- as.vector(matern.covariance(abs(x - 0.5),
  kappa = kappa, nu = nu, sigma = sigma))

  expect_equal(c1, c2, tolerance = 1e-10)
  expect_equal(c2, c3, tolerance = 1e-10)
  expect_equal(c3, c0, tolerance = 0.02)
})

test_that("Operator construction for non-fractional stationary Matern", {
  library(fmesher)
  x <- seq(from = 0, to = 1, length.out = 51)
  mesh_1d <- fm_mesh_1d(x)
  # fem <- rSPDE.fem1d(x)
  fem <- fm_fem(mesh_1d)

  d <- 1
  nu <- 1.5
  sigma <- 0.5
  kappa <- 20
  alpha <- nu + d/2
  range <- sqrt(8*nu)/kappa

  op1 <- matern.operators(
    range = range, sigma = sigma, nu = nu,
    loc_mesh = x, d = 1,
    type = "operator",
    parameterization = "matern"
  )

  tau <- sqrt(gamma(nu) / (sigma^2 * kappa^(2 * nu) *
  (4 * pi)^(d / 2) * gamma(nu + d / 2)))
  beta <- (nu + d / 2) / 2

  op2 <- spde.matern.operators(
    kappa = kappa, tau = tau, alpha = alpha,
    loc_mesh = x, d = d, type = "operator",
    parameterization = "spde"
    
  )

  L <- fem$g1 + kappa^2 * fem$c0
  op3 <- fractional.operators(
    L = L, scale.factor = kappa^2, tau = tau,
    beta = beta, C = fem$c0
  )
  # v <- t(rSPDE.A1d(x, 0.5))
  v <- t(fm_basis(mesh_1d, 0.5))
  c1 <- as.vector(Sigma.mult(op1, v))
  c2 <- as.vector(Sigma.mult(op2, v))
  c3 <- as.vector(Sigma.mult(op3, v))
  c0 <- as.vector(matern.covariance(abs(x - 0.5),
  kappa = kappa, nu = nu, sigma = sigma))

  expect_equal(c1, c2, tolerance = 1e-10)
  expect_equal(c2, c3, tolerance = 1e-10)
  expect_equal(c3, c0, tolerance = 0.02)
})


test_that("Operator construction for fractional
stationary Matern with beta>1", {
  library(fmesher)
  x <- seq(from = 0, to = 1, length.out = 51)
  mesh_1d <- fm_mesh_1d(x)
  # fem <- rSPDE.fem1d(x)
  fem <- fm_fem(mesh_1d)

  d <- 1
  nu <- 2
  sigma <- 0.5
  kappa <- 20
  alpha <- nu + d/2
  range <- sqrt(8*nu)/kappa

  op1 <- matern.operators(
    range = range, sigma = sigma, nu = nu,
    loc_mesh = x, d = 1,
    type = "operator",
    parameterization = "matern"
  )

  tau <- sqrt(gamma(nu) / (sigma^2 * kappa^(2 * nu) *
  (4 * pi)^(d / 2) * gamma(nu + d / 2)))
  beta <- (nu + d / 2) / 2

  op2 <- spde.matern.operators(
    kappa = kappa, tau = tau, alpha = alpha,
    loc_mesh = x, d = d, type = "operator",
    parameterization = "spde"
  )

  L <- fem$g1 + kappa^2 * fem$c0
  op3 <- fractional.operators(
    L = L, scale.factor = kappa^2, tau = tau,
    beta = beta, C = fem$c0
  )
  # v <- t(rSPDE.A1d(x, 0.5))
  v <- t(fm_basis(mesh_1d, 0.5))
  c1 <- as.vector(Sigma.mult(op1, v))
  c2 <- as.vector(Sigma.mult(op2, v))
  c3 <- as.vector(Sigma.mult(op3, v))
  c0 <- as.vector(matern.covariance(abs(x - 0.5),
  kappa = kappa, nu = nu, sigma = sigma))

  expect_equal(c1, c2, tolerance = 1e-10)
  expect_equal(c2, c3, tolerance = 1e-10)
  expect_equal(c3, c0, tolerance = 0.02)
})


test_that("Operator construction for non-stationary Matern", {
  library(fmesher)
  x <- seq(from = 0, to = 1, length.out = 51)
  mesh_1d <- fm_mesh_1d(x)
  # fem <- rSPDE.fem1d(x)
  fem <- fm_fem(mesh_1d)

  d <- 1
  nu <- 0.8
  kappa <- 10 * (1 + 2 * x^2)
  tau <- 0.1 * (1 - 0.7 * x^2)
  alpha <- nu + d/2
  op1 <- spde.matern.operators(
    kappa = kappa, tau = tau, alpha = alpha,
    loc_mesh = x, d = d, m = 1, type = "operator",
    parameterization = "spde"
  )

  beta <- (nu + d / 2) / 2

  L <- fem$g1 + fem$c0 %*% Matrix::Diagonal(dim(fem$c0)[1], kappa^2)
  op2 <- fractional.operators(
    L = L, scale.factor = min(kappa)^2, tau = tau,
    beta = beta, C = fem$c0
  )
  # v <- t(rSPDE.A1d(x, 0.5))
  v <- t(fm_basis(mesh_1d, 0.5))
  c1 <- as.vector(Sigma.mult(op1, v))
  c2 <- as.vector(Sigma.mult(op2, v))
  expect_equal(c1, c2, tolerance = 1e-10)
})
