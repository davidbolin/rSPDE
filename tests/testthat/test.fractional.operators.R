context("fractional.operators")

test_that("Operator construction for fractional stationary Matern", {
  x <- seq(from = 0, to = 1, length.out = 51)
  mesh_1d <- fmesher::fm_mesh_1d(x)
  # fem <- rSPDE.fem1d(x)
  fem <- fmesher::fm_fem(mesh_1d)

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
  v <- t(fmesher::fm_basis(mesh_1d, 0.5))
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
  x <- seq(from = 0, to = 1, length.out = 51)
  mesh_1d <- fmesher::fm_mesh_1d(x)
  # fem <- rSPDE.fem1d(x)
  fem <- fmesher::fm_fem(mesh_1d)

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
  v <- t(fmesher::fm_basis(mesh_1d, 0.5))
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
  x <- seq(from = 0, to = 1, length.out = 51)
  mesh_1d <- fmesher::fm_mesh_1d(x)
  # fem <- rSPDE.fem1d(x)
  fem <- fmesher::fm_fem(mesh_1d)

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
  v <- t(fmesher::fm_basis(mesh_1d, 0.5))
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
  x <- seq(from = 0, to = 1, length.out = 51)
  mesh_1d <- fmesher::fm_mesh_1d(x)
  # fem <- rSPDE.fem1d(x)
  fem <- fmesher::fm_fem(mesh_1d)

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
  v <- t(fmesher::fm_basis(mesh_1d, 0.5))
  c1 <- as.vector(Sigma.mult(op1, v))
  c2 <- as.vector(Sigma.mult(op2, v))
  expect_equal(c1, c2, tolerance = 1e-10)
})
test_that("the non-stationary covariance model has one block per pole", {
  ## Regression test. The bdiag that joins the blocks used to sit outside the
  ## loop over the poles, so every intermediate block was overwritten and the
  ## precision had three blocks whatever the order: for m >= 3 the model was
  ## silently the wrong one, missing the poles 2 to m-1.
  x <- seq(from = 0, to = 1, length.out = 61)
  fem <- rSPDE.fem1d(x)
  n <- length(x)
  d <- 1
  kappa <- 8
  tau <- 0.5
  alpha <- 1.3
  B.tau <- matrix(c(log(tau), 1, 0), 1, 3)
  B.kappa <- matrix(c(log(kappa), 0, 1), 1, 3)
  for (ty in c("brasil", "wl2")) {
    for (m in 1:4) {
      ns <- spde.matern.operators(
        C = fem$C, G = fem$G, d = d, alpha = alpha, m = m,
        B.tau = B.tau, B.kappa = B.kappa, theta = c(0, 0),
        parameterization = "spde", type = "covariance",
        check_stationarity = FALSE, loc_mesh = x,
        type_rational_approximation = ty
      )
      st <- matern.operators(
        C = fem$C, G = fem$G, d = d, alpha = alpha, m = m,
        tau = tau, kappa = kappa, parameterization = "spde",
        type = "covariance", loc_mesh = x,
        type_rational_approximation = ty
      )
      ## with constant coefficients the two constructions are the same model
      expect_equal(nrow(ns$Q), nrow(st$Q))
      expect_equal(nrow(ns$Q) / n, rSPDE:::rspde_n_blocks(m, ty))
      expect_equal(max(abs(ns$Q - st$Q)) / max(abs(st$Q)), 0, tolerance = 1e-12)
    }
  }
})

test_that("the non-stationary models use the rational type they are given", {
  ## type_rational_approximation used not to reach fractional.operators(), so
  ## every type gave the same operator-based model.
  x <- seq(from = 0, to = 1, length.out = 61)
  fem <- rSPDE.fem1d(x)
  n <- length(x)
  B.tau <- matrix(c(log(0.5), 1, 0), 1, 3)
  B.kappa <- matrix(c(log(8), 0, 1), 1, 3)
  sig <- function(ty) {
    op <- spde.matern.operators(
      C = fem$C, G = fem$G, d = 1, alpha = 1.3, m = 2,
      B.tau = B.tau, B.kappa = B.kappa, theta = c(0, 0),
      parameterization = "spde", type = "operator",
      check_stationarity = FALSE, loc_mesh = x,
      type_rational_approximation = ty
    )
    Sigma.mult(op, diag(n))
  }
  expect_gt(max(abs(sig("chebfunLB") - sig("wl2"))), 1e-6)
  ## and it agrees with the stationary builder given the same constant
  ## coefficients
  st <- matern.operators(
    C = fem$C, G = fem$G, d = 1, alpha = 1.3, m = 2, tau = 0.5, kappa = 8,
    parameterization = "spde", type = "operator", loc_mesh = x,
    type_rational_approximation = "chebfunLB"
  )
  expect_equal(
    max(abs(sig("chebfunLB") - Sigma.mult(st, diag(n)))) /
      max(abs(Sigma.mult(st, diag(n)))),
    0,
    tolerance = 1e-10
  )
})

test_that("non-stationary weighted-L2 models are valid", {
  ## genuinely varying kappa, all three covariance classes
  x <- seq(from = 0, to = 1, length.out = 81)
  fem <- rSPDE.fem1d(x)
  n <- length(x)
  B.tau <- cbind(rep(log(0.5), n), 1, 0, 0)
  B.kappa <- cbind(rep(log(8), n), 0, 1, x - 0.5)
  for (nu in c(0.4, 0.8, 2.1)) {
    op <- spde.matern.operators(
      C = fem$C, G = fem$G, d = 1, alpha = nu + 0.5, m = 2,
      B.tau = B.tau, B.kappa = B.kappa, theta = c(0, 0, 0.8),
      parameterization = "spde", type = "covariance",
      check_stationarity = FALSE, loc_mesh = x,
      type_rational_approximation = "wl2"
    )
    expect_equal(nrow(op$Q) / n, 2)
    expect_true(Matrix::isSymmetric(op$Q, tol = 1e-9))
    for (j in 1:2) {
      idx <- (j - 1) * n + seq_len(n)
      expect_gt(
        min(eigen(as.matrix(op$Q[idx, idx]), only.values = TRUE)$values), 0
      )
    }
  }
  ## and floor(alpha) beyond 2 is refused
  expect_error(
    spde.matern.operators(
      C = fem$C, G = fem$G, d = 1, alpha = 3.4, m = 2,
      B.tau = B.tau, B.kappa = B.kappa, theta = c(0, 0, 0.8),
      parameterization = "spde", type = "covariance",
      check_stationarity = FALSE, loc_mesh = x,
      type_rational_approximation = "wl2"
    ),
    "0, 1 or 2"
  )
})

test_that("type = 'operator' has two sets of coefficients, not four", {
  ## The operator-based construction factorises into Pl and Pr, and its
  ## tabulated roots come from a single table (get.roots()), so the three
  ## tabulated names select the same model there. The choice that matters is
  ## that one against "wl2".
  x <- seq(from = 0, to = 1, length.out = 101)
  fem <- rSPDE.fem1d(x)
  n <- length(x)
  sig <- function(ty, m = 3, ...) {
    op <- matern.operators(
      C = fem$C, G = fem$G, d = 1, alpha = 1.3, m = m, tau = 0.5, kappa = 8,
      parameterization = "spde", type = "operator", loc_mesh = x,
      type_rational_approximation = ty, ...
    )
    Sigma.mult(op, diag(n))
  }
  ## the tabulated roots come from one table, produced by "chebfunLB"; the
  ## other two tabulated names are refused rather than quietly given those
  ## roots, since they did not produce them
  base <- sig("chebfunLB")
  expect_error(sig("brasil"), "not available for the operator-based")
  expect_error(sig("chebfun"), "not available for the operator-based")
  expect_gt(max(abs(sig("wl2") - base)), 1e-6)
  ## a caller who does not name a type is not refused: the default is
  ## "brasil", which is right for type = "covariance", and the operator-based
  ## construction falls back to the table it does have
  expect_silent(
    matern.operators(
      C = fem$C, G = fem$G, d = 1, alpha = 1.3, m = 3, tau = 0.5, kappa = 8,
      parameterization = "spde", type = "operator", loc_mesh = x
    )
  )

  ## the tabulated roots are stored for m at most 4; wl2 fits them
  expect_error(sig("chebfunLB", m = 5), "order must be one of")
  expect_silent(sig("wl2", m = 5))

  ## The object has to store the type that was used, not the nominal default:
  ## the operator-based construction resolves an untouched default to
  ## "chebfunLB", and storing "brasil" would make update() refuse the rebuild.
  mk <- function(...) matern.operators(
    C = fem$C, G = fem$G, d = 1, alpha = 1.3, m = 2, tau = 0.5, kappa = 8,
    parameterization = "spde", loc_mesh = x, ...
  )
  op <- mk(type = "operator")
  expect_equal(op$type_rational_approximation, "chebfunLB")
  expect_silent(update(op, kappa = 9))
  op <- mk(type = "operator", type_rational_approximation = "wl2",
           x_min = 1 / 401)
  expect_equal(op$type_rational_approximation, "wl2")
  expect_silent(update(op, kappa = 9))
  ## the covariance type keeps the package default
  op <- mk(type = "covariance")
  expect_equal(op$type_rational_approximation, "brasil")
  expect_silent(update(op, kappa = 9))

  ## the s weight gives a different operator fit: it measures the error of the
  ## solution operator in H^s rather than in L_2
  tabs <- lapply(c(0, 1), function(s) {
    rspde.wl2.table(
      m = 3, d = 1, alpha = 1.3, x_min = 1 / 401, type = "operator", s = s
    )
  })
  expect_gt(
    max(abs(sig("wl2", wl2_table = tabs[[1]]) -
      sig("wl2", wl2_table = tabs[[2]]))),
    1e-8
  )
})
