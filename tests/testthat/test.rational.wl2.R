## Weighted-L2 rational coefficients: the fit itself, the coefficient tables,
## the caching, and the models built from them with matern.operators(),
## matern.rational() and rspde_lme(). The tests of the INLA interface are in
## test.rational.wl2.inla.R, which needs the compiled cgeneric library.

test_that("weighted-L2 coefficients reproduce the reference values", {
  ## Reference values from the numpy/scipy implementation the R code was
  ## ported from. The tolerance is 2% because the quadrature grids differ
  ## slightly; the largest pole is checked to within 10%.
  cases <- list(
    list(args = list(alpha = 0.75, d = 1, m = 2), err = 1.4265e-2, b = 5.187e1),
    list(args = list(alpha = 0.75, d = 1, m = 4), err = 1.3442e-3, b = 7.925e3),
    list(args = list(alpha = 1.5, d = 2, m = 2), err = 2.2597e-3, b = 8.525e1),
    list(args = list(alpha = 1.5, d = 2, m = 4), err = 8.5549e-5, b = 2.645e3),
    list(
      args = list(alpha = 0.75, d = 1, m = 4, x_min = 1 / 401),
      err = 9.9513e-5, b = 4.524e2
    ),
    list(
      args = list(alpha = 1.5, d = 2, m = 4, x_min = 1 / 801),
      err = 3.2750e-5, b = 1.098e3
    ),
    list(
      args = list(alpha = 1.25, d = 1, m = 3, x_min = 1e-4),
      err = 2.0336e-4, b = 1.593e3
    ),
    list(
      args = list(
        alpha = 1.5, d = 2, m = 2, x_min = 1 / 801, type = "operator"
      ),
      err = 2.7957e-3, b = 3.819e2
    ),
    list(
      args = list(
        alpha = 1.5, d = 2, m = 4, x_min = 1 / 801, type = "operator"
      ),
      err = 4.3432e-5, b = 1.580e3
    ),
    ## floor(alpha) = 2, from the same reference implementation
    ## (handoff/wl2_reference.py of the Rational project, kind "shifted2")
    list(args = list(alpha = 2.4, d = 1, m = 2), err = 7.3610e-5, b = 3.719e1),
    list(args = list(alpha = 2.4, d = 1, m = 4), err = 5.5299e-7, b = 4.044e2),
    list(args = list(alpha = 2.8, d = 2, m = 2), err = 5.8425e-5, b = 1.338e1),
    list(args = list(alpha = 2.8, d = 2, m = 4), err = 3.7251e-7, b = 1.385e2),
    list(
      args = list(alpha = 2.4, d = 2, m = 3, x_min = 1 / 801),
      err = 1.7144e-5, b = 1.820e2
    )
  )
  for (case in cases) {
    cf <- do.call(rational.coefficients.wl2, case$args)
    expect_equal(cf$rel_err, case$err, tolerance = 0.02)
    expect_equal(max(1 - cf$p), case$b, tolerance = 0.1)
  }
})

test_that("weighted-L2 coefficients are in the feasible set", {
  for (alpha in c(0.6, 0.95, 1.1, 1.75, 2.05, 2.6, 2.95)) {
    for (m in 1:4) {
      cf <- rational.coefficients.wl2(alpha, d = 1, m = m)
      expect_true(all(cf$r >= 0))
      expect_true(all(cf$p < 1))
      expect_equal(cf$k, 0)
      expect_equal(length(cf$r), m)
      if (floor(alpha) == 0) {
        expect_null(cf$p0)
      } else {
        expect_true(cf$p0 >= 0 && cf$p0 < 1)
      }
    }
  }
  expect_error(rational.coefficients.wl2(3.5, d = 1, m = 2), "0, 1 or 2")
  expect_error(
    rational.coefficients.wl2(2.5, d = 1, m = 2, type = "operator"),
    "smaller than 2"
  )
})

test_that("the symbol of the weighted-L2 fit approximates x^alpha", {
  cf <- rational.coefficients.wl2(1.5, d = 2, m = 4, x_min = 1e-4)
  x <- exp(seq(log(1e-4), 0, length.out = 200))
  expect_equal(rSPDE:::wl2_symbol(cf, x), x^1.5, tolerance = 1e-2)
})

test_that("the coefficient table interpolates as well as a direct fit", {
  for (m_alpha in 1:2) {
    tab <- rSPDE:::wl2_coefficient_table(d = 1, m = 4, m_alpha = m_alpha)
    expect_equal(rSPDE:::wl2_q(attr(tab, "kind")), m_alpha)
    for (alpha in m_alpha + c(0.113, 0.457, 0.802)) {
      ci <- rSPDE:::wl2_interp_coefficients(tab, alpha)
      direct <- rational.coefficients.wl2(alpha, d = 1, m = 4)
      expect_equal(ci$rel_err, direct$rel_err, tolerance = 0.05)
      expect_equal(ci$q, m_alpha)
      expect_true(all(ci$r > 0))
      expect_true(all(ci$p < 1))
    }
  }
})

test_that("the weighted-L2 blocks reproduce the rational symbol on a mesh", {
  ## Each block must be symmetric positive definite, and the covariance of the
  ## approximation must equal V diag(r(x_j)) V^T from the generalised
  ## eigendecomposition of (G, C).
  x <- seq(from = 0, to = 1, length.out = 101)
  n <- length(x)
  fem <- rSPDE.fem1d(x)
  Cl <- Matrix::Diagonal(n, rowSums(fem$C))
  kappa <- 15
  sigma <- 1
  d <- 1
  for (nu in c(0.25, 1, 1.9)) {
    alpha <- nu + d / 2
    m <- 2
    op <- matern.operators(
      loc_mesh = x, nu = nu, range = sqrt(8 * nu) / kappa, sigma = sigma,
      d = d, m = m, parameterization = "matern",
      type_rational_approximation = "wl2"
    )
    ## m blocks, not m + 1
    expect_equal(dim(op$Q)[1], m * n)
    expect_true(Matrix::isSymmetric(op$Q, tol = 1e-9))
    for (j in seq_len(m)) {
      idx <- (j - 1) * n + seq_len(n)
      Qj <- op$Q[idx, idx, drop = FALSE]
      expect_true(Matrix::isSymmetric(Qj, tol = 1e-9))
      expect_true(min(eigen(as.matrix(Qj), only.values = TRUE)$values) > 0)
    }

    L <- (fem$G + kappa^2 * Cl) / kappa^2
    ev <- eigen(as.matrix(Matrix::solve(Cl, L)))
    lambda <- Re(ev$values)
    V <- Re(ev$vectors)
    V <- sweep(V, 2, sqrt(diag(t(V) %*% as.matrix(Cl) %*% V)), "/")
    tau <- sqrt(gamma(nu) / (sigma^2 * kappa^(2 * nu) *
      (4 * pi)^(d / 2) * gamma(nu + d / 2)))
    cf <- rSPDE:::wl2_interp_coefficients(op$wl2_table, alpha)
    symbol <- rSPDE:::wl2_symbol(cf, 1 / lambda)
    Sigma_ref <- V %*% diag(symbol) %*% t(V) / (tau^2 * kappa^(2 * alpha))

    A <- kronecker(matrix(1, 1, m), Matrix::Diagonal(n))
    Sigma <- as.matrix(A %*% Matrix::solve(op$Q, t(A)))
    expect_equal(max(abs(Sigma - Sigma_ref)), 0, tolerance = 1e-8)
  }
})

test_that("the weighted-L2 blocks have the same sparsity as the tabulated ones", {
  x <- seq(from = 0, to = 1, length.out = 51)
  n <- length(x)
  for (nu in c(0.25, 1)) {
    op_wl2 <- matern.operators(
      loc_mesh = x, nu = nu, range = 0.2, sigma = 1, d = 1, m = 2,
      parameterization = "matern", type_rational_approximation = "wl2"
    )
    op_tab <- matern.operators(
      loc_mesh = x, nu = nu, range = 0.2, sigma = 1, d = 1, m = 2,
      parameterization = "matern", type_rational_approximation = "brasil"
    )
    nnz <- function(Q, j, n) {
      idx <- (j - 1) * n + seq_len(n)
      sum(Q[idx, idx, drop = FALSE] != 0)
    }
    expect_equal(nnz(op_wl2$Q, 1, n), nnz(op_tab$Q, 1, n))
  }
})

test_that("the weighted-L2 covariance beats the tabulated one", {
  ## d = 1, m = 4, kappa * h = 0.1, compared with the exact covariance of the
  ## discretised model, for one nu in each of the three classes.
  d <- 1
  m <- 4
  sigma <- 1
  kappa <- 20
  h <- 0.1 / kappa
  x <- seq(from = 0, to = 1, by = h)
  n <- length(x)
  fem <- rSPDE.fem1d(x)
  Cl <- Matrix::Diagonal(n, rowSums(fem$C))
  L <- (fem$G + kappa^2 * Cl) / kappa^2
  ev <- eigen(as.matrix(Matrix::solve(Cl, L)))
  lambda <- Re(ev$values)
  V <- Re(ev$vectors)
  V <- sweep(V, 2, sqrt(diag(t(V) %*% as.matrix(Cl) %*% V)), "/")

  for (nu in c(0.25, 0.9, 1.9)) {
    alpha <- nu + d / 2
    tau <- sqrt(gamma(nu) / (sigma^2 * kappa^(2 * nu) *
      (4 * pi)^(d / 2) * gamma(nu + d / 2)))
    Sigma_exact <- V %*% diag(lambda^(-alpha)) %*% t(V) /
      (tau^2 * kappa^(2 * alpha))

    err <- sapply(c("wl2", "brasil"), function(ty) {
      op <- matern.operators(
        loc_mesh = x, nu = nu, range = sqrt(8 * nu) / kappa, sigma = sigma,
        d = d, m = m, parameterization = "matern",
        type_rational_approximation = ty
      )
      A <- kronecker(
        matrix(1, 1, rSPDE:::rspde_n_blocks_obj(op)), Matrix::Diagonal(n)
      )
      Sigma <- as.matrix(A %*% Matrix::solve(op$Q, t(A)))
      max(abs(Sigma - Sigma_exact)) / sigma^2
    })
    expect_lt(err[["wl2"]], err[["brasil"]])
  }
})

test_that("the nodal variance correction removes the variance deficit", {
  ## The correction fills a deficit and never subtracts: away from the boundary
  ## the nodal variance becomes sigma^2 exactly, while within about one range
  ## of the boundary the model already has more variance than sigma^2 and is
  ## left alone (the exact covariance of the Neumann problem is the folded
  ## Matern, whose variance doubles there).
  x <- seq(from = 0, to = 4, length.out = 401)
  mk <- function(vc) {
    matern.operators(
      loc_mesh = x, nu = 0.3, range = 0.2, sigma = 1, d = 1, m = 2,
      parameterization = "matern", type_rational_approximation = "wl2",
      variance_correction = vc
    )
  }
  op <- mk("nodal")
  expect_equal(op$variance_correction, "nodal")
  expect_true(all(op$nodal_correction >= 0))
  inner <- which(x > 1 & x < 3)
  v0 <- as.vector(Matrix::diag(covariance_mesh(mk("none"))))
  v1 <- as.vector(Matrix::diag(covariance_mesh(op)))
  ## there is a deficit to remove, and it is removed in the interior
  expect_gt(max(1 - v0[inner]), 1e-3)
  expect_equal(v1[inner], rep(1, length(inner)), tolerance = 1e-8)
  ## and nothing is taken away anywhere
  expect_true(all(v1 >= v0 - 1e-12))
})

test_that("the variance correction is off by default and only where it means something", {
  x <- seq(from = 0, to = 4, length.out = 201)
  mk <- function(...) {
    matern.operators(
      loc_mesh = x, nu = 0.3, range = 0.2, sigma = 1, d = 1, m = 2,
      parameterization = "matern", type = "covariance", ...
    )
  }
  ## off unless asked for, and then the default model is exactly the
  ## uncorrected one
  default <- mk(type_rational_approximation = "wl2")
  expect_equal(default$variance_correction, "none")
  expect_null(default$nodal_correction)
  expect_equal(
    as.matrix(covariance_mesh(default)),
    as.matrix(covariance_mesh(mk(
      type_rational_approximation = "wl2", variance_correction = "none"
    )))
  )
  ## and it survives an update in either state
  expect_equal(update(default, nu = 0.35)$variance_correction, "none")
  on <- mk(type_rational_approximation = "wl2", variance_correction = "nodal")
  expect_equal(update(on, nu = 0.35)$variance_correction, "nodal")

  ## It is specific to the covariance-based weighted-L2 construction: the
  ## tabulated classes have a constant term already, and the operator-based
  ## construction forms the covariance as r(x)^2 from factors, which a
  ## diagonal addition does not fit. Both combinations warn rather than
  ## quietly doing nothing.
  expect_warning(
    mk(type_rational_approximation = "brasil", variance_correction = "nodal"),
    "only used with type_rational_approximation = 'wl2'"
  )
  expect_warning(
    matern.operators(
      loc_mesh = x, nu = 0.3, range = 0.2, sigma = 1, d = 1, m = 2,
      parameterization = "matern", type = "operator",
      type_rational_approximation = "wl2", variance_correction = "nodal"
    ),
    "only used with type = 'covariance'"
  )
})

test_that("the corrected covariance is still a covariance", {
  ## Targeting sigma^2 at the boundary, where the model correctly has about
  ## 2 sigma^2, used to subtract about 1 from the diagonal there and took the
  ## smallest eigenvalue from about +0.05 to about -0.9.
  x <- seq(from = 0, to = 4, length.out = 201)
  for (nu in c(0.3, 0.6, 1.2)) {
    for (m in c(1, 2, 4)) {
      mineig <- function(vc) {
        S <- as.matrix(covariance_mesh(matern.operators(
          loc_mesh = x, nu = nu, range = 0.2, sigma = 1, d = 1, m = m,
          parameterization = "matern", type_rational_approximation = "wl2",
          variance_correction = vc
        )))
        min(eigen((S + t(S)) / 2, symmetric = TRUE, only.values = TRUE)$values)
      }
      expect_gt(mineig("nodal"), 0)
      ## the correction only adds a non-negative diagonal, so it cannot make
      ## the smallest eigenvalue smaller
      expect_gte(mineig("nodal"), mineig("none") - 1e-10)
    }
  }
})

test_that("the nodal correction improves the covariance where there is a deficit", {
  ## Against the folded Matern on the interior of a domain long relative to the
  ## range, so that sigma^2 is the right nodal variance there. The gain is
  ## largest at low smoothness, where a mesh represents a rough field worst:
  ## the deficit is mostly finite element discretisation error, not rational
  ## approximation error, so it does not shrink with m.
  skip_on_cran()
  L <- 4
  n <- 401
  x <- seq(0, L, length.out = n)
  h <- x[2] - x[1]
  cl <- rep(h, n)
  cl[c(1, n)] <- h / 2
  inner <- which(x > 1 & x < L - 1)
  nu <- 0.3
  rho <- 0.2
  kappa <- sqrt(8 * nu) / rho
  target <- matrix(0, n, n)
  for (k in -8:8) {
    for (sg in c(1, -1)) {
      target <- target + matern.covariance(
        abs(outer(x, sg * x + 2 * L * k, "-")),
        kappa = kappa, nu = nu, sigma = 1
      )
    }
  }
  ci <- sqrt(cl[inner])
  err <- function(vc, m) {
    S <- as.matrix(covariance_mesh(matern.operators(
      loc_mesh = x, nu = nu, range = rho, sigma = 1, d = 1, m = m,
      parameterization = "matern", type_rational_approximation = "wl2",
      variance_correction = vc
    )))
    E <- S[inner, inner] - target[inner, inner]
    c(L2 = sqrt(sum((ci * t(ci * E))^2)), sup = max(abs(E)))
  }
  for (m in c(2, 4)) {
    a <- err("none", m)
    b <- err("nodal", m)
    ## measured: L2 3.5x (m = 2) and 10x (m = 4), sup 13.5x and 18.6x
    expect_gt(a[["L2"]] / b[["L2"]], 2.5)
    expect_gt(a[["sup"]] / b[["sup"]], 8)
  }
})

test_that("the operator-based conversion reproduces the rational symbol", {
  for (args in list(
    list(alpha = 1.5, d = 2, m = 2), list(alpha = 1.5, d = 2, m = 4),
    list(alpha = 0.75, d = 1, m = 2), list(alpha = 1.9, d = 1, m = 3)
  )) {
    cf <- do.call(rational.coefficients.wl2, c(args, list(
      x_min = 1 / 801, type = "operator"
    )))
    roots <- rSPDE:::wl2_roots(cf)
    expect_equal(length(roots$rb), args$m + 1)
    expect_equal(length(roots$rc), args$m)
    expect_true(all(roots$rb < 0))
    expect_true(all(roots$rc < 0))
    ## the numerator roots interlace the poles
    expect_true(all(diff(sort(c(1 / roots$rb, 1 / roots$rc))) != 0))
    lambda <- exp(seq(0, log(801), length.out = 100))
    symbol <- roots$factor *
      apply(outer(lambda, roots$rc, function(l, r) 1 - r * l), 1, prod) /
      apply(outer(lambda, roots$rb, function(l, r) 1 - r * l), 1, prod)
    expect_equal(symbol, lambda^(-args$alpha / 2), tolerance = 5e-2)
  }
})

test_that("operator-based models can be built with the weighted-L2 type", {
  x <- seq(from = 0, to = 1, length.out = 101)
  op <- matern.operators(
    loc_mesh = x, nu = 0.4, range = 0.2, sigma = 1, d = 1, m = 2,
    parameterization = "matern", type = "operator",
    type_rational_approximation = "wl2"
  )
  expect_equal(op$type_rational_approximation, "wl2")
  expect_true(op$x_min > 0 && op$x_min < 1)
  expect_equal(length(simulate(op)), length(x))
  ## a user-supplied lower bound for kappa tightens the interval
  op2 <- matern.operators(
    loc_mesh = x, nu = 0.4, range = 0.2, sigma = 1, d = 1, m = 2,
    parameterization = "matern", type = "operator",
    type_rational_approximation = "wl2", kappa_ref = sqrt(8 * 0.4) / 0.6
  )
  expect_true(op2$x_min > op$x_min)
  ## alpha >= 2 is out of scope for the operator type
  expect_error(matern.operators(
    loc_mesh = x, nu = 1.8, range = 0.2, sigma = 1, d = 1, m = 2,
    parameterization = "matern", type = "operator",
    type_rational_approximation = "wl2"
  ))
})

test_that("covariance-based methods work with the weighted-L2 type", {
  set.seed(123)
  x <- seq(from = 0, to = 1, length.out = 101)
  op <- matern.operators(
    loc_mesh = x, nu = 0.8, range = 0.2, sigma = 1, d = 1, m = 2,
    parameterization = "matern", type_rational_approximation = "wl2",
    compute_logdet = TRUE
  )
  A <- Matrix::Diagonal(length(x))
  Y <- as.vector(simulate(op)) + rnorm(length(x), sd = 0.1)
  expect_equal(length(simulate(op)), length(x))
  ll <- rSPDE.matern.loglike(op, Y, A = A, sigma.e = 0.1)
  expect_true(is.finite(ll))
  ## the log-determinant shortcut must agree with the direct computation
  op2 <- update(op, compute_higher_order = FALSE)
  op2$compute_logdet <- FALSE
  expect_equal(ll, rSPDE.matern.loglike(op2, Y, A = A, sigma.e = 0.1),
    tolerance = 1e-8
  )
  pred <- predict(op, A = A, Aprd = A, Y = Y, sigma.e = 0.1,
    compute.variances = TRUE
  )
  expect_equal(length(pred$mean), length(x))
  expect_true(all(pred$variance > 0))
  expect_equal(dim(covariance_mesh(op)), c(length(x), length(x)))
  expect_equal(
    dim(cov_function_mesh(op, p = matrix(0.5, 1, 1))), c(length(x), 1)
  )
  ## floor(alpha) >= 3 is out of scope
  expect_error(matern.operators(
    loc_mesh = x, nu = 2.8, range = 0.2, sigma = 1, d = 1, m = 2,
    parameterization = "matern", type_rational_approximation = "wl2"
  ), "0, 1 or 2")
})

test_that("rspde.xmin is a lower bound for the spectral interval", {
  x <- seq(from = 0, to = 1, length.out = 401)
  fem <- rSPDE.fem1d(x)
  Cl <- Matrix::Diagonal(length(x), rowSums(fem$C))
  kappa <- 10
  x_min <- rspde.xmin(C = Cl, G = fem$G, kappa_ref = kappa)
  L <- (fem$G + kappa^2 * Cl) / kappa^2
  lambda_max <- max(Re(eigen(as.matrix(Matrix::solve(Cl, L)),
    only.values = TRUE
  )$values))
  expect_lte(x_min, 1 / lambda_max)
  ## a smaller reference kappa gives a wider interval
  expect_lt(rspde.xmin(C = Cl, G = fem$G, kappa_ref = 1), x_min)
})

test_that("the weighted-L2 type is rejected where it is not implemented", {
  expect_error(rSPDE:::get_rational_coefficients(2, "wl2"), "computed at")
})

test_that("updating a model keeps the weighted-L2 settings", {
  x <- seq(from = 0, to = 1, length.out = 51)
  for (type in c("covariance", "operator")) {
    op <- matern.operators(
      loc_mesh = x, nu = 0.4, range = 0.2, sigma = 1, d = 1, m = 2,
      parameterization = "matern", type = type,
      type_rational_approximation = "wl2"
    )
    op2 <- update(op, range = 0.1)
    expect_equal(op2$type_rational_approximation, "wl2")
    expect_equal(op2$x_min, op$x_min)
    ## the two-stage helper tightens the interval
    op3 <- update_rational_coefficients(op)
    expect_equal(op3$type_rational_approximation, "wl2")
    expect_true(op3$x_min > (if (is.null(op$x_min)) 0 else op$x_min))
  }
  op_tab <- matern.operators(
    loc_mesh = x, nu = 0.4, range = 0.2, sigma = 1, d = 1, m = 2,
    parameterization = "matern", type = "operator"
  )
  expect_error(update_rational_coefficients(op_tab), "wl2")
})

test_that("rspde_lme works with the weighted-L2 type", {
  skip_on_cran()
  set.seed(42)
  x <- seq(from = 0, to = 1, length.out = 101)
  op_true <- matern.operators(
    loc_mesh = x, nu = 0.8, range = 0.3, sigma = 1, d = 1, m = 2,
    parameterization = "matern"
  )
  u <- as.vector(simulate(op_true))
  loc <- runif(80)
  Y <- as.vector(rSPDE.A1d(x, loc) %*% u) + rnorm(80, sd = 0.1)
  df <- data.frame(y = Y, x = loc)
  op <- matern.operators(
    loc_mesh = x, nu = 0.8, range = 0.3, sigma = 1, d = 1, m = 2,
    parameterization = "matern", type_rational_approximation = "wl2"
  )
  ## the smoothness bound is reduced so that alpha stays below 3
  expect_message(
    fit <- rspde_lme(y ~ 1, loc = "x", data = df, model = op, parallel = FALSE),
    "alpha < 3"
  )
  ## and nu stays inside what the classes support, floor(alpha) at most 2
  expect_lt(fit$coeff$random_effects[["nu"]] + 1 / 2, 3)
  expect_true(is.finite(fit$loglik))
  pred <- predict(fit, newdata = data.frame(x = c(0.25, 0.75)), loc = "x")
  expect_equal(length(as.vector(pred$mean)), 2)
})

test_that("the shifted term equals the tabulated one when p0 = 0", {
  ## The weighted-L2 term 1/((lambda - p0)(lambda - p)) reduces to the term
  ## lambda^-1/(lambda - p) of the tabulated classes when p0 = 0, so the
  ## generalised code must reproduce the existing one exactly.
  h <- seq(from = 0, to = 0.6, length.out = 9)
  kappa <- 7
  for (alpha in c(1.2, 1.7)) {
    for (p in c(-0.4, -3.2, -50)) {
      for (deriv in 0:2) {
        old <- sapply(h, function(hh) {
          rSPDE:::matern.p.deriv(hh, 0, kappa, p, alpha, deriv = deriv)
        })
        new <- sapply(h, function(hh) {
          rSPDE:::matern.p.deriv(hh, 0, kappa, p, alpha, deriv = deriv, p0 = 0)
        })
        expect_equal(new, old, tolerance = 1e-10)
      }
    }
  }
})

test_that("the one-dimensional models work with the weighted-L2 type", {
  set.seed(7)
  sigma <- 1.3
  kappa <- 9
  for (nu in c(0.3, 0.9, 1.7)) {
    ## The joint covariance of the field and its derivatives at two nearly
    ## coincident locations is singular to working precision once the field
    ## has two derivatives, for the tabulated types just as much as for "wl2",
    ## so the smoothest case is taken on a grid rather than on random points.
    loc <- if (nu > 1) seq(from = 0, to = 1, length.out = 40) else sort(runif(40))
    n <- length(loc)
    for (ty in c("brasil", "wl2")) {
      tmp <- rSPDE:::matern.rational.ldl(
        loc = loc, order = 2, nu = nu, kappa = kappa, sigma = sigma,
        type_rational = ty
      )
      ## m blocks for "wl2", m + 1 for the tabulated types; the blocks hold
      ## the process and its derivative when alpha > 1
      ## each pole block holds floor(alpha) + 1 entries per location, and the
      ## k block of the tabulated types holds max(floor(alpha), 1)
      fa <- floor(nu + 1 / 2) + 1
      expected <- 2 * fa * n
      if (ty != "wl2") {
        expected <- expected + max(floor(nu + 1 / 2), 1) * n
      }
      expect_equal(dim(tmp$L)[1], expected)

      ## the Markov construction must agree with the closed-form covariance
      Q <- Matrix::t(tmp$L) %*% tmp$D %*% tmp$L
      Sigma <- as.matrix(tmp$A %*% Matrix::solve(Q, Matrix::t(tmp$A)))
      Sigma_ref <- matern.rational.cov(
        as.matrix(dist(loc)), order = 2, kappa = kappa, nu = nu,
        sigma = sigma, type_rational = ty
      )
      expect_equal(max(abs(Sigma - Sigma_ref)), 0, tolerance = 1e-7)

      ## the location ordering is a permutation of the field ordering
      a <- rSPDE:::matern.rational.precision(
        loc = loc, order = 2, nu = nu, kappa = kappa, sigma = sigma,
        type_rational = ty, ordering = "field"
      )
      b <- rSPDE:::matern.rational.precision(
        loc = loc, order = 2, nu = nu, kappa = kappa, sigma = sigma,
        type_rational = ty, ordering = "location"
      )
      expect_equal(
        max(abs(as.matrix(a$A %*% Matrix::solve(a$Q, Matrix::t(a$A))) -
          as.matrix(b$A %*% Matrix::solve(b$Q, Matrix::t(b$A))))),
        0,
        tolerance = 1e-8
      )
    }
  }
})

test_that("matern.rational supports the weighted-L2 type", {
  s <- seq(from = 0, to = 1, length.out = 101)
  kappa <- 20
  sigma <- 2
  for (nu in c(0.3, 0.8, 1.8)) {
    err <- sapply(c("wl2", "brasil"), function(ty) {
      op <- matern.rational(
        loc = s, nu = nu, range = sqrt(8 * nu) / kappa, sigma = sigma, m = 2,
        parameterization = "matern", type_rational_approximation = ty
      )
      expect_equal(op$n_blocks, if (ty == "wl2") 2 else 3)
      expect_equal(length(simulate(op)), length(s))
      ## each pole block holds floor(alpha) + 1 entries per location, and the
      ## k block of the tabulated types holds max(floor(alpha), 1)
      alpha <- nu + 0.5
      fa_p <- floor(alpha) + 1
      expected <- 2 * fa_p * length(s)
      if (ty != "wl2") {
        expected <- expected + max(floor(alpha), 1) * length(s)
      }
      expect_equal(dim(precision(op)$Q)[1], expected)
      max(abs(op$covariance(ind = 1) -
        matern.covariance(abs(s - s[1]), kappa = kappa, sigma = sigma, nu = nu)))
    })
    expect_lt(err[["wl2"]], err[["brasil"]])
  }
  ## alpha >= 3 is out of scope
  expect_error(matern.rational(
    loc = s, nu = 2.8, range = 0.2, sigma = 1, m = 2,
    parameterization = "matern", type_rational_approximation = "wl2"
  ), "floor")
})

test_that("matern.rational.cov evaluates the covariance at the given lags", {
  ## Regression test: the covariance used to be evaluated at h[1] - h.
  s <- seq(from = 0, to = 1, length.out = 51)
  kappa <- 10
  op <- matern.rational(
    loc = s, nu = 0.8, range = sqrt(8 * 0.8) / kappa, sigma = 1, m = 2,
    parameterization = "matern"
  )
  true_cov <- function(ind) {
    matern.covariance(abs(s - s[ind]), kappa = kappa, sigma = 1, nu = 0.8)
  }
  expect_equal(as.vector(op$covariance(ind = 1)), true_cov(1), tolerance = 1e-2)
  expect_equal(as.vector(op$covariance(ind = 25)), true_cov(25),
    tolerance = 1e-2
  )
  full <- op$covariance()
  expect_equal(dim(full), c(length(s), length(s)))
  expect_equal(as.vector(full[, 25]), true_cov(25), tolerance = 1e-2)
})

test_that("the weighted-L2 error decreases with the order", {
  ## Without continuation in m, a plain multistart settles at the optimum of
  ## the next smaller order beyond m = 4, with one residue exactly zero.
  for (alpha in c(0.8, 1.2, 1.7)) {
    errs <- sapply(2:6, function(m) {
      cf <- rational.coefficients.wl2(alpha, d = 1, m = m)
      expect_true(min(cf$r) > 0)
      cf$rel_err
    })
    expect_true(all(diff(errs) < 0))
  }
  ## m = 6 collapses to the m = 5 fit without the continuation
  collapsed <- rational.coefficients.wl2(1.7, d = 1, m = 6, continuation = FALSE)
  full <- rational.coefficients.wl2(1.7, d = 1, m = 6)
  expect_lt(full$rel_err, collapsed$rel_err)
})

test_that("the NNLS solver is scale invariant and solves the problem", {
  set.seed(1)
  A <- matrix(abs(rnorm(200 * 5)), 200, 5)
  coef <- c(1, 2, 0.5, 3, 0.1)
  b <- as.vector(A %*% coef)
  expect_equal(rSPDE:::wl2_nnls(A, b), coef, tolerance = 1e-8)
  ## an absolute tolerance built from A alone would stop early for small b
  for (s in c(1e-8, 1e8)) {
    expect_equal(rSPDE:::wl2_nnls(A, s * b) / s, coef, tolerance = 1e-8)
  }
  ## the basis of an m = 6 fit spans six orders of magnitude in the poles; the
  ## solution must not flip when the right hand side is perturbed by a rounding
  ## error
  cells <- rSPDE:::wl2_weyl_cells(1, 1e-26, 900)
  theta <- c(3.75, log(c(1.95, 7.6, 42, 324, 4452, 712054)))
  Ab <- rSPDE:::wl2_basis("shifted", cells$x, theta) * sqrt(cells$w)
  f <- cells$x^1.11 * sqrt(cells$w)
  ## the true optimum, by enumerating the active sets
  brute <- function(A, b) {
    n <- ncol(A)
    best <- NULL
    best_res <- Inf
    for (k in 0:(2^n - 1)) {
      keep <- as.logical(bitwAnd(k, 2^(0:(n - 1))))
      x <- numeric(n)
      if (any(keep)) {
        ## a tight rank tolerance is essential: the columns are strongly
        ## correlated and the default 1e-7 declares them aliased
        cf <- tryCatch(qr.solve(A[, keep, drop = FALSE], b, tol = 1e-14),
          error = function(e) NULL
        )
        if (is.null(cf) || any(cf < 0)) next
        x[keep] <- cf
      }
      res <- sum((as.vector(A %*% x) - b)^2)
      if (res < best_res) {
        best_res <- res
        best <- x
      }
    }
    best
  }
  sse <- function(A, x, b) sum((as.vector(A %*% x) - b)^2)
  r1 <- rSPDE:::wl2_nnls(Ab, f)
  expect_equal(r1, brute(Ab, f), tolerance = 1e-6)
  ## and the solution must not flip when f is perturbed by a rounding error
  r2 <- rSPDE:::wl2_nnls(Ab, f * (1 + 1e-15))
  expect_equal(r1, r2 / (1 + 1e-15), tolerance = 1e-6)

  ## the same on bases drawn at random from the region the optimiser explores
  set.seed(3)
  for (trial in 1:10) {
    m <- sample(2:6, 1)
    shifted <- sample(c(TRUE, FALSE), 1)
    th <- sort(runif(m, 0, 2.2 * m + 6))
    if (shifted) th <- c(runif(1, -2, 4), th)
    A <- rSPDE:::wl2_basis(if (shifted) "shifted" else "plain", cells$x, th) *
      sqrt(cells$w)
    b <- cells$x^runif(1, 0.6, 1.9) * sqrt(cells$w)
    opt <- brute(A, b)
    expect_lte(
      sse(A, rSPDE:::wl2_nnls(A, b), b),
      sse(A, opt, b) * (1 + 1e-8) + 1e-300
    )
  }
})

test_that("the shifted class satisfies the stationarity diagnostics", {
  ## A genuine stationary point of the shifted-factor problem has all residues 
  ## positive, all poles negative, a shift in (0, 1), a vanishing derivative in 
  ## p0, and exactly 2m + 1 sign changes of the error on (0, 1) -- 2m for the 
  ## plain class.
  diagnose <- function(cf, alpha, d, m, x_min = NULL) {
    x_min <- if (is.null(x_min)) 1e-26 else x_min
    n_grid <- if (identical(x_min, 1e-26)) 900 else 300
    cells <- rSPDE:::wl2_weyl_cells(d, x_min, n_grid)
    e <- rSPDE:::wl2_symbol(cf, cells$x) - cells$x^alpha
    p0 <- if (is.null(cf$p0)) 0 else cf$p0
    wt <- cells$x^2 / (1 - p0 * cells$x)^2
    xs <- exp(seq(log(1e-24), log(1 - 1e-12), length.out = 50001))
    es <- rSPDE:::wl2_symbol(cf, xs) - xs^alpha
    list(
      stationarity = sum(cells$w * e * wt) / sum(cells$w * abs(e) * wt),
      sign_changes = sum(diff(sign(es)) != 0)
    )
  }

  ## the reference values of the note
  cases <- list(
    list(alpha = 1.5, d = 2, m = 2, x_min = NULL, err = 2.26e-3, b0 = 0.7951),
    list(alpha = 1.5, d = 2, m = 4, x_min = NULL, err = 8.55e-5, b0 = 0.8935),
    list(alpha = 1.25, d = 1, m = 3, x_min = NULL, err = 2.10e-4, b0 = 0.9160),
    list(alpha = 1.5, d = 2, m = 4, x_min = 1 / 801, err = 3.27e-5, b0 = 0.9047)
  )
  for (case in cases) {
    cf <- rational.coefficients.wl2(case$alpha, case$d, case$m,
      x_min = case$x_min
    )
    expect_equal(cf$rel_err, case$err, tolerance = 0.02)
    expect_equal(1 - cf$p0, case$b0, tolerance = 1e-3)
    ## the feasible set
    expect_true(all(cf$r > 0))
    expect_true(all(cf$p < 0))
    expect_true(cf$p0 > 0 && cf$p0 < 1)
    d9 <- diagnose(cf, case$alpha, case$d, case$m, case$x_min)
    expect_lt(abs(d9$stationarity), 1e-6)
    expect_equal(d9$sign_changes, 2 * case$m + 1)
  }

  ## the sign-change count over the whole class, both integer parts. The
  ## stationarity condition is a derivative in p0, so it only applies to the
  ## shifted class; for floor(alpha) = 0 there is no shift.
  for (alpha in c(0.7, 0.9, 1.2, 1.8)) {
    for (m in c(1, 2, 4)) {
      cf <- rational.coefficients.wl2(alpha, d = 1, m = m)
      d9 <- diagnose(cf, alpha, 1, m)
      expect_equal(d9$sign_changes, 2 * m + (alpha >= 1))
      if (!is.null(cf$p0)) {
        expect_lt(abs(d9$stationarity), 1e-6)
      }
    }
  }
})

test_that("the shift improves on the unshifted class for floor(alpha) = 1", {
  ## p0 = 0 is the plain class with the unshifted integer factor, i.e. the form
  ## the tabulated coefficients use. The note expects a gain of 1.6-2.8.
  fit_unshifted <- function(alpha, d, m) {
    cells <- rSPDE:::wl2_weyl_cells(d, 1e-26, 900)
    x <- cells$x
    sw <- sqrt(cells$w)
    f <- x^alpha * sw
    ## no analytic Jacobian for this ad hoc class, so wl2_lm falls back to
    ## finite differences
    res <- function(th, need_jac = FALSE) {
      A <- x * rSPDE:::wl2_basis("plain", x, th) * sw
      list(r = drop(A %*% rSPDE:::wl2_nnls(A, f)) - f)
    }
    best <- NULL
    for (t0 in rSPDE:::wl2_starts(m, 1e-26, 8L, "plain", NULL, 1L)) {
      o <- tryCatch(rSPDE:::wl2_lm(res, t0), error = function(e) NULL)
      if (!is.null(o) && (is.null(best) || o$cost < best$cost)) best <- o
    }
    sqrt(2 * best$cost / sum(cells$w * x^(2 * alpha)))
  }
  for (alpha in c(1.1, 1.5)) {
    for (m in c(1, 2, 4)) {
      shifted <- rational.coefficients.wl2(alpha, d = 1, m = m)
      expect_gt(fit_unshifted(alpha, 1, m) / shifted$rel_err, 1.5)
      expect_true(1 - shifted$p0 > 0.6 && 1 - shifted$p0 < 0.99)
    }
  }
})

test_that("the analytic Jacobian matches finite differences", {
  ## Golub-Pereyra derivative of the variable-projection residual, checked at a
  ## perturbed optimum so that the passive set is the realistic one.
  set.seed(1)
  for (kind in c("plain", "shifted", "shifted2")) {
    for (d in 1:2) {
      for (m in c(2, 4, 6)) {
        alpha <- switch(kind, plain = 0.8, shifted = 1.4, shifted2 = 2.4)
        cells <- rSPDE:::wl2_weyl_cells(d, 1e-26, 900)
        sw <- sqrt(cells$w)
        f <- cells$x^alpha * sw
        theta <- rational.coefficients.wl2(alpha, d, m)$theta +
          0.05 * rnorm(m + (rSPDE:::wl2_q(kind) > 0))
        a <- rSPDE:::wl2_resjac(kind, cells$x, sw, f, theta)
        expect_true(all(a$c > 0))
        Jn <- matrix(0, length(cells$x), length(theta))
        for (k in seq_along(theta)) {
          h <- 1e-7 * max(abs(theta[k]), 1)
          tk <- theta
          tk[k] <- tk[k] + h
          Jn[, k] <- (rSPDE:::wl2_resjac(
            kind, cells$x, sw, f, tk, FALSE
          )$r - a$r) / h
        }
        expect_lt(max(abs(a$J - Jn)) / max(abs(Jn)), 1e-4)
      }
    }
  }
})

test_that("the compiled solver agrees with the R implementation", {
  ## src/wl2_fit.cpp mirrors wl2_nnls(), wl2_resjac() and wl2_lm(); it is used
  ## whenever it is available, with the R code as the reference and fallback.
  skip_if_not(rSPDE:::wl2_have_cpp())
  old <- options(rSPDE.wl2.use.cpp = TRUE)
  on.exit(options(old), add = TRUE)
  cases <- list(
    list(alpha = 0.5, d = 1, m = 2, x_min = NULL, type = "covariance"),
    list(alpha = 0.75, d = 2, m = 4, x_min = NULL, type = "covariance"),
    list(alpha = 1.5, d = 2, m = 2, x_min = NULL, type = "covariance"),
    list(alpha = 1.25, d = 1, m = 3, x_min = NULL, type = "covariance"),
    list(alpha = 1.5, d = 2, m = 4, x_min = 1 / 801, type = "covariance"),
    list(alpha = 2.4, d = 1, m = 2, x_min = NULL, type = "covariance"),
    list(alpha = 2.8, d = 2, m = 4, x_min = NULL, type = "covariance"),
    list(alpha = 2.4, d = 2, m = 3, x_min = 1 / 801, type = "covariance"),
    list(alpha = 1.5, d = 2, m = 3, x_min = 1 / 801, type = "operator")
  )
  for (cs in cases) {
    options(rSPDE.wl2.use.cpp = TRUE)
    cf <- do.call(rational.coefficients.wl2, cs)
    options(rSPDE.wl2.use.cpp = FALSE)
    cr <- do.call(rational.coefficients.wl2, cs)
    ## The objective is flat near the optimum, so the coefficients may differ
    ## by more than the error they produce; require both.
    expect_equal(cf$rel_err, cr$rel_err, tolerance = 1e-8)
    expect_equal(cf$p, cr$p, tolerance = 1e-5)
    expect_equal(cf$r, cr$r, tolerance = 1e-5)
    if (!is.null(cr$p0)) {
      expect_equal(cf$p0, cr$p0, tolerance = 1e-5)
    }
  }
})

test_that("the non-negative least squares in C matches the R version", {
  skip_if_not(rSPDE:::wl2_have_cpp())
  ## Exercised indirectly above; here directly, on the badly scaled matrices
  ## the fit actually produces, through a one-iteration fit from a fixed start.
  cells <- rSPDE:::wl2_weyl_cells(2, 1e-26, 900)
  sw <- sqrt(cells$w)
  f <- cells$x^1.4 * sw
  theta <- c(1.2, 0.5, 2, 5, 9)
  a <- rSPDE:::wl2_resjac("shifted", cells$x, sw, f, theta, FALSE)
  expect_true(all(a$c >= 0))
  A <- rSPDE:::wl2_basis("shifted", cells$x, theta) * sw
  expect_equal(drop(A %*% a$c) - f, a$r, tolerance = 1e-12)
})

test_that("the shipped mesh-free tables cover the expected configurations", {
  ## Built by data-raw/wl2_tables.R and stored in R/sysdata.rda, because they
  ## depend on neither the mesh nor kappa and are what matern.operators() uses
  ## unless a spectral interval is asked for.
  tabs <- get0("wl2_meshfree_tables",
    envir = asNamespace("rSPDE"), ifnotfound = NULL
  )
  expect_false(is.null(tabs))
  for (d in 1:3) {
    for (m_alpha in 0:2) {
      for (m in 1:6) {
        n_alpha <- length(rSPDE:::wl2_alpha_grid(m_alpha, d, "covariance", 0.01))
        tab <- rSPDE:::wl2_shipped_table(d, m, m_alpha, 300)
        if (n_alpha == 0) {
          ## alpha = nu + d / 2 with nu > 0, and the fit needs a finite trace.
          expect_null(tab)
          next
        }
        expect_false(is.null(tab))
        expect_equal(nrow(tab), n_alpha)
        expect_true(all(c(
          "alpha", paste0("r", seq_len(m)), paste0("p", seq_len(m)),
          "k", "rel_err"
        ) %in% names(tab)))
        expect_equal("p0" %in% names(tab), m_alpha > 0)
        ## A valid model at every alpha: non-negative residues, poles below one.
        expect_true(all(tab[, paste0("r", seq_len(m))] >= 0))
        expect_true(all(tab[, paste0("p", seq_len(m))] < 1))
        expect_true(all(is.finite(tab$rel_err)))
        expect_equal(attr(tab, "x_min"), NULL)
        expect_equal(
          attr(tab, "kind"),
          c("plain", "shifted", "shifted2")[m_alpha + 1]
        )
        expect_false(is.null(attr(tab, "quad")))
      }
    }
  }
})

test_that("the shipped tables are what the fits produce", {
  ## Guards against the stored tables going stale: anything that changes the
  ## coefficients (wl2_fit, wl2_resjac, wl2_nnls, wl2_lm, wl2_weyl_cells, the
  ## sweep in wl2_coefficient_table, or wl2_tol) means data-raw/wl2_tables.R
  ## has to be run again. Skipped on CRAN: rebuilding takes half a minute, and
  ## the last digits of a fit depend on the platform's BLAS, while the stored
  ## table is the same everywhere.
  skip_on_cran()
  for (cs in list(c(1, 1, 0), c(2, 2, 1), c(3, 1, 1), c(1, 3, 2), c(2, 2, 2))) {
    d <- cs[1]
    m <- cs[2]
    m_alpha <- cs[3]
    rSPDE:::wl2_clear_cache()
    fresh <- rSPDE:::wl2_coefficient_table(
      d = d, m = m, m_alpha = m_alpha, type = "covariance",
      cache = FALSE, shipped = FALSE
    )
    stored <- rSPDE:::wl2_shipped_table(d, m, m_alpha, 300)
    expect_equal(stored$alpha, fresh$alpha)
    ## The error is the thing that matters and is insensitive to the last
    ## digits of the poles; the coefficients themselves are checked loosely,
    ## since the objective is flat at its minimum.
    expect_equal(stored$rel_err, fresh$rel_err, tolerance = 1e-6)
    for (j in seq_len(m)) {
      expect_equal(stored[[paste0("r", j)]], fresh[[paste0("r", j)]],
        tolerance = 1e-4
      )
      expect_equal(stored[[paste0("p", j)]], fresh[[paste0("p", j)]],
        tolerance = 1e-4
      )
    }
    if (m_alpha == 1) {
      expect_equal(stored$p0, fresh$p0, tolerance = 1e-4)
    }
  }
})

test_that("a model built from a shipped table equals one built from the fits", {
  ## The lookup must be transparent: same coefficients, same covariance.
  skip_on_cran()
  d <- 2
  m <- 2
  nu <- 0.6
  rSPDE:::wl2_clear_cache()
  built <- rSPDE:::wl2_coefficient_table(
    d = d, m = m, m_alpha = 1, type = "covariance",
    cache = FALSE, shipped = FALSE
  )
  stored <- rSPDE:::wl2_coefficient_table(
    d = d, m = m, m_alpha = 1, type = "covariance", cache = FALSE
  )
  a <- nu + d / 2
  cb <- rSPDE:::wl2_interp_coefficients(built, a)
  cs <- rSPDE:::wl2_interp_coefficients(stored, a)
  expect_equal(cs$r, cb$r, tolerance = 1e-6)
  expect_equal(cs$p, cb$p, tolerance = 1e-6)
  expect_equal(cs$p0, cb$p0, tolerance = 1e-6)
  x <- exp(seq(log(1e-12), 0, length.out = 500))
  expect_equal(rSPDE:::wl2_symbol(cs, x), rSPDE:::wl2_symbol(cb, x),
    tolerance = 1e-6
  )
})

## The columns of a table, without the attributes that record how it was
## asked for rather than what it contains.
wl2_cols <- function(tab) lapply(as.list(tab), as.numeric)

test_that("the persistent cache is off unless it is asked for", {
  ## A package must not write outside the session temporary directory without
  ## being told to.
  old <- options(rSPDE.cache.dir = NULL)
  on.exit(options(old), add = TRUE)
  withr::with_envvar(c(RSPDE_CACHE_DIR = NA), {
    expect_null(rspde.cache())
    expect_null(rSPDE:::wl2_cache_dir())
    expect_null(rSPDE:::wl2_cache_paths("covariance", 2, 4, 1, 1 / 801, 8, 300, 0.01))
    expect_null(rSPDE:::wl2_cache_read(
      "covariance", 2, 4, 1, 1 / 801, 8, 300, 0.01, 300
    ))
  })
  withr::with_envvar(c(RSPDE_CACHE_DIR = file.path(tempdir(), "env-cache")), {
    expect_equal(rSPDE:::wl2_cache_dir(), file.path(tempdir(), "env-cache"))
  })
})

test_that("x_min round-trips through a file name", {
  for (x in c(1e-26, 1 / 801, 0.5, 1.97428e-4, 0.999)) {
    enc <- rSPDE:::wl2_xmin_encode(x)
    expect_false(grepl("[^A-Za-z0-9]", enc))
    expect_equal(rSPDE:::wl2_xmin_decode(enc), x, tolerance = 1e-10)
  }
  expect_equal(rSPDE:::wl2_xmin_encode(NULL), "free")
  expect_null(rSPDE:::wl2_xmin_decode("free"))
})

test_that("a cached table is found again, and only when it is valid", {
  skip_on_cran()
  dir <- file.path(tempdir(), "rspde-cache-test")
  unlink(dir, recursive = TRUE)
  old <- options(rSPDE.cache.dir = NULL)
  on.exit(
    {
      options(old)
      unlink(dir, recursive = TRUE)
    },
    add = TRUE
  )
  rspde.cache(dir)
  x0 <- 1 / 801
  build <- function(x) {
    rSPDE:::wl2_clear_cache()
    rSPDE:::wl2_coefficient_table(
      d = 2, m = 2, m_alpha = 1, x_min = x, type = "covariance"
    )
  }
  first <- build(x0)
  again <- build(x0)
  expect_equal(wl2_cols(again), wl2_cols(first))
  expect_equal(attr(again, "x_min"), x0)

  ## A table fitted on a slightly wider interval still covers the whole
  ## spectrum, so it is reused, and says which interval it was fitted on.
  near <- build(x0 * 1.003)
  expect_equal(wl2_cols(near), wl2_cols(first))
  expect_equal(attr(near, "x_min"), x0)

  ## One fitted on a narrower interval would leave part of the spectrum
  ## unfitted, and must never be reused.
  sharper <- build(x0 * 0.97)
  expect_equal(attr(sharper, "x_min"), x0 * 0.97)

  ## Far enough away it is worth refitting rather than losing the accuracy.
  far <- build(x0 * 3)
  expect_equal(attr(far, "x_min"), x0 * 3)

  expect_gt(length(list.files(dir, recursive = TRUE)), 0)
  rspde.cache(clear = TRUE)
  expect_equal(length(list.files(dir, recursive = TRUE)), 0)
})

test_that("rspde.wl2.table gives the model the kappa_ref fit", {
  skip_on_cran()
  old <- options(rSPDE.cache.dir = NULL)
  on.exit(options(old), add = TRUE)
  mesh <- fmesher::fm_mesh_1d(seq(0, 1, length.out = 101))
  nu <- 0.8
  kappa_ref <- 5
  tab <- rspde.wl2.table(
    m = 2, d = 1, nu = nu, kappa_ref = kappa_ref, mesh = mesh
  )
  expect_equal(
    attr(tab, "x_min"),
    rspde.xmin(mesh = mesh, kappa_ref = kappa_ref, nu = nu)
  )
  expect_equal(attr(tab, "d"), 1)
  expect_equal(attr(tab, "m"), 2)
  expect_equal(attr(tab, "type"), "covariance")

  args <- list(
    nu = nu, range = 0.3, sigma = 1, mesh = mesh, m = 2, d = 1,
    parameterization = "matern", type = "covariance",
    type_rational_approximation = "wl2"
  )
  supplied <- do.call(matern.operators, c(args, list(wl2_table = tab)))
  computed <- do.call(matern.operators, c(args, list(kappa_ref = kappa_ref)))
  meshfree <- do.call(matern.operators, args)
  expect_equal(as.matrix(supplied$Q), as.matrix(computed$Q))
  ## and it is not simply the default table
  expect_false(isTRUE(all.equal(
    as.matrix(supplied$Q), as.matrix(meshfree$Q)
  )))

  ## x_min given directly needs no mesh at all. Only the coefficients are
  ## compared: a table built from kappa_ref also records it.
  direct <- rspde.wl2.table(m = 2, d = 1, nu = nu, x_min = attr(tab, "x_min"))
  expect_equal(wl2_cols(direct), wl2_cols(tab))
  expect_equal(attr(direct, "x_min"), attr(tab, "x_min"))
})

test_that("a table for the wrong model is refused", {
  skip_on_cran()
  mesh <- fmesher::fm_mesh_1d(seq(0, 1, length.out = 51))
  args <- list(
    nu = 0.8, range = 0.3, sigma = 1, mesh = mesh, m = 2, d = 1,
    parameterization = "matern", type = "covariance",
    type_rational_approximation = "wl2"
  )
  expect_error(
    do.call(matern.operators, c(args, list(
      wl2_table = rspde.wl2.table(m = 3, d = 1, nu = 0.8, x_min = 1e-3)
    ))),
    "m = 3"
  )
  expect_error(
    do.call(matern.operators, c(args, list(
      wl2_table = rspde.wl2.table(m = 2, d = 2, nu = 0.8, x_min = 1e-3)
    ))),
    "d = 2"
  )
  expect_error(
    do.call(matern.operators, c(args, list(wl2_table = data.frame(alpha = 1)))),
    "rspde.wl2.table"
  )
  expect_error(rspde.wl2.table(m = 2, d = 1, nu = 0.8), "x_min")
  expect_error(rspde.wl2.table(m = 2, d = 1), "alpha or nu")
})

test_that("the data weight x^s is the one it claims to be", {
  ## wl2_weyl_cells returns cell integrals of w_d(x) x^s, taken in log x, so
  ## the integrand there is x^(s - d/2) (1-x)^(d/2-1).
  for (d in 1:2) {
    for (s in c(0, 0.5, 1)) {
      n <- 200
      cl <- rSPDE:::wl2_weyl_cells(d, 1e-3, n, s = s)
      lx <- seq(log(1e-3), 0, length.out = n)
      ed <- c(lx[1], (lx[-1] + lx[-n]) / 2, 0)
      ## The last cells are dropped: for d = 1 the density has a (1-x)^(-1/2)
      ## singularity at x = 1, where the midpoint rule of wl2_weyl_cells and
      ## integrate() differ by more than the identity being checked.
      keep <- seq_len(n - 3)
      ref <- vapply(keep, function(i) {
        stats::integrate(function(u) {
          xx <- exp(u)
          xx^(s - d / 2) * (1 - xx)^(d / 2 - 1)
        }, ed[i], ed[i + 1], rel.tol = 1e-10)$value
      }, numeric(1))
      expect_equal(cl$w[keep], ref, tolerance = 1e-6)
    }
  }
  ## s = 0 is the default, and changes nothing
  expect_equal(
    rSPDE:::wl2_weyl_cells(2, 1e-4, 100),
    rSPDE:::wl2_weyl_cells(2, 1e-4, 100, s = 0)
  )
  expect_equal(
    rational.coefficients.wl2(1.5, 2, 3),
    rational.coefficients.wl2(1.5, 2, 3, s = 0, k_term = FALSE)
  )
})

test_that("the fit minimises the weighted objective it is given", {
  ## Checked against an independent quadrature ten times finer than the fit's.
  obj <- function(cf, alpha, d, x_min, s) {
    cl <- rSPDE:::wl2_weyl_cells(d, x_min, 3000, s = s)
    e <- rSPDE:::wl2_symbol(cf, cl$x) - cl$x^alpha
    sqrt(sum(cl$w * e^2) / sum(cl$w * cl$x^(2 * alpha)))
  }
  for (cs in list(
    list(a = 0.75, d = 2, m = 3, x = 1e-3, s = 0),
    list(a = 0.75, d = 2, m = 3, x = 1e-3, s = 1),
    list(a = 0.6, d = 2, m = 2, x = 1e-4, s = 0.5)
  )) {
    for (kt in c(FALSE, TRUE)) {
      cf <- rational.coefficients.wl2(
        cs$a, cs$d, cs$m, x_min = cs$x, s = cs$s, k_term = kt
      )
      expect_equal(cf$rel_err, obj(cf, cs$a, cs$d, cs$x, cs$s),
        tolerance = 5e-3
      )
      expect_equal(cf$k > 0, kt)
      expect_true(all(cf$r > 0))
      expect_true(all(cf$p < 1))
    }
  }
})

test_that("the constant term helps, and leaves on its own when it may not stay", {
  ## On (0, 1] with s = 0 the integral of k^2 w_d diverges, so the optimal
  ## constant is zero; on a finite interval it is not, and it is worth having.
  gain <- function(x_min) {
    a <- rational.coefficients.wl2(0.75, 2, 3, x_min = x_min, k_term = FALSE)
    b <- rational.coefficients.wl2(0.75, 2, 3, x_min = x_min, k_term = TRUE)
    c(gain = a$rel_err / b$rel_err, k = b$k)
  }
  wide <- gain(1e-2)
  mid <- gain(1e-4)
  free <- gain(1e-26)
  expect_gt(wide[["gain"]], 2)
  expect_gt(mid[["gain"]], 1.5)
  expect_equal(free[["gain"]], 1, tolerance = 1e-6)
  ## and the constant shrinks with the interval
  expect_gt(wide[["k"]], mid[["k"]])
  expect_lt(free[["k"]], 1e-12)
})

test_that("a constant term is refused where the model cannot represent it", {
  ## The operator factorisation has no constant term, so coefficients fitted
  ## with one must not be silently turned into an operator-based model.
  cf <- rational.coefficients.wl2(0.75, 2, 2, x_min = 1e-3, k_term = TRUE)
  expect_gt(cf$k, 0)
  expect_error(rSPDE:::wl2_roots(cf), "constant term")
  expect_silent(rSPDE:::wl2_roots(
    rational.coefficients.wl2(1.5, 2, 2, x_min = 1e-3, type = "operator")
  ))
})

test_that("the objective must be finite", {
  ## Near zero the integrand is e(x)^2 x^(s - 1 - d/2) with e ~ x^alpha, so
  ## 2 alpha + s > d / 2 is needed. The weight can make a fit possible that is
  ## not possible without it.
  expect_error(rational.coefficients.wl2(0.4, 2, 2, x_min = 1e-3), "finite")
  expect_silent(rational.coefficients.wl2(0.4, 2, 2, x_min = 1e-3, s = 0.5))
  expect_error(rational.coefficients.wl2(0.75, 2, 2, s = "a"), "single finite")
})

test_that("the compiled solver handles the constant term too", {
  skip_if_not(rSPDE:::wl2_have_cpp())
  old <- options(rSPDE.wl2.use.cpp = TRUE)
  on.exit(options(old), add = TRUE)
  for (cs in list(
    list(a = 0.75, d = 2, m = 3, x = 1e-3, s = 0),
    list(a = 0.75, d = 2, m = 4, x = 1e-3, s = 1),
    list(a = 0.9, d = 2, m = 3, x = 1e-5, s = 1)
  )) {
    options(rSPDE.wl2.use.cpp = TRUE)
    cc <- rational.coefficients.wl2(
      cs$a, cs$d, cs$m, x_min = cs$x, s = cs$s, k_term = TRUE
    )
    options(rSPDE.wl2.use.cpp = FALSE)
    rr <- rational.coefficients.wl2(
      cs$a, cs$d, cs$m, x_min = cs$x, s = cs$s, k_term = TRUE
    )
    expect_equal(cc$rel_err, rr$rel_err, tolerance = 1e-8)
    expect_equal(cc$k, rr$k, tolerance = 1e-5)
    expect_equal(cc$p, rr$p, tolerance = 1e-5)
  }
})

test_that("a table can carry the weight and the constant term", {
  ## The options of rational.coefficients.wl2() have to be reachable at the
  ## table level too, otherwise they can only be used one alpha at a time.
  old <- options(rSPDE.cache.dir = NULL)
  on.exit(options(old), add = TRUE)
  for (cs in list(
    list(s = 0, k = FALSE), list(s = 0, k = TRUE), list(s = 1, k = TRUE)
  )) {
    rSPDE:::wl2_clear_cache()
    tab <- rSPDE:::wl2_coefficient_table(
      d = 1, m = 2, m_alpha = 0, x_min = 1e-3, type = "covariance",
      s = cs$s, k_term = cs$k, by = 0.05
    )
    expect_equal(attr(tab, "s"), cs$s)
    expect_equal(attr(tab, "k_term"), cs$k)
    expect_equal(any(tab$k > 0), cs$k)
    ## the interpolation must carry the constant term through and stay as good
    ## as a direct fit at the same alpha
    cf <- rSPDE:::wl2_interp_coefficients(tab, 0.775)
    direct <- rational.coefficients.wl2(
      0.775, 1, 2, x_min = 1e-3, s = cs$s, k_term = cs$k
    )
    expect_equal(cf$k > 0, cs$k)
    expect_equal(cf$rel_err, direct$rel_err, tolerance = 1e-3)
  }
})

test_that("stored tables are never reused for a different weight", {
  old <- options(rSPDE.cache.dir = NULL)
  on.exit(options(old), add = TRUE)
  rSPDE:::wl2_clear_cache()
  a <- rSPDE:::wl2_coefficient_table(d = 2, m = 2, m_alpha = 1,
                                     type = "covariance")
  b <- rSPDE:::wl2_coefficient_table(d = 2, m = 2, m_alpha = 1,
                                     type = "covariance", s = 1)
  expect_equal(attr(a, "s"), 0)
  expect_equal(attr(b, "s"), 1)
  expect_false(isTRUE(all.equal(a$rel_err, b$rel_err)))
  ## and the cache paths differ
  p0 <- rSPDE:::wl2_cache_paths("covariance", 2, 2, 1, 1e-3, 8, 300, 0.01)
  p1 <- rSPDE:::wl2_cache_paths("covariance", 2, 2, 1, 1e-3, 8, 300, 0.01,
                                s = 1, k_term = TRUE)
  expect_null(p0)
  expect_null(p1)
  withr::with_options(
    list(rSPDE.cache.dir = file.path(tempdir(), "wcache")), {
      q0 <- rSPDE:::wl2_cache_paths("covariance", 2, 2, 1, 1e-3, 8, 300, 0.01)
      q1 <- rSPDE:::wl2_cache_paths("covariance", 2, 2, 1, 1e-3, 8, 300, 0.01,
                                    s = 1, k_term = TRUE)
      expect_false(identical(q0$dir, q1$dir))
    }
  )
})

test_that("the largest eigenvalue is computed without RSpectra", {
  ## rspde.xmin(eigenvalue = "exact") used RSpectra::eigs_sym(), which fails to
  ## converge on exactly the meshes it is wanted for -- a 1d mesh with 2001
  ## nodes returned nothing at all, and the zero-length result propagated
  ## silently into x_min. Lanczos handles it in a hundredth of a second.
  for (cfg in list(
    list(mesh = fmesher::fm_mesh_1d(seq(0, 1, length.out = 401))),
    list(mesh = fmesher::fm_mesh_1d(seq(0, 1, length.out = 2001)))
  )) {
    fe <- fmesher::fm_fem(cfg$mesh)
    cl <- Matrix::rowSums(fe$c0)
    Ci2 <- Matrix::Diagonal(length(cl), 1 / sqrt(cl))
    S <- Ci2 %*% fe$g1 %*% Ci2
    got <- rSPDE:::wl2_lambda_max(S)
    expect_length(got, 1)
    expect_true(is.finite(got))
    ref <- max(eigen(as.matrix((S + Matrix::t(S)) / 2),
      symmetric = TRUE, only.values = TRUE
    )$values)
    ## converges from below, and slowly in one dimension where the top of the
    ## spectrum is clustered: 1.6e-3 at the default twenty steps
    expect_lt(abs(got / ref - 1), 5e-3)
    expect_lte(got, ref * (1 + 1e-10))
  }
  ## and through the user-facing function, where the result must be a single
  ## number in (0, 1) and at most the Gershgorin bound's interval
  for (n in c(401, 2001)) {
    mesh <- fmesher::fm_mesh_1d(seq(0, 1, length.out = n))
    xe <- rspde.xmin(mesh = mesh, kappa_ref = 5, nu = 0.5, eigenvalue = "exact")
    xb <- rspde.xmin(mesh = mesh, kappa_ref = 5, nu = 0.5)
    expect_length(xe, 1)
    expect_true(is.finite(xe) && xe > 0 && xe < 1)
    ## the bound over-estimates the largest eigenvalue, so it under-estimates
    ## x_min; the exact value is never smaller, and never larger than the truth
    ## allows, since rspde.xmin caps it at the bound after a 1% margin
    expect_gte(xe, xb * (1 - 1e-8))
  }
  ## a small matrix falls through to a dense decomposition
  S <- Matrix::Diagonal(5, c(1, 2, 3, 4, 5))
  expect_equal(rSPDE:::wl2_lambda_max(S), 5)
})

test_that("the intrinsic scaling is computed without RSpectra", {
  ## The scaling of the intrinsic models is the smallest non-zero eigenvalue of
  ## C^{-1/2} G C^{-1/2}. Its null vector is known rather than computed, since
  ## G is a stiffness matrix and G 1 = 0, so deflating it turns the second
  ## smallest eigenvalue into the smallest and shift-and-invert Lanczos reaches
  ## it in about a dozen steps.
  skip_if_not_installed("fmesher")
  set.seed(1)
  mesh <- fmesher::fm_mesh_2d(
    loc = cbind(runif(150), runif(150)), max.edge = c(0.08, 0.24)
  )
  fe <- fmesher::fm_fem(mesh)
  C <- fe$c0
  G <- fe$g1
  n <- dim(C)[1]
  ## the assumption the deflation rests on
  expect_lt(max(abs(as.vector(G %*% rep(1, n)))), 1e-10 * max(abs(Matrix::diag(G))))
  Cd <- Matrix::Diagonal(n, 1 / sqrt(Matrix::diag(C)))
  Gg <- as(Cd %*% G %*% Cd, "CsparseMatrix")
  got <- rSPDE:::rspde_lambda_min_nonzero(Gg, sqrt(Matrix::diag(C)))
  ref <- sort(eigen(as.matrix((Gg + Matrix::t(Gg)) / 2),
    symmetric = TRUE, only.values = TRUE
  )$values)[2]
  expect_equal(got, ref, tolerance = 1e-8)
  ## and through the model, where it is the scaling
  op <- intrinsic.operators(C = C, G = G, beta = 1, tau = 1, d = 2, m = 2)
  expect_equal(op$scaling, ref, tolerance = 1e-8)
  ## a supplied scaling still wins, and opts is honoured rather than ignored
  expect_equal(
    intrinsic.operators(C = C, G = G, beta = 1, tau = 1, d = 2, m = 2,
                        scaling = 4.6)$scaling, 4.6
  )
  expect_equal(
    intrinsic.operators(C = C, G = G, beta = 1, tau = 1, d = 2, m = 2,
                        opts = list(tol = 1e-6, maxitr = 100))$scaling,
    ref, tolerance = 1e-5
  )
  ## a vector that is not a null vector is refused rather than silently used
  expect_true(is.na(rSPDE:::rspde_lambda_min_nonzero(Gg, seq_len(n))))
})
