context("hybrid.spde")

# ---------------------------------------------------------------------------
# Helper used throughout the tests: a simple 1d FEM discretisation.
# ---------------------------------------------------------------------------
make_simple_setup <- function(n = 101, kappa = 10, tau = 0.5,
                              alpha = 1.3, beta_X = 1.5,
                              X_field = function(s) sin(2 * pi * s)) {
  x <- seq(from = 0, to = 1, length.out = n)
  X <- matrix(X_field(x), ncol = 1)
  list(
    x = x,
    X = X,
    kappa = kappa,
    tau = tau,
    alpha = alpha,
    beta_X = beta_X
  )
}


test_that("hybrid.spde constructor stores X, beta_X and a usable mean", {
  s <- make_simple_setup()

  op <- hybrid.spde(
    X = s$X, beta_X = s$beta_X,
    kappa = s$kappa, tau = s$tau, alpha = s$alpha,
    loc_mesh = s$x, d = 1, parameterization = "spde",
    type = "covariance"
  )

  expect_s3_class(op, "hybrid_spde")
  expect_equal(as.numeric(op$beta_X), s$beta_X)
  expect_equal(dim(as.matrix(op$X)), c(length(s$x), 1))
  expect_equal(length(op$mu), length(s$x))
  expect_true(all(is.finite(op$mu)))

  # The mean must be zero whenever beta_X is zero.
  op0 <- hybrid.spde(
    X = s$X, beta_X = 0,
    kappa = s$kappa, tau = s$tau, alpha = s$alpha,
    loc_mesh = s$x, d = 1, parameterization = "spde",
    type = "covariance"
  )
  expect_equal(as.numeric(op0$mu), rep(0, length(s$x)))
})


test_that("hybrid.spde rejects mismatched dimensions and inputs", {
  s <- make_simple_setup(n = 51)
  X_wrong <- matrix(1, nrow = length(s$x) + 1, ncol = 1)

  expect_error(
    hybrid.spde(
      X = X_wrong, beta_X = s$beta_X,
      kappa = s$kappa, tau = s$tau, alpha = s$alpha,
      loc_mesh = s$x, d = 1, parameterization = "spde"
    ),
    regexp = "number of mesh nodes"
  )

  expect_error(
    hybrid.spde(
      X = s$X, beta_X = c(1, 2),
      kappa = s$kappa, tau = s$tau, alpha = s$alpha,
      loc_mesh = s$x, d = 1, parameterization = "spde"
    ),
    regexp = "length"
  )

  expect_error(
    hybrid.spde(
      X = NULL, beta_X = 1,
      kappa = s$kappa, tau = s$tau, alpha = s$alpha,
      loc_mesh = s$x, d = 1, parameterization = "spde"
    ),
    regexp = "X"
  )
})


test_that("Mean computed via the hybrid model approximates beta * L^{-alpha/2} X", {
  # The hybrid model is Y = beta L^{-alpha/2} X + tau^{-1} L^{-alpha/2} W,
  # so the mean does not depend on tau. Compare the mean computed by
  # hybrid.spde to a direct evaluation of L^{-alpha/2} X via an
  # independent fractional.operators object built with tau = 1.
  s <- make_simple_setup(n = 81, kappa = 7, tau = 0.8,
                         alpha = 1.4, beta_X = 2)
  op_h <- hybrid.spde(
    X = s$X, beta_X = s$beta_X,
    kappa = s$kappa, tau = s$tau, alpha = s$alpha,
    loc_mesh = s$x, d = 1, parameterization = "spde",
    type = "covariance"
  )

  mesh_1d <- fmesher::fm_mesh_1d(s$x)
  fem <- fmesher::fm_fem(mesh_1d)
  op_ref <- fractional.operators(
    L = fem$g1 + s$kappa^2 * fem$c0,
    beta = s$alpha / 2,
    C = fem$c0,
    scale.factor = s$kappa^2,
    m = op_h$m,
    tau = 1
  )

  # FEM-consistent reference: mu_h = Pr Pl^{-1} (C X) * beta, where
  # the C factor comes from the variational form of L mu = beta X.
  C_diag_ref <- Matrix::diag(op_ref$C)
  CX <- C_diag_ref * s$X
  LinvCX_ref <- as.vector(Pr.mult(op_ref, Pl.solve(op_ref, CX)))
  mu_ref <- s$beta_X * LinvCX_ref

  expect_equal(as.numeric(op_h$mu), as.numeric(mu_ref), tolerance = 1e-10)

  # Mean must NOT change when tau changes (the whole point of the new
  # parameterisation).
  op_h_tau <- hybrid.spde(
    X = s$X, beta_X = s$beta_X,
    kappa = s$kappa, tau = 3 * s$tau, alpha = s$alpha,
    loc_mesh = s$x, d = 1, parameterization = "spde",
    type = "covariance"
  )
  expect_equal(as.numeric(op_h$mu), as.numeric(op_h_tau$mu),
               tolerance = 1e-10)
})


test_that("Separate kappa_mu: constructor and update preserve mode and parameters", {
  set.seed(1)
  n_mesh <- 51
  x <- seq(0, 1, length.out = n_mesh)
  X <- matrix(sin(2 * pi * x), n_mesh, 1)

  # Linked (default)
  op_lk <- hybrid.spde(X = X, beta_X = 1.0,
                      kappa = 5, tau = 1, alpha = 2,
                      loc_mesh = x, d = 1, parameterization = "spde",
                      type = "covariance")
  expect_false(isTRUE(op_lk$separate_kappa_mu))
  expect_equal(as.numeric(op_lk$kappa_mu), as.numeric(op_lk$kappa))

  # Separate (kappa_mu = 10, kappa = 5)
  op_sep <- hybrid.spde(X = X, beta_X = 1.0,
                       kappa = 5, kappa_mu = 10, tau = 1, alpha = 2,
                       loc_mesh = x, d = 1, parameterization = "spde",
                       type = "covariance")
  expect_true(isTRUE(op_sep$separate_kappa_mu))
  expect_equal(op_sep$kappa, 5)
  expect_equal(op_sep$kappa_mu, 10)
  # The mu values should differ from the linked case at kappa = 5.
  expect_true(max(abs(op_sep$mu - op_lk$mu)) > 1e-6)

  # Updating kappa on linked: kappa_mu follows.
  op_lk2 <- update(op_lk, kappa = 8, check_stationarity = FALSE)
  expect_false(isTRUE(op_lk2$separate_kappa_mu))
  expect_equal(as.numeric(op_lk2$kappa_mu), as.numeric(op_lk2$kappa))
  expect_equal(as.numeric(op_lk2$kappa), 8)

  # Update kappa on separate: kappa changes but kappa_mu stays.
  op_sep2 <- update(op_sep, kappa = 7, check_stationarity = FALSE)
  expect_true(isTRUE(op_sep2$separate_kappa_mu))
  expect_equal(as.numeric(op_sep2$kappa), 7)
  expect_equal(as.numeric(op_sep2$kappa_mu), 10)

  # Update kappa_mu on linked: switches to separate.
  op_lk3 <- update(op_lk, kappa_mu = 12, check_stationarity = FALSE)
  expect_true(isTRUE(op_lk3$separate_kappa_mu))
  expect_equal(as.numeric(op_lk3$kappa), 5)
  expect_equal(as.numeric(op_lk3$kappa_mu), 12)
})


test_that("rspde_lme exposes kappa_mu as a parameter only when separate", {
  set.seed(11)
  n_mesh <- 31
  x <- seq(0, 1, length.out = n_mesh)
  X <- matrix(sin(2 * pi * x), n_mesh, 1)

  # Linked: no kappa_mu in possible_params
  op_lk <- hybrid.spde(X = X, beta_X = 1.0,
                      kappa = 5, tau = 1, alpha = 2,
                      loc_mesh = x, d = 1, parameterization = "spde",
                      type = "covariance")
  params_lk <- rSPDE:::extract_possible_parameters(op_lk)
  expect_false("kappa_mu" %in% params_lk)

  # Separate: kappa_mu appears as a parameter
  op_sep <- hybrid.spde(X = X, beta_X = 1.0,
                       kappa = 5, kappa_mu = 10, tau = 1, alpha = 2,
                       loc_mesh = x, d = 1, parameterization = "spde",
                       type = "covariance")
  params_sep <- rSPDE:::extract_possible_parameters(op_sep)
  expect_true("kappa_mu" %in% params_sep)
})


test_that("rspde_lme with fix_kappa_mu recovers parameters", {
  skip_on_cran()
  set.seed(42)

  true_kappa <- 5
  true_tau <- 1.0
  true_kappa_mu <- 10
  true_beta_X <- 1.5
  true_sigma_e <- 0.05

  n_mesh <- 121
  x <- seq(0, 1, length.out = n_mesh)
  X <- matrix(sin(2 * pi * x), n_mesh, 1)
  op_true <- hybrid.spde(X = X, beta_X = true_beta_X,
                        kappa = true_kappa, kappa_mu = true_kappa_mu,
                        tau = true_tau, alpha = 2,
                        loc_mesh = x, d = 1, parameterization = "spde",
                        type = "covariance")

  n_rep <- 10
  n_obs <- 200
  obs.loc <- runif(n_obs, 0, 1)
  A_obs <- make_A(op_true, obs.loc)
  Y <- matrix(NA, n_obs, n_rep)
  for (r in seq_len(n_rep)) {
    u <- as.vector(simulate(op_true, nsim = 1, seed = 100 + r))
    Y[, r] <- as.numeric(A_obs %*% u) + true_sigma_e * rnorm(n_obs)
  }
  dat <- data.frame(
    y = as.vector(Y),
    loc = rep(obs.loc, n_rep),
    repl = rep(seq_len(n_rep), each = n_obs)
  )

  op_start <- hybrid.spde(X = X, beta_X = true_beta_X,
                         kappa = 5, kappa_mu = true_kappa_mu,
                         tau = 1, alpha = 2,
                         loc_mesh = x, d = 1, parameterization = "spde",
                         type = "covariance")

  fit <- rspde_lme(
    y ~ 1, loc = "loc", data = dat,
    model = op_start, repl = "repl", parallel = FALSE,
    model_options = list(fix_alpha = 2,
                         fix_kappa_mu = true_kappa_mu,
                         start_beta_x = true_beta_X)
  )
  re <- fit$coeff$random_effects
  expect_lt(abs(re[["tau"]]    - true_tau)    / true_tau,    0.15)
  expect_lt(abs(re[["kappa"]]  - true_kappa)  / true_kappa,  0.15)
  expect_lt(abs(re[["beta_x1"]] - true_beta_X) / abs(true_beta_X), 0.2)
  expect_lt(abs(as.numeric(fit$coeff$measurement_error) - true_sigma_e) /
              true_sigma_e, 0.1)
})


test_that("hybrid mean is FEM-consistent: mu = beta/kappa^2 for X = 1", {
  # Analytical: L mu = beta X with X(s) = 1 has interior solution
  # mu(s) = beta / kappa^2 (ignoring boundary effects).
  n_mesh <- 81
  x <- seq(0, 1, length.out = n_mesh)
  X <- matrix(1, n_mesh, 1)
  kappa_v <- 10
  beta_v  <- 2

  op <- hybrid.spde(
    X = X, beta_X = beta_v,
    kappa = kappa_v, tau = 1, alpha = 2,
    loc_mesh = x, d = 1, parameterization = "spde",
    type = "covariance"
  )
  mu <- compute_hybrid_mean(op)
  # Check at an interior node, away from the boundary.
  expect_equal(mu[40], beta_v / kappa_v^2, tolerance = 1e-12)
})


test_that("rspde_lme subtracts A %*% mu correctly: hybrid likelihood at truth equals Matern on centered data", {
  # If all hybrid parameters are fixed at the truth, the rspde_lme
  # log-likelihood must coincide with the log-likelihood of the
  # corresponding mean-zero Matern fit on y - A %*% mu, with the SPDE
  # parameters fixed at the same values. Any drift in the mu
  # treatment inside create_likelihood would break this equality.
  set.seed(42)
  true_kappa  <- 8; true_tau <- 1.0
  true_beta_X <- 50; true_sigma_e <- 0.05

  n_mesh <- 121
  x <- seq(0, 1, length.out = n_mesh)
  mesh <- fmesher::fm_mesh_1d(x)
  X <- matrix(sin(2 * pi * x), n_mesh, 1)
  op_true <- hybrid.spde(
    X = X, beta_X = true_beta_X,
    kappa = true_kappa, tau = true_tau, alpha = 2.0,
    loc_mesh = x, d = 1, parameterization = "spde",
    type = "covariance"
  )

  n_obs <- 300
  obs.loc <- runif(n_obs, 0, 1)
  A_obs <- make_A(op_true, obs.loc)
  u <- as.vector(simulate(op_true, nsim = 1, seed = 99))
  y <- as.numeric(A_obs %*% u) + true_sigma_e * rnorm(n_obs)
  mu_at_obs <- as.numeric(A_obs %*% compute_hybrid_mean(op_true))

  # Hybrid likelihood at truth (with everything fixed).
  dat <- data.frame(y = y, loc = obs.loc)
  fit_h <- rspde_lme(
    y ~ 0, loc = "loc", data = dat,
    model = op_true, parallel = FALSE,
    model_options = list(
      fix_alpha = 2, fix_tau = true_tau, fix_kappa = true_kappa,
      fix_beta_x = true_beta_X, fix_sigma_e = true_sigma_e
    )
  )

  # Pure Matern likelihood on centered data (y - A %*% mu).
  op_m <- spde.matern.operators(
    kappa = true_kappa, tau = true_tau, alpha = 2,
    loc_mesh = x, d = 1, parameterization = "spde",
    type = "covariance"
  )
  dat2 <- data.frame(y = y - mu_at_obs, loc = obs.loc)
  fit_b <- rspde_lme(
    y ~ 0, loc = "loc", data = dat2,
    model = op_m, parallel = FALSE,
    model_options = list(
      fix_alpha = 2, fix_tau = true_tau, fix_kappa = true_kappa,
      fix_sigma_e = true_sigma_e
    )
  )

  expect_equal(fit_h$loglik, fit_b$loglik, tolerance = 1e-8)
})


test_that("Operator-based and covariance-based hybrids share the same mean", {
  s <- make_simple_setup(n = 81)
  op_cov <- hybrid.spde(
    X = s$X, beta_X = s$beta_X,
    kappa = s$kappa, tau = s$tau, alpha = s$alpha,
    loc_mesh = s$x, d = 1, parameterization = "spde",
    type = "covariance"
  )
  op_op <- hybrid.spde(
    X = s$X, beta_X = s$beta_X,
    kappa = s$kappa, tau = s$tau, alpha = s$alpha,
    loc_mesh = s$x, d = 1, parameterization = "spde",
    type = "operator"
  )
  expect_equal(as.numeric(op_cov$mu), as.numeric(op_op$mu),
               tolerance = 1e-10)
})


test_that("simulate.hybrid_spde with beta_X = 0 matches the base simulate", {
  s <- make_simple_setup(n = 61)

  op_h <- hybrid.spde(
    X = s$X, beta_X = 0,
    kappa = s$kappa, tau = s$tau, alpha = s$alpha,
    loc_mesh = s$x, d = 1, parameterization = "spde",
    type = "covariance"
  )
  set.seed(42)
  Y_h <- simulate(op_h, nsim = 1)

  op_base <- spde.matern.operators(
    kappa = s$kappa, tau = s$tau, alpha = s$alpha,
    loc_mesh = s$x, d = 1, parameterization = "spde",
    type = "covariance"
  )
  set.seed(42)
  Y_base <- simulate(op_base, nsim = 1)

  expect_equal(as.numeric(Y_h), as.numeric(Y_base), tolerance = 1e-10)
})


test_that("simulate.hybrid_spde adds the mean to simulated fields", {
  s <- make_simple_setup(n = 81, beta_X = 2)

  op_h <- hybrid.spde(
    X = s$X, beta_X = s$beta_X,
    kappa = s$kappa, tau = s$tau, alpha = s$alpha,
    loc_mesh = s$x, d = 1, parameterization = "spde",
    type = "covariance"
  )
  op_zero <- hybrid.spde(
    X = s$X, beta_X = 0,
    kappa = s$kappa, tau = s$tau, alpha = s$alpha,
    loc_mesh = s$x, d = 1, parameterization = "spde",
    type = "covariance"
  )

  set.seed(1)
  Y_h <- simulate(op_h, nsim = 1)
  set.seed(1)
  Y_zero <- simulate(op_zero, nsim = 1)

  diff <- as.numeric(Y_h) - as.numeric(Y_zero)
  expect_equal(diff, as.numeric(op_h$mu), tolerance = 1e-10)
})


test_that("update.hybrid_spde refreshes the mean and rebuilds op_mean", {
  s <- make_simple_setup(n = 81, beta_X = 1.5)
  op_h <- hybrid.spde(
    X = s$X, beta_X = s$beta_X,
    kappa = s$kappa, tau = s$tau, alpha = s$alpha,
    loc_mesh = s$x, d = 1, parameterization = "spde",
    type = "covariance"
  )

  op_h2 <- update(op_h, beta_X = 3, check_stationarity = FALSE)
  expect_equal(as.numeric(op_h2$beta_X), 3)
  # The mean scales linearly with beta_X.
  expect_equal(as.numeric(op_h2$mu),
               2 * as.numeric(op_h$mu),
               tolerance = 1e-10)

  # tau does NOT enter the mean in this parameterisation.
  op_h3 <- update(op_h, tau = 2 * s$tau, check_stationarity = FALSE)
  expect_equal(as.numeric(op_h3$mu),
               as.numeric(op_h$mu),
               tolerance = 1e-10)
})


test_that("predict.hybrid_spde recovers the data mean for many samples", {
  set.seed(7)
  s <- make_simple_setup(n = 81, beta_X = 1.2, alpha = 1.3,
                         kappa = 8, tau = 0.5)

  op_h <- hybrid.spde(
    X = s$X, beta_X = s$beta_X,
    kappa = s$kappa, tau = s$tau, alpha = s$alpha,
    loc_mesh = s$x, d = 1, parameterization = "spde",
    type = "covariance"
  )

  obs.loc <- runif(40, min = 0, max = 1)
  sigma.e <- 0.1
  A_obs <- make_A(op_h, obs.loc)
  Aprd <- make_A(op_h, s$x)

  u <- as.vector(simulate(op_h, nsim = 1, seed = 11))
  y <- as.numeric(A_obs %*% u) + sigma.e * rnorm(length(obs.loc))

  pred <- predict(op_h, A = A_obs, Aprd = Aprd, Y = y,
                  sigma.e = sigma.e, compute.variances = TRUE)

  expect_equal(length(as.numeric(pred$mean)), length(s$x))
  expect_true(all(is.finite(as.numeric(pred$mean))))
  expect_true(all(pred$variance >= 0))

  # Sanity check: the prediction should be closer to the truth than the
  # prior mean alone (which already incorporates the deterministic
  # component beta_X L^{-alpha/2} X / tau).
  rmse_pred <- sqrt(mean((as.numeric(pred$mean) - u)^2))
  rmse_prior <- sqrt(mean((as.numeric(op_h$mu) - u)^2))
  expect_lt(rmse_pred, rmse_prior)
})


test_that("rspde_lme fits a hybrid.spde model with the mean built in", {
  skip_on_cran()
  set.seed(123)

  n_mesh <- 81
  s <- make_simple_setup(n = n_mesh, kappa = 8, tau = 0.5,
                         alpha = 1.3, beta_X = 1.5)

  op_h <- hybrid.spde(
    X = s$X, beta_X = s$beta_X,
    kappa = s$kappa, tau = s$tau, alpha = s$alpha,
    loc_mesh = s$x, d = 1, parameterization = "spde",
    type = "covariance"
  )

  # Simulate data with several replicates for stability.
  n_rep <- 5
  n_obs <- 80
  obs.loc <- runif(n_obs, min = 0, max = 1)
  A_obs <- make_A(op_h, obs.loc)
  sigma.e <- 0.1
  Y <- matrix(NA, nrow = n_obs, ncol = n_rep)
  for (r in seq_len(n_rep)) {
    u <- as.vector(simulate(op_h, nsim = 1, seed = 100 + r))
    Y[, r] <- as.numeric(A_obs %*% u) + sigma.e * rnorm(n_obs)
  }
  y_vec <- as.vector(Y)
  repl <- rep(seq_len(n_rep), each = n_obs)
  loc_vec <- rep(obs.loc, times = n_rep)

  dat <- data.frame(y = y_vec, loc = loc_vec, repl = repl)
  fit <- rspde_lme(
    y ~ 1, loc = "loc", data = dat,
    model = op_h, repl = "repl",
    parallel = FALSE
  )

  expect_s3_class(fit, "rspde_lme")
  expect_true(is.finite(fit$loglik))

  # The estimated SPDE parameters should be in a sensible range.
  re <- fit$coeff$random_effects
  expect_true(all(re > 0))
  expect_true(is.finite(fit$coeff$measurement_error))
})


test_that("rspde_lme recovers all parameters of a simulated hybrid.spde model", {
  skip_on_cran()
  set.seed(2026)

  # True parameters. beta_X here is interpreted in the FEM-consistent
  # convention where mu_h = beta * L_d^{-1} C X (with C the mass-lumped
  # FEM mass matrix), so the deterministic mean has magnitude
  # ~ beta / kappa^2 for X of order 1. We pick beta large enough that
  # the deterministic signal is comparable to the measurement noise
  # (i.e. recoverable from data).
  true_kappa   <- 8
  true_tau     <- 0.5
  true_alpha   <- 1.3
  true_beta_X  <- c(150, -70)
  true_sigma_e <- 0.1

  # 1d mesh and two covariate fields at the mesh nodes
  n_mesh <- 121
  x      <- seq(from = 0, to = 1, length.out = n_mesh)
  X_mat  <- cbind(sin(2 * pi * x), cos(2 * pi * x))

  op_true <- hybrid.spde(
    X = X_mat, beta_X = true_beta_X,
    kappa = true_kappa, tau = true_tau, alpha = true_alpha,
    loc_mesh = x, d = 1, parameterization = "spde",
    type = "covariance", m = 2
  )

  # Simulate several replicates of the field with measurement noise
  n_rep <- 20
  n_obs <- 150
  obs.loc <- runif(n_obs, min = 0, max = 1)
  A_obs   <- make_A(op_true, obs.loc)
  Y <- matrix(NA, nrow = n_obs, ncol = n_rep)
  for (r in seq_len(n_rep)) {
    u_r    <- as.vector(simulate(op_true, nsim = 1, seed = 1000 + r))
    Y[, r] <- as.numeric(A_obs %*% u_r) + true_sigma_e * rnorm(n_obs)
  }
  dat <- data.frame(
    y    = as.vector(Y),
    loc  = rep(obs.loc, times = n_rep),
    repl = rep(seq_len(n_rep), each = n_obs)
  )

  # Fit from a deliberately incorrect SPDE starting point. With the
  # FEM-consistent mu formula (mu_h = beta * L_d^{-1} C X) the
  # gradient of the likelihood with respect to beta at beta = 0 is
  # mild, so we seed the optimiser at a sensible scale via
  # start_beta_x. In practice users would derive such a starting
  # value from an OLS regression of y on the projected covariates.
  op_start <- hybrid.spde(
    X = X_mat, beta_X = true_beta_X,  # value used only by simulate/predict
    kappa = 5, tau = 1, alpha = 1.5,
    loc_mesh = x, d = 1, parameterization = "spde",
    type = "covariance", m = 2
  )
  fit <- rspde_lme(
    y ~ 1, loc = "loc", data = dat,
    model = op_start, repl = "repl",
    parallel = FALSE,
    model_options = list(start_beta_x = true_beta_X)
  )

  re <- fit$coeff$random_effects
  est <- c(tau = re[["tau"]], kappa = re[["kappa"]], alpha = re[["alpha"]],
           beta_x1 = re[["beta_x1"]], beta_x2 = re[["beta_x2"]],
           sigma_e = as.numeric(fit$coeff$measurement_error))
  truth <- c(tau = true_tau, kappa = true_kappa, alpha = true_alpha,
             beta_x1 = true_beta_X[1], beta_x2 = true_beta_X[2],
             sigma_e = true_sigma_e)
  rel_err <- abs(est - truth) / abs(truth)

  # All six estimated parameters should be recovered to within ~10%.
  expect_lt(rel_err[["tau"]],     0.1)
  expect_lt(rel_err[["kappa"]],   0.1)
  expect_lt(rel_err[["alpha"]],   0.1)
  expect_lt(rel_err[["beta_x1"]], 0.1)
  expect_lt(rel_err[["beta_x2"]], 0.1)
  expect_lt(rel_err[["sigma_e"]], 0.1)
})


test_that("fix_beta_x in model_options holds beta_X fixed during estimation", {
  skip_on_cran()
  set.seed(2027)

  true_kappa   <- 8
  true_tau     <- 0.5
  true_alpha   <- 1.3
  true_beta_X  <- c(1.5, -0.7)
  true_sigma_e <- 0.1

  n_mesh <- 121
  x      <- seq(0, 1, length.out = n_mesh)
  X_mat  <- cbind(sin(2 * pi * x), cos(2 * pi * x))

  op_true <- hybrid.spde(
    X = X_mat, beta_X = true_beta_X,
    kappa = true_kappa, tau = true_tau, alpha = true_alpha,
    loc_mesh = x, d = 1, parameterization = "spde",
    type = "covariance", m = 2
  )
  n_rep <- 10
  n_obs <- 120
  obs.loc <- runif(n_obs, 0, 1)
  A_obs <- make_A(op_true, obs.loc)
  Y <- matrix(NA, n_obs, n_rep)
  for (r in seq_len(n_rep)) {
    u_r <- as.vector(simulate(op_true, nsim = 1, seed = 2000 + r))
    Y[, r] <- as.numeric(A_obs %*% u_r) + true_sigma_e * rnorm(n_obs)
  }
  dat <- data.frame(y = as.vector(Y),
                    loc = rep(obs.loc, n_rep),
                    repl = rep(seq_len(n_rep), each = n_obs))

  # Vector-valued fix_beta_x freezes all beta_X coefficients at the
  # supplied values.
  fit_fixed <- rspde_lme(
    y ~ 1, loc = "loc", data = dat,
    model = op_true, repl = "repl",
    parallel = FALSE,
    model_options = list(fix_beta_x = true_beta_X)
  )
  re <- fit_fixed$coeff$random_effects
  # The estimated beta_x* should be exactly the values we fixed.
  expect_equal(as.numeric(re[["beta_x1 (fixed)"]]), true_beta_X[1],
               tolerance = 1e-12)
  expect_equal(as.numeric(re[["beta_x2 (fixed)"]]), true_beta_X[2],
               tolerance = 1e-12)

  # Per-component fixing also works.
  fit_partial <- rspde_lme(
    y ~ 1, loc = "loc", data = dat,
    model = op_true, repl = "repl",
    parallel = FALSE,
    model_options = list(fix_beta_x1 = true_beta_X[1])
  )
  re2 <- fit_partial$coeff$random_effects
  expect_equal(as.numeric(re2[["beta_x1 (fixed)"]]), true_beta_X[1],
               tolerance = 1e-12)
  # beta_x2 was free, so it carries no "(fixed)" tag
  expect_true("beta_x2" %in% names(re2))
})


test_that("rspde_lme on hybrid model with beta_X fixed at 0 matches the underlying SPDE", {
  skip_on_cran()
  set.seed(321)

  n_mesh <- 81
  s <- make_simple_setup(n = n_mesh, alpha = 1.3, beta_X = 0)

  op_h <- hybrid.spde(
    X = s$X, beta_X = 0,
    kappa = s$kappa, tau = s$tau, alpha = s$alpha,
    loc_mesh = s$x, d = 1, parameterization = "spde",
    type = "covariance"
  )
  op_base <- spde.matern.operators(
    kappa = s$kappa, tau = s$tau, alpha = s$alpha,
    loc_mesh = s$x, d = 1, parameterization = "spde",
    type = "covariance"
  )

  n_obs <- 60
  obs.loc <- runif(n_obs, min = 0, max = 1)
  A_obs <- make_A(op_base, obs.loc)
  sigma.e <- 0.1
  u <- as.vector(simulate(op_base, nsim = 1, seed = 99))
  y <- as.numeric(A_obs %*% u) + sigma.e * rnorm(n_obs)

  dat <- data.frame(y = y, loc = obs.loc)
  fit_h <- rspde_lme(y ~ 1, loc = "loc", data = dat,
                     model = op_h, parallel = FALSE,
                     model_options = list(fix_beta_x = 0))
  fit_b <- rspde_lme(y ~ 1, loc = "loc", data = dat,
                     model = op_base, parallel = FALSE)

  expect_equal(fit_h$loglik, fit_b$loglik, tolerance = 1e-6)
})
