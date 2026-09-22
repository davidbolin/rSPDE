context("rspde.hybrid.matern (true INLA cgeneric, alpha = 2)")

# ---------------------------------------------------------------------------
# Helper to set up the simulation in d = 1 with alpha = 2 (so nu = 1.5).
# ---------------------------------------------------------------------------
make_inla_hybrid_data <- function(seed       = 2026,
                                  true_kappa = 5,
                                  true_tau   = 1.0,
                                  true_beta_X = c(1.5, -0.7),
                                  true_sigma_e = 0.05,
                                  n_mesh = 101,
                                  n_obs  = 300) {
  set.seed(seed)
  x <- seq(0, 1, length.out = n_mesh)
  mesh <- fmesher::fm_mesh_1d(x)
  X <- cbind(sin(2 * pi * x), cos(2 * pi * x))

  op_true <- hybrid.spde(
    X = X, beta_X = true_beta_X,
    kappa = true_kappa, tau = true_tau, alpha = 2.0,
    loc_mesh = x, d = 1, parameterization = "spde",
    type = "covariance"
  )

  obs.loc <- runif(n_obs, 0, 1)
  A_obs <- make_A(op_true, obs.loc)
  u <- as.vector(simulate(op_true, nsim = 1, seed = seed + 100))
  y <- as.numeric(A_obs %*% u) + true_sigma_e * rnorm(n_obs)

  list(mesh = mesh, X = X, x = x, obs.loc = obs.loc, y = y,
       op_true = op_true,
       true_kappa = true_kappa, true_tau = true_tau,
       true_beta_X = true_beta_X, true_sigma_e = true_sigma_e,
       n_mesh = n_mesh)
}


test_that("rspde.hybrid.matern builds with expected structure", {
  local_rspde_safe_inla(
    required_symbol = "inla_cgeneric_rspde_hybrid_alpha2_model"
  )
  s <- make_inla_hybrid_data(n_obs = 30, n_mesh = 31)
  hyb <- rspde.hybrid.matern(mesh = s$mesh, X = s$X)
  expect_s3_class(hyb, "inla_rspde_hybrid_alpha2")
  expect_s3_class(hyb, "inla_rspde")
  expect_equal(hyb$p, 2L)
  expect_equal(hyb$alpha, 2)
  expect_equal(hyb$n.spde, s$n_mesh)
  expect_equal(hyb$rspde.order, 0L)
  expect_true(hyb$integer.nu)
  expect_equal(hyb$d, 1L)
  expect_equal(hyb$nu, 1.5)
})


test_that("Q from the hybrid cgeneric matches the FEM formula for alpha = 2", {
  skip_on_cran()
  local_rspde_safe_inla()

  s <- make_inla_hybrid_data(n_obs = 10, n_mesh = 21)
  # rspde.matern with fixed nu = 1.5 (alpha = 2) — known correct Q
  reg <- rspde.matern(
    mesh = s$mesh, nu = 1.5, parameterization = "spde",
    start.ltau = log(s$true_tau), start.lkappa = log(s$true_kappa)
  )
  Q_reg <- INLA::inla.cgeneric.q(reg)$Q

  # Manual: Q = tau^2 * (G + kappa^2 C) * Ci * (G + kappa^2 C)
  fem <- fmesher::fm_fem(s$mesh, order = 2)
  C_lumped <- Matrix::Diagonal(s$n_mesh, x = Matrix::rowSums(fem$c0))
  Cinv     <- Matrix::Diagonal(s$n_mesh, x = 1 / Matrix::rowSums(fem$c0))
  L_manual <- fem$g1 + s$true_kappa^2 * C_lumped
  Q_manual <- s$true_tau^2 * L_manual %*% Cinv %*% L_manual

  # rspde.matern's Q should match the manual L*Ci*L formula.
  expect_lt(max(abs(as.matrix(Q_reg) - as.matrix(Q_manual))), 1e-8)
})


test_that("Hybrid cgeneric mu formula matches compute_hybrid_mean", {
  # Verify that the C++ Eigen helper computes the same mean as the
  # operator-based R routine compute_hybrid_mean. We compare against a
  # manual L^{-1} (X * beta) computation that uses the same Q-formula
  # the C++ aux uses.
  set.seed(1)
  n_mesh <- 21
  x <- seq(0, 1, length.out = n_mesh)
  mesh <- fmesher::fm_mesh_1d(x)
  X <- cbind(sin(2 * pi * x), cos(2 * pi * x))
  kappa_v <- 5
  beta_v  <- c(1.5, -0.7)

  op <- hybrid.spde(
    X = X, beta_X = beta_v,
    kappa = kappa_v, tau = 1.0, alpha = 2.0,
    loc_mesh = x, d = 1, parameterization = "spde",
    type = "covariance"
  )
  mu_R <- compute_hybrid_mean(op)

  # Manual reference for the FEM variational form:
  #   mu_h = L_d^{-1} (C * X * beta), with C the mass-lumped diagonal.
  fem <- fmesher::fm_fem(mesh, order = 2)
  C_diag <- Matrix::rowSums(fem$c0)
  C_lumped <- Matrix::Diagonal(n_mesh, x = C_diag)
  L_manual <- fem$g1 + kappa_v^2 * C_lumped
  rhs <- as.vector(C_diag * as.vector(X %*% beta_v))
  mu_manual <- as.vector(Matrix::solve(L_manual, rhs))

  expect_equal(as.numeric(mu_R), as.numeric(mu_manual), tolerance = 1e-12)
})


test_that("INLA recovers beta_X from simulated data (alpha = 2)", {
  skip_on_cran()
  local_rspde_safe_inla(
    required_symbol = "inla_cgeneric_rspde_hybrid_alpha2_model"
  )
  suppressMessages(library(INLA))

  # Simulate from a hybrid SPDE with non-zero beta_X. With diffuse
  # priors on (tau, kappa) the joint posterior of (tau, kappa, beta_X)
  # has a degenerate mode in which the random field absorbs the
  # deterministic signal; that is an intrinsic identifiability
  # property of the model class. Using moderately informative priors
  # on (tau, kappa) centred at sensible values breaks the degeneracy
  # and INLA recovers beta_X (and the other parameters) accurately.
  set.seed(2026)

  # The deterministic mean is mu_h = beta * L_d^{-1} C X, so the
  # mean magnitude scales as |beta| / kappa^2. We pick beta large
  # enough that the signal is comparable to the measurement noise.
  true_kappa   <- 5
  true_tau     <- 1.0
  true_beta_X  <- c(40, -20)
  true_sigma_e <- 0.05

  n_mesh <- 101
  x <- seq(0, 1, length.out = n_mesh)
  mesh <- fmesher::fm_mesh_1d(x)
  X <- cbind(sin(2 * pi * x), cos(2 * pi * x))

  op_true <- hybrid.spde(
    X = X, beta_X = true_beta_X,
    kappa = true_kappa, tau = true_tau, alpha = 2,
    loc_mesh = x, d = 1, parameterization = "spde",
    type = "covariance"
  )
  # n_obs = 1000, not 300. At 300 observations the degenerate mode described
  # above is not merely a competing mode, it is the better one: INLA lands on
  # tau = 0.10, beta_X = (29.4, -8.7) with a marginal likelihood of 458.2
  # against 453.7 for the mode at the truth, so this test asserted something
  # false about that data set and failed whenever the optimiser found the
  # better fit. That is a question of information, not of numerics -- pinning
  # num.threads makes it deterministic but deterministically wrong. With 1000
  # observations the mode at the truth wins outright and the recovery is
  # comfortable: over six runs the relative errors were tau 0.06-0.22 (the
  # tolerance is 0.30) and at most 0.013 for everything else.
  n_obs <- 1000
  obs.loc <- runif(n_obs, 0, 1)
  A_obs <- make_A(op_true, obs.loc)
  u <- as.vector(simulate(op_true, nsim = 1, seed = 200))
  y <- as.numeric(A_obs %*% u) + true_sigma_e * rnorm(n_obs)

  hyb <- rspde.hybrid.matern(
    mesh = mesh, X = X,
    start.ltau   = log(true_tau),
    start.lkappa = log(true_kappa),
    start.beta_x = true_beta_X,
    # Informative priors on (tau, kappa); diffuse prior on beta_X.
    prior.tau    = list(mean = log(true_tau),   prec = 10),
    prior.kappa  = list(mean = log(true_kappa), prec = 10),
    prior.beta_x = list(mean = c(0, 0), prec = c(0.001, 0.001))
  )
  A <- fmesher::fm_basis(mesh, obs.loc)
  stk <- inla.stack(data = list(y = y), A = list(A),
                    effects = list(idx = 1:n_mesh))
  fit <- inla(y ~ -1 + f(idx, model = hyb),
              data = inla.stack.data(stk),
              family = "gaussian",
              control.predictor = list(A = inla.stack.A(stk)),
              control.inla = list(int.strategy = "eb"))

  hp <- fit$summary.hyperpar
  prec_idx  <- grep("Precision", rownames(hp))
  theta_idx <- grep("Theta", rownames(hp))
  est <- c(
    tau     = as.numeric(exp(hp[theta_idx[1], "mean"])),
    kappa   = as.numeric(exp(hp[theta_idx[2], "mean"])),
    beta_x1 = as.numeric(hp[theta_idx[3], "mean"]),
    beta_x2 = as.numeric(hp[theta_idx[4], "mean"]),
    sigma_e = as.numeric(1 / sqrt(hp[prec_idx, "mean"]))
  )
  truth <- c(
    tau = true_tau, kappa = true_kappa,
    beta_x1 = true_beta_X[1], beta_x2 = true_beta_X[2],
    sigma_e = true_sigma_e
  )
  rel_err <- abs(est - truth) / abs(truth)

  expect_lt(rel_err[["tau"]],     0.30)   # tau is the loosest
  expect_lt(rel_err[["kappa"]],   0.10)
  expect_lt(rel_err[["beta_x1"]], 0.15)
  expect_lt(rel_err[["beta_x2"]], 0.15)
  expect_lt(rel_err[["sigma_e"]], 0.10)
})


test_that("separate_kappa_mu adds an extra hyperparameter and recovers kappa_mu", {
  skip_on_cran()
  local_rspde_safe_inla(
    required_symbol = "inla_cgeneric_rspde_hybrid_alpha2_model"
  )
  suppressMessages(library(INLA))

  # Simulate from a model with separate kappa_mu = 10 (while kappa = 5),
  # fit with the cgeneric set to separate mode and tight priors at
  # truth on (tau, kappa, kappa_mu); check that kappa_mu is recovered
  # near 10 and the model exposes 6 hyperparameters (Gaussian prec +
  # ltau + lkappa + lkappa_mu + 2 beta_x).
  set.seed(2026)

  true_kappa    <- 5
  true_tau      <- 1.0
  true_kappa_mu <- 10
  true_beta_X   <- c(1.5, -0.7)
  true_sigma_e  <- 0.05

  n_mesh <- 121
  x <- seq(0, 1, length.out = n_mesh)
  mesh <- fmesher::fm_mesh_1d(x)
  X <- cbind(sin(2 * pi * x), cos(2 * pi * x))

  op_true <- hybrid.spde(X = X, beta_X = true_beta_X,
                        kappa = true_kappa, kappa_mu = true_kappa_mu,
                        tau = true_tau, alpha = 2,
                        loc_mesh = x, d = 1, parameterization = "spde",
                        type = "covariance")
  n_obs <- 300
  obs.loc <- runif(n_obs, 0, 1)
  A_obs <- make_A(op_true, obs.loc)
  u <- as.vector(simulate(op_true, nsim = 1, seed = 200))
  y <- as.numeric(A_obs %*% u) + true_sigma_e * rnorm(n_obs)
  A <- fmesher::fm_basis(mesh, obs.loc)
  stk <- inla.stack(data = list(y = y), A = list(A),
                    effects = list(idx = 1:n_mesh))

  hyb <- rspde.hybrid.matern(
    mesh = mesh, X = X,
    separate_kappa_mu = TRUE,
    start.ltau = log(true_tau),
    start.lkappa = log(true_kappa),
    start.beta_x = true_beta_X,
    start.lkappa_mu = log(true_kappa_mu),
    prior.tau      = list(mean = log(true_tau),      prec = 100),
    prior.kappa    = list(mean = log(true_kappa),    prec = 100),
    prior.kappa_mu = list(mean = log(true_kappa_mu), prec = 100),
    prior.beta_x   = list(mean = c(0, 0), prec = c(1, 1))
  )
  expect_true(isTRUE(hyb$separate_kappa_mu))

  fit <- inla(y ~ -1 + f(idx, model = hyb),
              data = inla.stack.data(stk),
              family = "gaussian",
              control.predictor = list(A = inla.stack.A(stk)),
              control.inla = list(int.strategy = "eb"))

  # 6 hyperparameters: Gaussian prec + ltau + lkappa + lkappa_mu + 2 beta_x
  expect_equal(nrow(fit$summary.hyperpar), 6L)
  hp <- fit$summary.hyperpar
  theta_idx <- grep("Theta", rownames(hp))
  expect_equal(length(theta_idx), 5L)

  # kappa_mu should be recovered within 15% of truth.
  est_kappa_mu <- as.numeric(exp(hp[theta_idx[3], "mean"]))
  expect_lt(abs(est_kappa_mu - true_kappa_mu) / true_kappa_mu, 0.15)
})


test_that("Hybrid with beta_X = 0 fixed gives same fit as rspde.matern (nu=1.5)", {
  skip_on_cran()
  local_rspde_safe_inla(
    required_symbol = "inla_cgeneric_rspde_hybrid_alpha2_model"
  )
  suppressMessages(library(INLA))

  set.seed(11)
  n_mesh <- 81
  x <- seq(0, 1, length.out = n_mesh)
  mesh <- fmesher::fm_mesh_1d(x)

  # Generate from a pure (no mean) stationary Matern alpha=2 model.
  kappa_true <- 5
  tau_true   <- 1.0
  sigma_e    <- 0.1
  op <- spde.matern.operators(
    kappa = kappa_true, tau = tau_true, alpha = 2,
    loc_mesh = x, d = 1, parameterization = "spde",
    type = "covariance"
  )
  n_obs <- 200
  obs.loc <- runif(n_obs, 0, 1)
  A_obs <- make_A(op, obs.loc)
  u <- as.vector(simulate(op, nsim = 1, seed = 99))
  y <- as.numeric(A_obs %*% u) + sigma_e * rnorm(n_obs)

  A <- fmesher::fm_basis(mesh, obs.loc)
  stk <- inla.stack(data = list(y = y), A = list(A),
                    effects = list(idx = 1:n_mesh))

  # Common (ltau, lkappa) prior and starting values.
  prior.tau   <- list(mean = 0, prec = 0.1)
  prior.kappa <- list(mean = log(sqrt(8 * 0.5) / 0.2), prec = 0.1)
  ltau0       <- log(tau_true)
  lkappa0     <- log(kappa_true)

  X <- matrix(0, n_mesh, 1)  # dummy; beta_X frozen at 0
  hyb <- rspde.hybrid.matern(
    mesh = mesh, X = X,
    start.ltau   = ltau0,
    start.lkappa = lkappa0,
    start.beta_x = 0,
    prior.tau    = prior.tau,
    prior.kappa  = prior.kappa,
    prior.beta_x = list(mean = 0, prec = 1e12)
  )
  reg <- rspde.matern(
    mesh = mesh, nu = 1.5, parameterization = "spde",
    start.ltau = ltau0, start.lkappa = lkappa0,
    theta.prior.mean = c(prior.tau$mean, prior.kappa$mean),
    theta.prior.prec = diag(c(prior.tau$prec, prior.kappa$prec))
  )

  fit_h <- inla(y ~ -1 + f(idx, model = hyb),
                data = inla.stack.data(stk),
                family = "gaussian",
                control.predictor = list(A = inla.stack.A(stk)),
                control.inla = list(int.strategy = "eb"))
  fit_r <- inla(y ~ -1 + f(idx, model = reg),
                data = inla.stack.data(stk),
                family = "gaussian",
                control.predictor = list(A = inla.stack.A(stk)),
                control.inla = list(int.strategy = "eb"))

  hp_h <- fit_h$summary.hyperpar
  hp_r <- fit_r$summary.hyperpar
  prec_idx_h  <- grep("Precision", rownames(hp_h))
  theta_idx_h <- grep("Theta", rownames(hp_h))
  prec_idx_r  <- grep("Precision", rownames(hp_r))
  theta_idx_r <- grep("Theta", rownames(hp_r))

  # The claim is that these are the same model, and the way to test it is to
  # evaluate both at the same theta. Comparing two separate fits does not test
  # it: with alpha = 2 in one dimension log(tau) and log(kappa) are only weakly
  # identified (they trade off along a flat ridge), and the hybrid carries an
  # extra frozen hyperparameter that makes INLA stop at a different point on
  # that ridge from one run to the next. The individual hyperparameters then
  # differed by 0.05 to 0.29 at random, which is what used to fail here, while
  # the fitted values never differed by more than 3 per cent. Tightening the
  # prior does not help and neither does more data, because the ridge is a
  # property of the parameterisation.
  #
  # Held at a common theta the two are identical to roundoff.
  th0 <- c(log(1 / sigma_e^2), ltau0, lkappa0)
  fixed_fit <- function(mod, theta) {
    inla(y ~ -1 + f(idx, model = mod),
         data = inla.stack.data(stk), family = "gaussian",
         control.predictor = list(A = inla.stack.A(stk)),
         control.inla = list(int.strategy = "eb"),
         ## restart = FALSE is implied by fixed = TRUE; stating it keeps INLA
         ## from warning about having set it.
         control.mode = list(theta = theta, fixed = TRUE, restart = FALSE))
  }
  fix_h <- fixed_fit(hyb, c(th0, 0))
  fix_r <- fixed_fit(reg, th0)
  idx_o <- seq_len(n_obs)
  expect_equal(fix_h$summary.fitted.values[idx_o, "mean"],
               fix_r$summary.fitted.values[idx_o, "mean"], tolerance = 1e-8)
  expect_equal(fix_h$summary.random$idx$mean,
               fix_r$summary.random$idx$mean, tolerance = 1e-8)
  expect_lt(abs(fix_h$mlik[1, 1] - fix_r$mlik[1, 1]), 1e-5)

  # Fitted freely, only what does not depend on where the optimiser stopped.
  fv_h <- fit_h$summary.fitted.values[idx_o, "mean"]
  fv_r <- fit_r$summary.fitted.values[idx_o, "mean"]
  expect_lt(max(abs(fv_h - fv_r)) / max(abs(fv_r)), 0.06)
  sig_h <- 1 / sqrt(hp_h[prec_idx_h, "mean"])
  sig_r <- 1 / sqrt(hp_r[prec_idx_r, "mean"])
  expect_lt(abs(sig_h - sig_r) / sig_r, 0.10)
})


test_that("Hybrid cgeneric fits a full INLA model and recovers sigma_e", {
  skip_on_cran()
  local_rspde_safe_inla(
    required_symbol = "inla_cgeneric_rspde_hybrid_alpha2_model"
  )
  suppressMessages(library(INLA))

  s <- make_inla_hybrid_data()
  hyb <- rspde.hybrid.matern(
    mesh = s$mesh, X = s$X,
    start.ltau   = log(s$true_tau),
    start.lkappa = log(s$true_kappa),
    start.beta_x = c(0, 0)
  )

  A <- fmesher::fm_basis(s$mesh, s$obs.loc)
  stk <- inla.stack(data = list(y = s$y), A = list(A),
                    effects = list(idx = 1:s$n_mesh))
  fit <- inla(y ~ -1 + f(idx, model = hyb),
              data = inla.stack.data(stk),
              family = "gaussian",
              control.predictor = list(A = inla.stack.A(stk)),
              control.inla = list(int.strategy = "eb"))

  expect_true(is.finite(fit$mlik[1, 1]))
  expect_true(nrow(fit$summary.hyperpar) == 5L)   # prec, ltau, lkappa, b1, b2
  # sigma_e is well-identified by the noise structure even when the
  # latent decomposition is multi-modal -- it should be close to truth.
  prec_idx <- grep("Precision", rownames(fit$summary.hyperpar))
  sigma_e_est <- 1 / sqrt(fit$summary.hyperpar[prec_idx, "mean"])
  expect_lt(abs(sigma_e_est - s$true_sigma_e) / s$true_sigma_e, 0.5)
})


test_that("Hybrid model works with inlabru via the inla_rspde mapper", {
  skip_on_cran()
  local_rspde_safe_inla(
    required_symbol = "inla_cgeneric_rspde_hybrid_alpha2_model"
  )
  skip_if_not_installed("inlabru")
  suppressMessages({library(INLA); library(inlabru)})

  s <- make_inla_hybrid_data(n_obs = 60, n_mesh = 41)
  hyb <- rspde.hybrid.matern(
    mesh = s$mesh, X = s$X,
    start.ltau   = log(s$true_tau),
    start.lkappa = log(s$true_kappa),
    start.beta_x = c(0, 0)
  )
  df <- data.frame(y = s$y, loc = s$obs.loc)

  cmp <- y ~ -1 + field(loc, model = hyb)
  fit <- bru(cmp, data = df,
             family = "gaussian",
             options = list(num.threads = "1:1",
                            control.inla = list(int.strategy = "eb")))

  expect_true(is.finite(fit$mlik[1, 1]))
  expect_true(nrow(fit$summary.hyperpar) == 5L)
})
