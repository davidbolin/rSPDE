context("posterior_crossvalidation for rspde_lme")

# Helper: fit a small non-graph rspde_lme on a 1D mesh using a CBrSPDEobj
# (covariance-based) operator, which is the standard path the cross-
# validation code is exercised against.
.fit_rspde_lme_1d <- function(seed = 1L, n = 40L) {
  set.seed(seed)
  x <- seq(0, 10, length.out = n)
  model <- rSPDE::matern.operators(range = 2, sigma = 1, nu = 0.8,
                                   loc_mesh = seq(0, 10, length.out = 150),
                                   d = 1, type = "covariance",
                                   parameterization = "matern")
  u    <- as.vector(simulate(model, nsim = 1, seed = seed))
  Aobs <- rSPDE::make_A(model, x)
  y    <- as.vector(Aobs %*% u) + 0.3 * stats::rnorm(n)
  dat  <- data.frame(y = y, x = x)
  fit  <- suppressWarnings(
    rspde_lme(y ~ 1, loc = "x", data = dat, model = model,
              optim_controls = list(maxit = 20),
              model_options = list(start_sigma_e = 0.3)))
  list(fit = fit, data = dat)
}


test_that("k-fold pseudo CV returns a tibble with the requested scores", {
  o <- .fit_rspde_lme_1d(seed = 1L, n = 30L)
  cv <- posterior_crossvalidation(o$fit, mode = "k-fold", k = 3,
                                  scores = c("logscore", "crps", "mae", "rmse"),
                                  seed = 42)
  expect_named(cv, c("mu", "var", "scores"), ignore.order = TRUE)
  expect_length(cv$mu,  o$fit$nobs)
  expect_length(cv$var, o$fit$nobs)
  expect_true(all(is.finite(cv$mu)))
  expect_true(all(is.finite(cv$var)))
  expect_true(all(c("logscore", "crps", "mae", "rmse") %in% names(cv$scores)))
  expect_false("scrps" %in% names(cv$scores))
})


test_that("LOO via posterior_crossvalidation matches manual masked predict", {
  o <- .fit_rspde_lme_1d(seed = 1L, n = 30L)
  fit <- o$fit
  i <- 5L
  manual <- predict(fit,
                    newdata = o$data[i, , drop = FALSE],
                    loc = matrix(o$data$x[i], ncol = 1),
                    compute_variances = TRUE,
                    advanced_options = list(na_test_idx = i))
  sigma_e <- as.numeric(fit$coeff$measurement_error)
  cv <- posterior_crossvalidation(
    fit,
    train_test_indices = list(list(train = setdiff(seq_len(fit$nobs), i),
                                   test = i)),
    scores = c("logscore", "crps", "mae", "rmse"),
    tibble = FALSE)
  expect_equal(cv$mu[i], as.numeric(manual$mean), tolerance = 1e-10)
  expect_equal(cv$var[i], as.numeric(manual$variance) + sigma_e^2,
               tolerance = 1e-10)
})


test_that("LPO mode runs and respects number_folds", {
  o <- .fit_rspde_lme_1d(seed = 2L, n = 30L)
  cv <- posterior_crossvalidation(o$fit, mode = "lpo", percentage = 70,
                                  number_folds = 4, seed = 1,
                                  return_indices = TRUE)
  expect_length(cv$indices, 4L)
})


test_that("true_CV refits the model and returns finite scores", {
  o <- .fit_rspde_lme_1d(seed = 3L, n = 30L)
  cv <- suppressWarnings(
    posterior_crossvalidation(o$fit, mode = "k-fold", k = 3, true_CV = TRUE,
                              seed = 1))
  expect_true(all(is.finite(unlist(cv$scores))))
})


test_that("list-of-models input returns one row per model", {
  o <- .fit_rspde_lme_1d(seed = 4L, n = 30L)
  res <- posterior_crossvalidation(list(A = o$fit, B = o$fit),
                                   mode = "k-fold", k = 3, seed = 1)
  expect_true("Model" %in% names(res$scores))
  expect_equal(nrow(res$scores), 2L)
  expect_equal(res$scores$Model, c("A", "B"))
})


test_that("invalid scores are warned about and dropped", {
  o <- .fit_rspde_lme_1d(seed = 5L, n = 30L)
  expect_warning(
    cv <- posterior_crossvalidation(o$fit, mode = "k-fold", k = 3,
                                    scores = c("mae", "nonsense"),
                                    seed = 1),
    "Invalid scores")
  expect_true("mae" %in% names(cv$scores))
  expect_false("nonsense" %in% names(cv$scores))
})


# --- Graph-based tests ------------------------------------------------------

skip_if_no_metric_graph <- function() {
  testthat::skip_if_not_installed("MetricGraph")
}


.fit_graph_rspde_lme <- function(seed = 11L, n_per_edge = 6L, mesh_h = 0.05) {
  skip_if_no_metric_graph()
  set.seed(seed)
  V <- rbind(c(0, 0), c(1, 0), c(1, 1), c(0, 1))
  E <- rbind(c(1, 2), c(2, 3), c(3, 4), c(4, 1))
  graph <- MetricGraph::metric_graph$new(V = V, E = E, verbose = 0)
  graph$build_mesh(h = mesh_h)
  PtE <- NULL
  for (i in seq_len(graph$nE)) {
    PtE <- rbind(PtE, cbind(rep(i, n_per_edge), stats::runif(n_per_edge)))
  }
  u <- MetricGraph::sample_spde(kappa = 5, tau = 1, alpha = 1,
                                graph = graph, PtE = PtE)
  y <- u + 0.3 * stats::rnorm(length(u))
  graph$add_observations(
    data = data.frame(y = y,
                      edge_number = PtE[, 1],
                      distance_on_edge = PtE[, 2]),
    normalized = TRUE, verbose = 0)
  fit <- MetricGraph::graph_lme(y ~ -1, graph = graph,
                                model = list(type = "WhittleMatern",
                                             alpha = 1, fem = TRUE))
  class(fit) <- "rspde_lme"
  fit
}


test_that("pseudo-CV LOO for a pure-OLS fit equals X_full %*% beta_full", {
  set.seed(123)
  n <- 40L
  dat <- data.frame(elev = stats::rnorm(n), z = stats::rnorm(n))
  dat$y <- 0.5 + 1.5 * dat$elev - 0.7 * dat$z + 0.4 * stats::rnorm(n)
  fit <- rspde_lme(y ~ elev + z, data = dat)
  expect_true(isTRUE(fit$null_model))

  cv <- posterior_crossvalidation(fit, mode = "loo")
  X    <- stats::model.matrix(~ elev + z, data = dat)
  beta <- as.numeric(fit$coeff$fixed_effects)
  mu_manual <- as.numeric(X %*% beta)
  expect_equal(unname(cv$mu), unname(mu_manual), tolerance = 1e-12)
})


test_that("true_CV LOO for OLS matches refitting lm() on training rows", {
  set.seed(7)
  n <- 30L
  dat <- data.frame(elev = stats::rnorm(n))
  dat$y <- 1 - 0.4 * dat$elev + 0.5 * stats::rnorm(n)
  fit <- rspde_lme(y ~ elev, data = dat)

  cv <- posterior_crossvalidation(fit, mode = "loo", true_CV = TRUE)

  mu_lm <- numeric(n)
  for (i in seq_len(n)) {
    m <- stats::lm(y ~ elev, data = dat[-i, , drop = FALSE])
    mu_lm[i] <- as.numeric(stats::predict(m,
                                          newdata = dat[i, , drop = FALSE]))
  }
  expect_equal(unname(cv$mu), unname(mu_lm), tolerance = 1e-10)
})


test_that("pseudo-CV variance for OLS equals sigma^2 (1 + x'(X'X)^{-1}x)", {
  set.seed(42)
  n <- 25L
  dat <- data.frame(x1 = stats::rnorm(n), x2 = stats::rnorm(n))
  dat$y <- 0.2 + 0.6 * dat$x1 + 0.9 * dat$x2 + 0.3 * stats::rnorm(n)
  fit <- rspde_lme(y ~ x1 + x2, data = dat)

  X     <- stats::model.matrix(~ x1 + x2, data = dat)
  sigma <- as.numeric(fit$coeff$measurement_error)
  XtX_inv <- chol2inv(chol(crossprod(X)))
  expected_var <- sigma^2 * (1 + rowSums((X %*% XtX_inv) * X))

  cv <- posterior_crossvalidation(fit, mode = "loo",
                                  scores = c("logscore", "mae", "rmse"))
  expect_equal(unname(cv$var), unname(expected_var), tolerance = 1e-10)
})


test_that("k-fold CV on a 2D FEM alpha=2 fit matches a from-scratch kriging loop", {
  # Setup: simulate data on the unit square from spde.matern.operators with
  # alpha = 2 (integer alpha, standard FEM path), fit with rspde_lme. The
  # *manual reference* below does the kriging from scratch using only the
  # FEM mass / stiffness matrices and the fitted (tau, kappa, sigma_e,
  # beta) — it never calls predict.rspde_lme. The point is that a bug in
  # predict.rspde_lme would not be hidden by using predict on both sides
  # of the comparison.
  #
  # For alpha = 2 in 2D the SPDE precision is
  #     Q = tau^2 * (kappa^2 C + G)^T C^{-1} (kappa^2 C + G)
  # and a held-out kriging prediction at locations `loc` is
  #     mu_pred = X_pred beta + A_pred (A_obs^T A_obs/sigma_e^2 + Q)^{-1}
  #                                       A_obs^T y_centered / sigma_e^2
  #     var_pred = diag(A_pred (A_obs^T A_obs/sigma_e^2 + Q)^{-1} A_pred^T)
  #                  + sigma_e^2
  # where A is the FEM projection matrix and y_centered = y - X beta.

  set.seed(1)
  mesh <- fmesher::fm_mesh_2d(
    loc.domain = cbind(c(0, 1, 1, 0), c(0, 0, 1, 1)),
    max.edge   = c(0.08, 0.25),
    offset     = c(0.05, 0.20))
  fem_mats <- fmesher::fm_fem(mesh)
  Cmat <- fem_mats$c0
  Gmat <- fem_mats$g1
  Cinv <- Matrix::Diagonal(x = 1 / Matrix::diag(Cmat))

  op <- spde.matern.operators(kappa = 7, tau = 1, alpha = 2,
                              mesh = mesh,
                              parameterization = "spde",
                              type = "covariance")
  u <- as.vector(simulate(op, nsim = 1, seed = 1))
  n <- 80L
  loc <- cbind(stats::runif(n), stats::runif(n))
  A   <- fmesher::fm_basis(mesh, loc)              # n x nmesh
  sigma_e_true <- 0.3
  covar <- stats::rnorm(n)
  y <- 0.4 + 0.6 * covar + as.vector(A %*% u) +
       sigma_e_true * stats::rnorm(n)
  dat <- data.frame(y = y, x1 = loc[, 1], x2 = loc[, 2], covar = covar)

  # Fix alpha at 2 (the value used to simulate) so the random-effect set
  # is just (tau, kappa) and the manual Q construction below is exact.
  fit <- suppressWarnings(rspde_lme(
    y ~ covar, loc = c("x1", "x2"), data = dat, model = op,
    optim_controls = list(maxit = 25),
    model_options  = list(start_sigma_e = sigma_e_true,
                          fix_alpha     = 2)))

  # Pull the fitted parameters back out for the manual kriging.
  beta_hat    <- as.numeric(fit$coeff$fixed_effects)
  sigma_e_hat <- as.numeric(fit$coeff$measurement_error)
  tau_hat     <- as.numeric(fit$coeff$random_effects[["tau"]])
  kappa_hat   <- as.numeric(fit$coeff$random_effects[["kappa"]])

  # Sanity check that we know which parameterisation the fit reported.
  expect_true(all(c("tau", "kappa") %in%
                  names(fit$coeff$random_effects)))

  # Build Q manually and verify it agrees with op$Q at these parameters.
  L_op <- kappa_hat^2 * Cmat + Gmat
  Q_manual <- tau_hat^2 * Matrix::t(L_op) %*% Cinv %*% L_op
  op_hat <- update(op, kappa = kappa_hat, tau = tau_hat,
                   check_stationarity = FALSE)
  expect_equal(as.numeric(Matrix::norm(op_hat$Q - Q_manual, "F")), 0,
               tolerance = 1e-8)

  # Build a fixed k-fold partition used by BOTH the manual loop and
  # posterior_crossvalidation, so any difference is purely in the
  # scoring / kriging code, not in fold assignment.
  set.seed(11)
  k <- 5L
  fold_id <- sample(rep(seq_len(k), length.out = n))
  tt <- lapply(seq_len(k), function(i) {
    list(train = which(fold_id != i), test = which(fold_id == i))
  })

  # X design (including intercept) — these covariates do not vary with the
  # fold, so we build them once. y_centered is the same as well; the only
  # thing that changes per fold is the rows we keep in A and y_centered.
  X <- stats::model.matrix(~ covar, data = dat)
  y_centered_full <- as.numeric(dat$y - X %*% beta_hat)

  # Manual pseudo-CV loop --------------------------------------------------
  mu_manual  <- rep(NA_real_, n)
  var_manual <- rep(NA_real_, n)
  for (i in seq_len(k)) {
    ti <- tt[[i]]$test
    tr <- tt[[i]]$train

    A_train <- A[tr, , drop = FALSE]
    A_test  <- A[ti, , drop = FALSE]
    y_train <- y_centered_full[tr]

    Q_xgiveny <- Matrix::crossprod(A_train) / sigma_e_hat^2 + Q_manual
    rhs       <- Matrix::crossprod(A_train, y_train) / sigma_e_hat^2
    mu_krig   <- as.numeric(Matrix::solve(Q_xgiveny, rhs))

    # Posterior variance at the test locations:
    #   var = diag(A_test * Q_xgiveny^{-1} * A_test^T) + sigma_e^2.
    # Compute one column at a time via solve(Q_xgiveny, t(A_test)).
    post_cov <- A_test %*% Matrix::solve(Q_xgiveny, Matrix::t(A_test))
    var_post <- pmax(Matrix::diag(post_cov), 0)

    mu_manual[ti]  <- as.numeric(X[ti, , drop = FALSE] %*% beta_hat) +
                        as.numeric(A_test %*% mu_krig)
    var_manual[ti] <- as.numeric(var_post) + sigma_e_hat^2
  }

  sd_manual <- sqrt(pmax(var_manual, 0))
  resid     <- dat$y - mu_manual

  # CRPS / SCRPS helpers — identical formulas to the package internals,
  # but rewritten here so this reference does not depend on rSPDE's
  # implementation either.
  Efnorm <- function(mu, sigma) {
    sigma * sqrt(2 / pi) * exp(-(mu^2) / (2 * sigma^2)) +
      mu * (1 - 2 * stats::pnorm(-mu / sigma))
  }
  Exx <- function(mu, sigma) Efnorm(0, sqrt(2) * sigma)
  Exy <- function(mu, sigma, y) Efnorm(mu - y, sigma)
  CRPS_manual  <- -Exy(mu_manual, sd_manual, dat$y) +
                  0.5 * Exx(mu_manual, sd_manual)
  SCRPS_manual <- -Exy(mu_manual, sd_manual, dat$y) /
                  Exx(mu_manual, sd_manual) -
                  0.5 * log(Exx(mu_manual, sd_manual))

  expected <- list(
    logscore = -mean(stats::dnorm(dat$y, mu_manual, sd_manual, log = TRUE)),
    crps     = -mean(CRPS_manual),
    scrps    = -mean(SCRPS_manual),
    mae      =  mean(abs(resid)),
    rmse     =  sqrt(mean(resid^2))
  )

  # posterior_crossvalidation call -----------------------------------------
  cv <- posterior_crossvalidation(fit,
                                  train_test_indices = tt,
                                  scores = c("logscore", "crps", "scrps",
                                             "mae", "rmse"),
                                  tibble = FALSE)

  expect_equal(unname(cv$mu),  unname(mu_manual),  tolerance = 1e-8)
  expect_equal(unname(cv$var), unname(var_manual), tolerance = 1e-8)
  expect_equal(as.numeric(cv$scores$logscore), expected$logscore,
               tolerance = 1e-8)
  expect_equal(as.numeric(cv$scores$crps),  expected$crps,  tolerance = 1e-8)
  expect_equal(as.numeric(cv$scores$scrps), expected$scrps, tolerance = 1e-8)
  expect_equal(as.numeric(cv$scores$mae),   expected$mae,   tolerance = 1e-8)
  expect_equal(as.numeric(cv$scores$rmse),  expected$rmse,  tolerance = 1e-8)
})


test_that("graph-based posterior_crossvalidation matches MetricGraph's", {
  fit <- .fit_graph_rspde_lme(seed = 11L)
  # MetricGraph's posterior_crossvalidation needs the graph_lme class.
  fit_mg <- fit
  class(fit_mg) <- c("graph_lme", class(fit_mg))

  cv_mg    <- MetricGraph::posterior_crossvalidation(fit_mg, mode = "loo")
  cv_rspde <- posterior_crossvalidation(fit, mode = "loo")

  expect_equal(unname(cv_rspde$mu),  unname(cv_mg$mu),  tolerance = 1e-10)
  expect_equal(unname(cv_rspde$var), unname(cv_mg$var), tolerance = 1e-10)
  expect_equal(as.numeric(cv_rspde$scores$crps),
               as.numeric(cv_mg$scores$crps),  tolerance = 1e-10)
  expect_equal(as.numeric(cv_rspde$scores$rmse),
               as.numeric(cv_mg$scores$rmse),  tolerance = 1e-10)
})
