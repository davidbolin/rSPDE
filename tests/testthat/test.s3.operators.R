context("S3 operator methods")

test_that("S3 methods for matern operators work on small 1d meshes", {
  x <- seq(0, 1, length.out = 6)
  loc <- matrix(c(0.2, 0.4), ncol = 1)

  op_operator <- matern.operators(
    range = 1, sigma = 1, nu = 0.7,
    loc_mesh = x, d = 1,
    type = "operator",
    parameterization = "matern"
  )

  A1 <- make_A(op_operator, loc = loc)
  expect_equal(nrow(A1), nrow(loc))
  expect_equal(ncol(A1), length(x))

  cov_vec <- cov_function_mesh(op_operator, p = matrix(0.2, ncol = 1))
  expect_equal(length(cov_vec), length(x))

  cov_mat <- covariance_mesh(op_operator)
  expect_equal(dim(cov_mat), c(length(x), length(x)))

  op_cov <- matern.operators(
    range = 1, sigma = 1, nu = 0.7,
    loc_mesh = x, d = 1,
    type = "covariance",
    parameterization = "matern"
  )

  A2 <- make_A(op_cov, loc = loc)
  expect_equal(nrow(A2), nrow(loc))
  expect_equal(ncol(A2), length(x))

  cov_vec2 <- cov_function_mesh(op_cov, p = matrix(0.2, ncol = 1))
  expect_equal(dim(cov_vec2), c(length(x), 1))

  cov_mat2 <- covariance_mesh(op_cov)
  expect_equal(dim(cov_mat2), c(length(x), length(x)))
})

test_that("S3 methods for spde matern operators work on small 1d meshes", {
  x <- seq(0, 1, length.out = 6)
  loc <- matrix(c(0.3, 0.6), ncol = 1)

  op_spde <- spde.matern.operators(
    kappa = 2, tau = 1, alpha = 1.2,
    loc_mesh = x, d = 1,
    type = "operator",
    parameterization = "spde"
  )

  A <- make_A(op_spde, loc = loc)
  expect_equal(nrow(A), nrow(loc))
  expect_equal(ncol(A), length(x))

  cov_vec <- cov_function_mesh(op_spde, p = matrix(0.3, ncol = 1))
  expect_equal(length(cov_vec), length(x))

  cov_mat <- covariance_mesh(op_spde)
  expect_equal(dim(cov_mat), c(length(x), length(x)))
})

test_that("S3 methods for matern2d operators work on small 2d meshes", {
  locs <- matrix(
    c(0, 0,
      1, 0,
      0, 1,
      1, 1),
    ncol = 2,
    byrow = TRUE
  )
  mesh <- fmesher::fm_mesh_2d(locs, max.edge = 2)

  op2d <- matern2d.operators(
    mesh = mesh, sigma = 1, nu = 1,
    hx = 0.2, hy = 0.2, hxy = 0,
    m = 1
  )

  A <- make_A(op2d, loc = matrix(c(0.2, 0.2), ncol = 2))
  expect_equal(nrow(A), 1)
  expect_true(ncol(A) > 0)

  cov_vec <- cov_function_mesh(op2d, p = matrix(c(0.2, 0.2), ncol = 2))
  expect_equal(dim(cov_vec)[1], ncol(A))

  cov_mat <- covariance_mesh(op2d)
  expect_equal(nrow(cov_mat), ncol(A))
  expect_equal(ncol(cov_mat), ncol(A))
})

test_that("S3 make_A for intrinsic operators works on small 1d meshes", {
  x <- seq(0, 1, length.out = 6)
  op_intrinsic <- intrinsic.matern.operators(
    kappa = 1, tau = 1, alpha = 1.2, beta = 1,
    loc_mesh = x, d = 1, m_alpha = 1, m_beta = 1
  )

  A <- make_A(op_intrinsic, loc = matrix(0.3, ncol = 1))
  expect_equal(nrow(A), 1)
  expect_equal(ncol(A), length(x) * op_intrinsic$m)
})

test_that("S3 make_A for spacetime operators works on small meshes", {
  s <- seq(0, 1, length.out = 4)
  t <- seq(0, 1, length.out = 3)
  op_st <- spacetime.operators(
    space_loc = s, time_loc = t,
    kappa = 1, sigma = 1, gamma = 0.1,
    rho = 0, alpha = 1, beta = 1
  )

  A <- make_A(op_st, loc = matrix(c(0.2, 0.7), ncol = 1), time = c(0.2, 0.8))
  expect_equal(nrow(A), 2)
  expect_true(ncol(A) > 0)
})

test_that("rspde_lme uses make_A S3 methods with a small model", {
  set.seed(1)
  x <- seq(0, 1, length.out = 6)
  data <- data.frame(
    y = rnorm(length(x)),
    x = x
  )

  model <- matern.operators(
    range = 1, sigma = 1, nu = 0.7,
    loc_mesh = x, d = 1,
    type = "operator",
    parameterization = "matern"
  )

  fit <- rspde_lme(
    y ~ 1,
    loc = "x",
    data = data,
    model = model,
    optim_controls = list(maxit = 0),
    model_options = list(
      fix_range = 1,
      fix_sigma = 1,
      fix_nu = 0.7,
      start_sigma_e = 0.1
    )
  )

  expect_true(inherits(fit, "rspde_lme"))
  expect_true(inherits(fit$latent_model, "matern_operator"))
})

test_that("spde.matern.operators delegates to matern.operators for constant parameters", {
  x <- seq(0, 1, length.out = 6)
  tau <- 1
  kappa <- 2
  alpha <- 1.2

  op_spde <- spde.matern.operators(
    kappa = kappa, tau = tau, alpha = alpha,
    loc_mesh = x, d = 1,
    type = "operator",
    parameterization = "spde"
  )

  expect_true(inherits(op_spde, "matern_operator"))
  expect_equal(op_spde$kappa, kappa)
  expect_equal(op_spde$tau, tau)
})

test_that("spde.matern.operators computes tau/kappa from theta with spde parameterization", {
  x <- seq(0, 1, length.out = 6)
  B.tau <- matrix(c(log(2), 0, 0), 1, 3)
  B.kappa <- matrix(c(log(3), 0, 0), 1, 3)

  theta <- c(0, 0)
  tau_exp <- as.numeric(exp(B.tau %*% c(1, theta)))
  kappa_exp <- as.numeric(exp(B.kappa %*% c(1, theta)))
  op_theta <- spde.matern.operators(
    theta = theta,
    B.tau = B.tau,
    B.kappa = B.kappa,
    alpha = 1.2,
    loc_mesh = x, d = 1,
    type = "operator",
    parameterization = "spde"
  )

  expect_true(inherits(op_theta, "matern_operator"))
  expect_equal(op_theta$tau, tau_exp)
  expect_equal(op_theta$kappa, kappa_exp)
})

test_that("spde.matern.operators computes tau/kappa from theta with matern parameterization", {
  x <- seq(0, 1, length.out = 6)
  B.sigma <- matrix(c(log(1.5), 0, 0), 1, 3)
  B.range <- matrix(c(log(2.5), 0, 0), 1, 3)

  theta <- c(0, 0)
  sigma_exp <- as.numeric(exp(B.sigma %*% c(1, theta)))
  range_exp <- as.numeric(exp(B.range %*% c(1, theta)))
  op_theta <- spde.matern.operators(
    theta = theta,
    B.sigma = B.sigma,
    B.range = B.range,
    nu = 0.7,
    loc_mesh = x, d = 1,
    type = "operator",
    parameterization = "matern"
  )

  expect_true(inherits(op_theta, "matern_operator"))
  expect_equal(op_theta$parameterization, "matern")
  expect_equal(op_theta$range, range_exp)
  expect_equal(op_theta$sigma, sigma_exp)
  expect_equal(op_theta$kappa, sqrt(8 * 0.7) / range_exp)
})
