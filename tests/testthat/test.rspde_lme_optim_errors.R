test_that("rspde_lme reports optim errors when all methods fail", {
  set.seed(1)
  x <- seq(0, 1, length.out = 6)
  data <- data.frame(
    y = rnorm(length(x)),
    x = x
  )

  model <- rSPDE::matern.operators(
    range = 1, sigma = 1, nu = 0.7,
    loc_mesh = x, d = 1,
    type = "operator",
    parameterization = "matern"
  )

  expect_error(
    rSPDE::rspde_lme(
      y ~ 1,
      loc = "x",
      data = data,
      model = model,
      optim_method = "NOT_A_METHOD",
      possible_methods = "NOT_A_METHOD",
      optim_controls = list(maxit = 1),
      model_options = list(start_sigma_e = 0.1)
    ),
    "Errors: NOT_A_METHOD:",
    fixed = TRUE
  )
})

test_that("rspde_lme stores likelihood issues on successful fit", {
  set.seed(2)
  x <- seq(0, 1, length.out = 6)
  data <- data.frame(
    y = rnorm(length(x)),
    x = x
  )

  model <- rSPDE::matern.operators(
    range = 1, sigma = 1, nu = 0.7,
    loc_mesh = x, d = 1,
    type = "operator",
    parameterization = "matern"
  )

  fit <- suppressWarnings(rSPDE::rspde_lme(
    y ~ 1,
    loc = "x",
    data = data,
    model = model,
    optim_controls = list(maxit = 0),
    model_options = list(start_sigma_e = 0.1)
  ))

  expect_true(is.list(fit$likelihood_issues))
})
