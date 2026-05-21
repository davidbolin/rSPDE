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

test_that("rspde_lme falls back to OLS when no model is given", {
  # Regression test: when called without a model argument, rspde_lme
  # should do a plain OLS regression. Several variables (`loc_df`,
  # `coeff_alt_par_result`, ...) are only defined in the random-effect
  # branch; they must be initialised to NULL in the OLS branch so that
  # the final object construction works.
  set.seed(1)
  n <- 30
  data <- data.frame(elev = rnorm(n))
  data$y <- 0.5 + 1.5 * data$elev + 0.3 * rnorm(n)

  fit <- rspde_lme(y ~ elev, data = data)
  expect_s3_class(fit, "rspde_lme")
  expect_true(fit$null_model)
  expect_null(fit$latent_model)
  expect_named(fit$coeff$fixed_effects, c("(Intercept)", "elev"))
  expect_equal(unname(fit$coeff$fixed_effects[2]), 1.5, tolerance = 0.2)

  # Print and summary should describe a linear regression, not a
  # latent SPDE model.
  printed <- capture.output(print(fit))
  expect_true(any(grepl("Linear regression model", printed, fixed = TRUE)))
  summ <- capture.output(print(summary(fit)))
  expect_true(any(grepl("Linear regression model", summ, fixed = TRUE)))
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
