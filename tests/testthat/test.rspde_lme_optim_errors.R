test_that("rspde_lme reports optim errors when all methods fail", {
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

  ns <- asNamespace("rSPDE")
  old_optim <- get("optim", envir = ns, inherits = FALSE)
  old_create_likelihood <- get("create_likelihood", envir = ns, inherits = FALSE)

  unlockBinding("create_likelihood", ns)
  assign("create_likelihood", function(model, model_options, y_resp, X_cov, A_list,
                                       repl, start_values, mean_correction,
                                       smoothness_upper_bound, loc_df) {
    list(
      likelihood = function(theta) {
        warning("likelihood warn")
        stop("likelihood err")
      },
      estimate_params = seq_along(start_values),
      n_coeff_nonfixed = 0L
    )
  }, envir = ns)
  lockBinding("create_likelihood", ns)

  unlockBinding("optim", ns)
  assign("optim", function(par, fn, ...) {
    fn(par)
    stop("optim kaboom")
  }, envir = ns)
  lockBinding("optim", ns)
  on.exit({
    unlockBinding("optim", ns)
    assign("optim", old_optim, envir = ns)
    lockBinding("optim", ns)
    unlockBinding("create_likelihood", ns)
    assign("create_likelihood", old_create_likelihood, envir = ns)
    lockBinding("create_likelihood", ns)
  }, add = TRUE)

  expect_error(
    rspde_lme(
      y ~ 1,
      loc = "x",
      data = data,
      model = model,
      optim_method = "L-BFGS-B",
      possible_methods = "L-BFGS-B",
      optim_controls = list(maxit = 1),
      model_options = list(start_sigma_e = 0.1)
    ),
    "Errors: L-BFGS-B: optim kaboom",
    fixed = TRUE
  )

  expect_error(
    rspde_lme(
      y ~ 1,
      loc = "x",
      data = data,
      model = model,
      optim_method = "L-BFGS-B",
      possible_methods = "L-BFGS-B",
      optim_controls = list(maxit = 1),
      model_options = list(start_sigma_e = 0.1)
    ),
    "Last likelihood issues: warning: likelihood warn \\| error: likelihood err"
  )
})

test_that("rspde_lme stores likelihood issues on successful fit", {
  set.seed(2)
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

  expect_true(is.list(fit$likelihood_issues))
})
