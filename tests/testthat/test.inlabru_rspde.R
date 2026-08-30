context("inlabru_rspde")

test_that("Inlabru predict method works for rspde.matern1d with nu = 0.7", {
    testthat::skip_on_cran()
    local_rspde_safe_inla()
    skip_if_not_installed("inlabru")

    set.seed(1)
    library(INLA)
    library(inlabru)

    # Create observations at two locations (0 and 1)
    # Use fixed values instead of sampling
    y <- c(2.0, 4.0) # Observations at locations 0 and 1
    obs_loc <- c(0, 1)

    # Create data frame
    df_bru <- data.frame(y = y, loc = obs_loc)

    # Create SPDE model for the interval [0, 1]
    bru_model <- rspde.matern1d(loc = obs_loc, nu = 0.7)

    # Create component
    cmp <- y ~ -1 + field(loc, model = bru_model)

    # Fit the model
    bru_fit <- bru(cmp,
        data = df_bru,
        options = list(
            num.threads = "1:1",
            verbose = FALSE,
            control.inla = list(int.strategy = "eb")
        )
    )

    # Predict at midpoint (location 0.5)
    pred_loc <- data.frame(loc = 0.5)

    field_pred <- predict(bru_model, cmp, bru_fit,
        newdata = pred_loc,
        formula = ~field
    )

    # Get predicted mean
    pred_mean <- field_pred$mean

    # Get average of observations
    obs_mean <- mean(y)

    # Check that prediction is close to average (with reasonable tolerance)
    expect_equal(pred_mean, obs_mean, tolerance = 0.5)
})


test_that("Inlabru predict method works for rspde.matern1d with Intercept and nu = 1.1", {
    testthat::skip_on_cran()
    local_rspde_safe_inla()
    skip_if_not_installed("inlabru")

    set.seed(2)
    library(INLA)
    library(inlabru)

    # Create observations at two locations (0 and 1)
    # Use fixed values instead of sampling
    y <- c(2, 2.2) # Observations at locations 0 and 1
    obs_loc <- c(0, 1)

    # Create data frame
    df_bru <- data.frame(y = y, loc = obs_loc)

    # Create SPDE model for the interval [0, 1]
    bru_model <- rspde.matern1d(loc = obs_loc, nu = 1.1)

    # Create component without intercept
    cmp <- y ~ -1 + Intercept(1) + field(loc, model = bru_model)

    # Fit the model
    bru_fit <- bru(cmp,
        data = df_bru,
        options = list(
            num.threads = "1:1",
            verbose = FALSE
        )
    )

    # Predict at midpoint (location 0.5)
    pred_loc <- data.frame(loc = 0.5)

    field_pred <- predict(bru_model, cmp, bru_fit,
        newdata = pred_loc,
        formula = ~ Intercept + field
    )

    # Get predicted mean
    pred_mean <- field_pred$mean

    # Get average of observations
    obs_mean <- mean(y)

    # Check that prediction is close to average (with reasonable tolerance)
    expect_equal(pred_mean, obs_mean, tolerance = 0.5)
})
