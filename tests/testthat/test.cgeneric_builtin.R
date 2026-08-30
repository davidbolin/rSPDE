test_that("built-in Cgeneric models retain the shared-library argument slot", {
  local_rspde_safe_inla()

  model <- INLA::inla.cgeneric.define(
    model = "inla_cgeneric_rspde_stat_int_model",
    shlib = NULL,
    n = 1L,
    parameterization = "spde"
  )

  expect_null(model$f$cgeneric$shlib)
  expect_false("shlib" %in% names(model$f$cgeneric$data$characters))

  model <- rspde_prepare_cgeneric_model(model)
  expect_identical(
    names(model$f$cgeneric$data$characters),
    c("model", "shlib", "parameterization")
  )
  expect_identical(model$f$cgeneric$data$characters$shlib, "")
  expect_silent(rspde_check_cgeneric_symbol(model))

  expect_identical(rspde_prepare_cgeneric_model(model), model)
})
