test_that("cross_validation routes exact graph fits when supported", {
  skip_if_not_installed("MetricGraph")
  fake_fit <- function(model_class, family = "gaussian") {
    structure(list(
      bru_info = list(model = list(effects = list(
        field = list(main = list(model = structure(list(), class = model_class)))
      ))),
      .args = list(family = family)
    ), class = "bru")
  }

  # Unsupported families stay in rSPDE, so this check does not need INLA.
  expect_error(rSPDE::cross_validation(fake_fit("inla_metric_graph_spde",
                                               "stochvol"), cv_type = "bogus"),
               "possible options for cv_type")
})

test_that("cross_validation delegates exact graph fits to MetricGraph", {
  skip_if_not_installed("MetricGraph")
  # MetricGraph::cross_validation() checks for INLA before validating cv_type.
  local_rspde_safe_inla()
  fake_fit <- function(model_class, family = "gaussian") {
    structure(list(
      bru_info = list(model = list(effects = list(
        field = list(main = list(model = structure(list(), class = model_class)))
      ))),
      .args = list(family = family)
    ), class = "bru")
  }

  graph_fit <- fake_fit("inla_metric_graph_spde")
  rspde_fit <- fake_fit("inla_rspde")
  # Validation wording identifies the delegated engine without requiring INLA fits.
  expect_error(rSPDE::cross_validation(graph_fit, cv_type = "bogus"),
               "cv_type must be one of")
  expect_error(rSPDE::cross_validation(list(rspde_fit, graph_fit),
                                       cv_type = "bogus"),
               "cv_type must be one of")
})
