## posterior_crossvalidation() is an S3 generic so that other packages
## (MetricGraph's graph_lme method) can share it. These tests use stand-in
## classes so that the dispatch does not depend on MetricGraph.

test_that("a list of models dispatches on each element", {
  registerS3method("posterior_crossvalidation", "fake_cv_fit",
                   function(object, ...) {
                     list(mu = object$mu, var = 1,
                          scores = data.frame(mae = object$mu,
                                              rmse = 2 * object$mu))
                   },
                   envir = asNamespace("rSPDE"))
  fits <- list(A = structure(list(mu = 1), class = "fake_cv_fit"),
               B = structure(list(mu = 3), class = "fake_cv_fit"))

  res <- posterior_crossvalidation(fits)
  expect_equal(res$scores$Model, c("A", "B"))
  expect_equal(res$scores$mae, c(1, 3))
  expect_equal(res$scores$rmse, c(2, 6))
  expect_equal(res$mu, list(A = 1, B = 3))
})

test_that("a list of models only forwards the arguments that were supplied", {
  seen <- NULL
  registerS3method("posterior_crossvalidation", "fake_cv_args",
                   function(object, ...) {
                     seen <<- names(list(...))
                     list(mu = 0, var = 0, scores = data.frame(mae = 0))
                   },
                   envir = asNamespace("rSPDE"))

  posterior_crossvalidation(list(a = structure(list(), class = "fake_cv_args")),
                            k = 5)
  expect_setequal(seen, c("k", "tibble", "return_indices"))
})

test_that("unsupported objects give an informative error", {
  expect_error(posterior_crossvalidation(1),
               "No posterior_crossvalidation\\(\\) method")
  expect_error(posterior_crossvalidation(list(a = 1)),
               "must be a fitted model")
})
