context("predict.rspde_lme: precompute and na_test_idx")

# These tests cover three behaviours added to predict.rspde_lme:
#   1. advanced_options$precompute_data returns / accepts cached
#      parameter-dependent structures (new_rspde_obj, Q, mu_corr) so a
#      cross-validation loop can skip the expensive update() + Q
#      assembly on every fold.
#   2. advanced_options$na_test_idx masks Y AND the corresponding rows
#      of A_list so leave-out predictions are *true* LOO, not
#      "predict at the test location while still using it as training".
#   3. $variance is appended across replicates, not overwritten with
#      the last replicate's values.

skip_if_no_metric_graph <- function() {
  skip_if_not_installed("MetricGraph")
}

.fit_rspde_lme <- function(seed = 1L, mesh_h = 0.05, n_per_edge = 6L) {
  skip_if_no_metric_graph()
  set.seed(seed)
  V <- rbind(c(0, 0), c(1, 0), c(1, 1), c(0, 1))
  E <- rbind(c(1, 2), c(2, 3), c(3, 4), c(4, 1))
  graph <- MetricGraph::metric_graph$new(V = V, E = E, verbose = 0)
  graph$build_mesh(h = mesh_h)
  PtE <- NULL
  for (i in seq_len(graph$nE)) {
    PtE <- rbind(PtE, cbind(rep(i, n_per_edge), runif(n_per_edge)))
  }
  u <- MetricGraph::sample_spde(kappa = 5, tau = 1, alpha = 1,
                                graph = graph, PtE = PtE)
  y <- u + 0.3 * rnorm(length(u))
  graph$add_observations(
    data = data.frame(y = y,
                      edge_number = PtE[, 1],
                      distance_on_edge = PtE[, 2]),
    normalized = TRUE, verbose = 0)
  fit <- MetricGraph::graph_lme(y ~ -1, graph = graph,
                                model = list(type = "WhittleMatern",
                                             alpha = 1, fem = TRUE))
  # These tests target predict.rspde_lme directly; strip the graph_lme
  # wrapper so we don't go through MetricGraph's dispatch (which may
  # not yet forward advanced_options in the user's installed version).
  class(fit) <- "rspde_lme"
  fit
}

test_that("precompute_data round-trips and accelerates a LOO loop", {
  fit <- .fit_rspde_lme(seed = 11L)
  gd <- fit$graph$.__enclos_env__$private$data
  nd_one <- lapply(gd, function(x) x[1])

  # The "store" call returns precomputed_data alongside the prediction.
  pred_store <- predict(fit, newdata = nd_one,
                        edge_number = ".edge_number",
                        distance_on_edge = ".distance_on_edge",
                        normalized = TRUE,
                        compute_variances = TRUE,
                        advanced_options = list(precompute_data = TRUE))
  expect_true(is.list(pred_store$precomputed_data))
  expect_named(pred_store$precomputed_data,
               c("new_rspde_obj", "Q", "mu_corr"), ignore.order = TRUE)

  # The "reuse" call must give the same prediction as a fresh call.
  pred_fresh <- predict(fit, newdata = nd_one,
                        edge_number = ".edge_number",
                        distance_on_edge = ".distance_on_edge",
                        normalized = TRUE, compute_variances = TRUE)
  pred_reuse <- predict(fit, newdata = nd_one,
                        edge_number = ".edge_number",
                        distance_on_edge = ".distance_on_edge",
                        normalized = TRUE,
                        compute_variances = TRUE,
                        advanced_options = list(
                          precompute_data = pred_store$precomputed_data))
  expect_equal(unname(pred_reuse$mean),     unname(pred_fresh$mean),
               tolerance = 1e-10)
  expect_equal(unname(pred_reuse$variance), unname(pred_fresh$variance),
               tolerance = 1e-10)
})

test_that("na_test_idx produces a true leave-out (masks both Y and A_list)", {
  fit <- .fit_rspde_lme(seed = 7L)
  gd <- fit$graph$.__enclos_env__$private$data
  i_test <- 5L
  nd <- lapply(gd, function(x) x[i_test])

  # Reference: do the masking manually outside predict.
  cv_model <- fit
  mm <- as.matrix(cv_model$model_matrix)
  mm[i_test, 1] <- NA
  cv_model$model_matrix <- mm
  cv_model$A_list <- fit$A_list
  cv_model$A_list[[1]] <- cv_model$A_list[[1]][-i_test, , drop = FALSE]
  pred_manual <- predict(cv_model, newdata = nd,
                         edge_number = ".edge_number",
                         distance_on_edge = ".distance_on_edge",
                         normalized = TRUE, compute_variances = TRUE)

  # New path: ask predict.rspde_lme to mask via na_test_idx.
  pred_new <- predict(fit, newdata = nd,
                      edge_number = ".edge_number",
                      distance_on_edge = ".distance_on_edge",
                      normalized = TRUE, compute_variances = TRUE,
                      advanced_options = list(na_test_idx = i_test))

  expect_equal(unname(pred_new$mean),     unname(pred_manual$mean),
               tolerance = 1e-10)
  expect_equal(unname(pred_new$variance), unname(pred_manual$variance),
               tolerance = 1e-10)

  # Sanity: a call without the mask uses the held-out point as training
  # and produces a noticeably smaller kriging variance.
  pred_leak <- predict(fit, newdata = nd,
                       edge_number = ".edge_number",
                       distance_on_edge = ".distance_on_edge",
                       normalized = TRUE, compute_variances = TRUE)
  expect_lt(pred_leak$variance, pred_new$variance)
})

test_that("variance is reported per replicate (not overwritten by the last loop iteration)", {
  skip_if_no_metric_graph()
  set.seed(2)
  V <- rbind(c(0, 0), c(1, 0), c(1, 1), c(0, 1))
  E <- rbind(c(1, 2), c(2, 3), c(3, 4), c(4, 1))
  graph <- MetricGraph::metric_graph$new(V = V, E = E, verbose = 0)
  graph$build_mesh(h = 0.1)
  PtE <- NULL
  for (i in seq_len(graph$nE)) {
    PtE <- rbind(PtE, cbind(rep(i, 4), runif(4)))
  }
  df_all <- NULL
  for (r in 1:2) {
    u <- MetricGraph::sample_spde(kappa = 5, tau = 1, alpha = 1,
                                  graph = graph, PtE = PtE)
    y <- u + 0.3 * rnorm(length(u))
    df_all <- rbind(df_all,
      data.frame(y = y, edge_number = PtE[, 1],
                 distance_on_edge = PtE[, 2], repl = r))
  }
  graph$add_observations(data = df_all, normalized = TRUE, verbose = 0,
                         group = "repl")
  fit <- MetricGraph::graph_lme(y ~ -1, graph = graph,
                                model = list(type = "WhittleMatern",
                                             alpha = 1, fem = TRUE))
  class(fit) <- "rspde_lme"  # see helper comment above

  # Mask one observation in replicate 1 (global index 5) and predict at
  # that same location. Replicate 1's kriging variance must be sizeable
  # (it's now LOO at that point); replicate 2's must be near zero (the
  # location is still observed there). The bug overwrote both with the
  # same value.
  gd <- fit$graph$.__enclos_env__$private$data
  nd <- list(
    .edge_number      = gd[[".edge_number"]][5],
    .distance_on_edge = gd[[".distance_on_edge"]][5],
    .group            = "1",
    y                 = NA_real_
  )
  pred <- predict(fit, newdata = nd,
                  edge_number = ".edge_number",
                  distance_on_edge = ".distance_on_edge",
                  normalized = TRUE, compute_variances = TRUE,
                  advanced_options = list(na_test_idx = 5L))

  v1 <- pred$variance[pred$repl == "1"]
  v2 <- pred$variance[pred$repl == "2"]
  expect_length(v1, 1L)
  expect_length(v2, 1L)
  expect_gt(v1, 1e-3)
  expect_gt(v1 / v2, 2)
})
