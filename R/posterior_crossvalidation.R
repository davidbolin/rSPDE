#' @name posterior_crossvalidation
#' @title Posterior cross-validation for \code{rspde_lme} models
#' @description Performs cross-validation on objects fitted with
#'   \code{\link{rspde_lme}} and computes common predictive scores. The interface
#'   mirrors \code{MetricGraph::posterior_crossvalidation} so that the two can be
#'   used interchangeably for graph-based and non-graph rSPDE fits. Pure OLS
#'   fits (\code{rspde_lme(formula, data)} with no \code{model}) are also
#'   supported: predictions reduce to evaluating the covariates at the test
#'   locations, with predictive variance
#'   \eqn{\sigma_\epsilon^2 \, (1 + x^\top (X^\top X)^{-1} x)}.
#'
#'   For \code{true_CV = FALSE} (default, "pseudo" cross-validation) the model
#'   parameters are kept fixed at the estimates from the full fit and only the
#'   held-out points are masked at prediction time. When \code{use_precomputed
#'   = TRUE} the parameter-dependent structures (updated rSPDE operator and the
#'   precision matrix Q) are computed once and reused across folds, mirroring
#'   the \code{advanced_options$precompute_data} / \code{na_test_idx} path of
#'   \code{predict.rspde_lme}.
#'
#'   For \code{true_CV = TRUE} the model is refit on each training fold (with
#'   the held-out response values set to \code{NA}). The previous fit is
#'   forwarded via \code{previous_fit} so that the optimisation starts from the
#'   full-data estimates.
#'
#' @param object A fitted object of class \code{rspde_lme}, or a (preferably
#'   named) list of such objects. When a list is supplied, the function is
#'   applied to each element and the scores are returned in a single
#'   \code{data.frame} / \code{tibble}.
#' @param scores Character vector of scores to compute. Possible values are
#'   \code{"logscore"}, \code{"crps"}, \code{"scrps"}, \code{"mae"} and
#'   \code{"rmse"}.
#' @param mode One of \code{"k-fold"}, \code{"loo"} or \code{"lpo"}.
#' @param k Number of folds for k-fold cross-validation. Default is 10.
#' @param percentage Percentage (1-99) of observations used for training in
#'   leave-percentage-out (\code{"lpo"}) cross-validation. Default is 20.
#' @param number_folds Number of folds for \code{"lpo"}. Default is 10.
#' @param train_test_indices Optional pre-specified list of folds. Each element
#'   must be a list with integer vectors \code{train} and \code{test} giving
#'   positions into the fitted observations (length \code{object$nobs}). When
#'   supplied, \code{mode}, \code{k}, \code{percentage} and \code{number_folds}
#'   are ignored.
#' @param true_CV If \code{TRUE} the model is refit on every training fold;
#'   if \code{FALSE} (default) the original fitted parameters are reused
#'   (pseudo cross-validation).
#' @param factor Multiplier applied to the mean of each score (default 1).
#' @param tibble If \code{TRUE} (default), the returned scores are coerced to
#'   a \code{tibble} (requires the \pkg{tidyr} package).
#' @param parallel_folds Process folds in parallel. Default \code{FALSE}.
#' @param parallel_fitting Run the per-fold model refit in parallel (only
#'   relevant when \code{true_CV = TRUE}). Default \code{FALSE}.
#'   \code{parallel_folds} and \code{parallel_fitting} cannot both be
#'   \code{TRUE}.
#' @param n_cores Number of cores for parallel computation. Default
#'   \code{parallel::detectCores() - 1}.
#' @param print Print progress messages.
#' @param seed Random seed for fold creation reproducibility.
#' @param return_indices If \code{TRUE} the train/test indices used in each
#'   fold are also returned.
#' @param use_precomputed If \code{TRUE} (default) the parameter-dependent
#'   structures used by \code{predict.rspde_lme} are precomputed once and
#'   reused across folds. Only relevant when \code{true_CV = FALSE}.
#' @param data Optional \code{data.frame} containing the original data used to
#'   fit the model. Required for non-graph fits with non-trivial covariates
#'   (e.g. factors or transformations), where the design matrix stored in
#'   \code{object$model_matrix} is insufficient to rebuild the covariates at
#'   the test locations. Ignored for graph-based fits, where the data is taken
#'   from the metric graph.
#' @param nelder_mead_init Logical, forwarded to \code{\link{rspde_lme}} when
#'   \code{true_CV = TRUE}. Default is \code{FALSE}: per-fold refits warm-start
#'   from the full-fit parameters via \code{previous_fit}, so the Nelder-Mead
#'   pre-pass adds little except cost. Set to \code{TRUE} to re-enable it for
#'   robustness if the full fit is itself suspected of being a poor local
#'   optimum.
#'
#' @return A list with elements
#'   \describe{
#'     \item{mu}{vector of posterior means, one per observation}
#'     \item{var}{vector of posterior variances (including measurement
#'       error), one per observation}
#'     \item{scores}{a \code{data.frame} (or \code{tibble}) with the requested
#'       aggregated scores}
#'     \item{indices}{list of train/test indices used (if
#'       \code{return_indices = TRUE})}
#'   }
#' @export
posterior_crossvalidation <- function(object,
                                      scores = c("logscore", "crps", "scrps",
                                                 "mae", "rmse"),
                                      mode = "k-fold",
                                      k = 10,
                                      percentage = 20,
                                      number_folds = 10,
                                      train_test_indices = NULL,
                                      true_CV = FALSE,
                                      factor = 1,
                                      tibble = TRUE,
                                      parallel_folds = FALSE,
                                      parallel_fitting = FALSE,
                                      n_cores = parallel::detectCores() - 1,
                                      print = FALSE,
                                      seed = NULL,
                                      return_indices = FALSE,
                                      use_precomputed = TRUE,
                                      data = NULL,
                                      nelder_mead_init = FALSE) {

  if (!inherits(object, "rspde_lme") && !is.list(object)) {
    stop("object should be of class 'rspde_lme' or a list of 'rspde_lme' objects.")
  }

  scores <- tolower(scores)
  valid_scores <- c("logscore", "crps", "scrps", "mae", "rmse")
  invalid_scores <- setdiff(scores, valid_scores)
  if (length(invalid_scores) > 0) {
    warning(paste("Invalid scores:", paste(invalid_scores, collapse = ", "),
                  "- will be ignored. Valid options are:",
                  paste(valid_scores, collapse = ", ")))
    scores <- intersect(scores, valid_scores)
  }
  if (length(scores) == 0) {
    stop("No valid scores requested.")
  }

  if (isTRUE(parallel_folds) && isTRUE(parallel_fitting)) {
    stop("'parallel_folds' and 'parallel_fitting' cannot both be TRUE.")
  }

  if (!inherits(object, "rspde_lme")) {
    is_rspde_list <- vapply(object, inherits, logical(1), what = "rspde_lme")
    if (!all(is_rspde_list)) {
      stop("All elements of 'object' must be of class 'rspde_lme'.")
    }
    if (is.null(names(object))) {
      warning("The list of fitted models has no names; results will not be properly named.")
      names(object) <- paste0("Model ", seq_along(object))
    }

    results_list <- lapply(object, function(obj) {
      posterior_crossvalidation(obj,
                                scores = scores, mode = mode, k = k,
                                percentage = percentage,
                                number_folds = number_folds,
                                train_test_indices = train_test_indices,
                                true_CV = true_CV, factor = factor,
                                tibble = FALSE,
                                parallel_folds = parallel_folds,
                                parallel_fitting = parallel_fitting,
                                n_cores = n_cores, print = print,
                                seed = seed,
                                return_indices = return_indices,
                                use_precomputed = use_precomputed,
                                data = data,
                                nelder_mead_init = nelder_mead_init)
    })

    res <- list()
    res[["mu"]]  <- lapply(results_list, `[[`, "mu")
    res[["var"]] <- lapply(results_list, `[[`, "var")
    scores_df <- do.call(rbind, lapply(results_list, `[[`, "scores"))
    rownames(scores_df) <- names(results_list)
    if (isTRUE(tibble) && requireNamespace("tidyr", quietly = TRUE)) {
      scores_df[["Model"]] <- rownames(scores_df)
      score_cols <- intersect(c("logscore", "crps", "scrps", "mae", "rmse"),
                              names(scores_df))
      scores_df <- tidyr::as_tibble(scores_df[, c("Model", score_cols),
                                              drop = FALSE])
    }
    res[["scores"]] <- scores_df
    if (return_indices && !is.null(results_list[[1]][["indices"]])) {
      res[["indices"]] <- results_list[[1]][["indices"]]
    }
    return(res)
  }

  cv_info <- .rspde_cv_build_data(object, data)
  n_obs <- cv_info$n_obs
  repl_vec <- cv_info$repl_vec
  y_vec    <- cv_info$y_vec

  if (is.null(train_test_indices)) {
    if (!is.null(seed)) {
      set.seed(seed)
    }
    if (mode == "loo") {
      train_test_indices <- lapply(seq_len(n_obs), function(i) {
        list(train = setdiff(seq_len(n_obs), i), test = i)
      })
      k <- n_obs
    } else if (mode == "k-fold") {
      fold_indices <- sample(rep(seq_len(k), length.out = n_obs))
      train_test_indices <- lapply(seq_len(k), function(i) {
        test <- which(fold_indices == i)
        list(train = setdiff(seq_len(n_obs), test), test = test)
      })
    } else if (mode == "lpo") {
      if (percentage <= 0 || percentage >= 100) {
        stop("'percentage' must be a number between 1 and 99.")
      }
      n_train <- round(n_obs * percentage / 100)
      train_test_indices <- lapply(seq_len(number_folds), function(i) {
        train <- sample(seq_len(n_obs), n_train)
        test  <- setdiff(seq_len(n_obs), train)
        list(train = train, test = test)
      })
      k <- number_folds
    } else {
      stop("Invalid 'mode'. Choose 'k-fold', 'loo' or 'lpo'.")
    }
  } else {
    k <- length(train_test_indices)
  }

  precomputed_data <- NULL
  need_variances <- any(scores %in% c("logscore", "crps", "scrps"))

  if (isTRUE(use_precomputed) && !isTRUE(true_CV) && !isTRUE(cv_info$is_null)) {
    if (isTRUE(print)) {
      cat("Pre-computing parameter-dependent structures...\n")
    }
    dummy_idx <- train_test_indices[[1]]$test[1]
    pred_args <- .rspde_cv_predict_args(cv_info, dummy_idx,
                                        compute_variances = need_variances,
                                        which_repl = repl_vec[dummy_idx])
    pred_args$advanced_options <- list(precompute_data = TRUE)
    precomp_pred <- tryCatch(do.call(stats::predict, c(list(object), pred_args)),
                             error = function(e) e)
    if (inherits(precomp_pred, "error")) {
      warning("Precomputation failed; falling back to standard prediction. ",
              "Reason: ", conditionMessage(precomp_pred))
    } else {
      precomputed_data <- precomp_pred$precomputed_data
      if (is.null(precomputed_data)) {
        warning("Precomputation returned no structures; falling back to standard prediction.")
      } else if (isTRUE(print)) {
        cat("Precomputation successful.\n")
      }
    }
  }

  process_fold <- function(fold_idx) {
    if (isTRUE(print) && !isTRUE(parallel_folds)) {
      cat("Processing fold", fold_idx, "of", k, "\n")
    }

    train_indices <- train_test_indices[[fold_idx]]$train
    test_indices  <- train_test_indices[[fold_idx]]$test
    n_test <- length(test_indices)

    local_results <- list(
      test_indices = test_indices,
      mu.p     = rep(NA_real_, n_test),
      var.p    = rep(NA_real_, n_test),
      logscore = rep(NA_real_, n_test),
      crps     = rep(NA_real_, n_test),
      scrps    = rep(NA_real_, n_test),
      mae      = rep(NA_real_, n_test),
      rmse     = rep(NA_real_, n_test)
    )

    # OLS (null-model) fast path: no latent kriging, just X %*% beta. 
    if (isTRUE(cv_info$is_null)) {
      ols_res <- .rspde_cv_ols_fold(cv_info,
                                    train_indices = train_indices,
                                    test_indices  = test_indices,
                                    true_CV       = true_CV)
      mu_vec  <- ols_res$mu
      var_vec <- ols_res$var

      local_results$mu.p <- mu_vec
      if (need_variances) {
        local_results$var.p <- var_vec
        sd_vec <- sqrt(pmax(var_vec, 0))
      }

      y_test_vec <- y_vec[test_indices]
      resid_vec  <- y_test_vec - mu_vec
      if ("logscore" %in% scores) {
        local_results$logscore <- stats::dnorm(y_test_vec, mean = mu_vec,
                                               sd = sd_vec, log = TRUE)
      }
      if ("crps" %in% scores) {
        local_results$crps <- .rspde_cv_CRPS(y_test_vec, mu_vec, sd_vec)
      }
      if ("scrps" %in% scores) {
        local_results$scrps <- .rspde_cv_SCRPS(y_test_vec, mu_vec, sd_vec)
      }
      if ("mae" %in% scores) {
        local_results$mae <- abs(resid_vec)
      }
      if ("rmse" %in% scores) {
        local_results$rmse <- resid_vec^2
      }
      return(local_results)
    }

    if (isTRUE(true_CV)) {
      cv_model <- .rspde_cv_refit(object, cv_info, test_indices,
                                  parallel_fitting = parallel_fitting,
                                  n_cores = n_cores,
                                  nelder_mead_init = nelder_mead_init)
    } else {
      cv_model <- object
    }

    test_repls <- repl_vec[test_indices]
    unique_test_repls <- unique(test_repls)
    mu_vec  <- rep(NA_real_, n_test)
    var_vec <- rep(NA_real_, n_test)

    for (r in unique_test_repls) {
      pos_r <- which(test_repls == r)
      ti_r  <- test_indices[pos_r]

      pred_args <- .rspde_cv_predict_args(cv_info, ti_r,
                                          compute_variances = need_variances,
                                          which_repl = r)
      adv <- list(na_test_idx = test_indices)
      if (!isTRUE(true_CV) && !is.null(precomputed_data)) {
        adv$precompute_data <- precomputed_data
      }
      pred_args$advanced_options <- adv

      pred_r <- do.call(stats::predict, c(list(cv_model), pred_args))

      mu_vec[pos_r] <- as.numeric(pred_r$mean)
      if (need_variances) {
        var_vec[pos_r] <- as.numeric(pred_r$variance)
      }
    }

    local_results$mu.p <- mu_vec
    if (need_variances) {
      sigma_e <- as.numeric(cv_model$coeff$measurement_error[[1]])
      var_total <- var_vec + sigma_e^2
      local_results$var.p <- var_total
      sd_vec <- sqrt(pmax(var_total, 0))
    }

    y_test_vec <- y_vec[test_indices]
    resid_vec  <- y_test_vec - mu_vec

    if ("logscore" %in% scores) {
      local_results$logscore <- stats::dnorm(y_test_vec, mean = mu_vec,
                                             sd = sd_vec, log = TRUE)
    }
    if ("crps" %in% scores) {
      local_results$crps <- .rspde_cv_CRPS(y_test_vec, mu_vec, sd_vec)
    }
    if ("scrps" %in% scores) {
      local_results$scrps <- .rspde_cv_SCRPS(y_test_vec, mu_vec, sd_vec)
    }
    if ("mae" %in% scores) {
      local_results$mae <- abs(resid_vec)
    }
    if ("rmse" %in% scores) {
      local_results$rmse <- resid_vec^2
    }

    local_results
  }

  if (isTRUE(parallel_folds)) {
    if (!requireNamespace("parallel", quietly = TRUE)) {
      warning("Package 'parallel' is not available. Using sequential processing instead.")
      parallel_folds <- FALSE
    }
  }

  if (isTRUE(parallel_folds)) {
    cl <- parallel::makeCluster(n_cores)
    on.exit(try(parallel::stopCluster(cl), silent = TRUE), add = TRUE)
    parallel::clusterEvalQ(cl, suppressPackageStartupMessages(library(Matrix)))
    parallel::clusterExport(cl,
                            varlist = c("object", "cv_info",
                                        "train_test_indices",
                                        "precomputed_data", "need_variances",
                                        "scores", "true_CV", "repl_vec",
                                        "y_vec", "k", "print", "parallel_folds",
                                        "parallel_fitting", "n_cores"),
                            envir = environment())
    results <- parallel::parLapply(cl, seq_len(k), process_fold)
  } else {
    results <- lapply(seq_len(k), process_fold)
  }

  # combine results across folds 
  mu.p <- var.p <- rep(NA_real_, n_obs)
  logscore <- crps <- scrps <- mae <- rmse <- rep(NA_real_, n_obs)

  for (i in seq_len(k)) {
    test_indices <- results[[i]]$test_indices
    mu.p[test_indices]     <- results[[i]]$mu.p
    var.p[test_indices]    <- results[[i]]$var.p
    logscore[test_indices] <- results[[i]]$logscore
    crps[test_indices]     <- results[[i]]$crps
    scrps[test_indices]    <- results[[i]]$scrps
    mae[test_indices]      <- results[[i]]$mae
    rmse[test_indices]     <- results[[i]]$rmse
  }

  res <- list(mu = mu.p, var = var.p)
  res[["scores"]] <- data.frame(dummy = NA)
  if ("logscore" %in% scores) {
    res[["scores"]]$logscore <- -factor * mean(logscore, na.rm = TRUE)
  }
  if ("crps" %in% scores) {
    res[["scores"]]$crps <- -factor * mean(crps, na.rm = TRUE)
  }
  if ("scrps" %in% scores) {
    res[["scores"]]$scrps <- -factor * mean(scrps, na.rm = TRUE)
  }
  if ("mae" %in% scores) {
    res[["scores"]]$mae <- factor * mean(mae, na.rm = TRUE)
  }
  if ("rmse" %in% scores) {
    res[["scores"]]$rmse <- factor * sqrt(mean(rmse, na.rm = TRUE))
  }
  res[["scores"]]$dummy <- NULL
  attr(res[["scores"]], "factor") <- factor

  if (return_indices) {
    res[["indices"]] <- train_test_indices
  }

  if (isTRUE(tibble) && requireNamespace("tidyr", quietly = TRUE)) {
    res[["scores"]] <- tidyr::as_tibble(res[["scores"]])
  }

  res
}

# Internal helpers

#' @noRd
.rspde_cv_Efnorm <- function(mu, sigma) {
  sigma * sqrt(2 / pi) * exp(-(mu^2) / (2 * sigma^2)) +
    mu * (1 - 2 * stats::pnorm(-mu / sigma))
}

#' @noRd
.rspde_cv_Exx <- function(mu, sigma) {
  .rspde_cv_Efnorm(0, sqrt(2) * sigma)
}

#' @noRd
.rspde_cv_Exy <- function(mu, sigma, y) {
  .rspde_cv_Efnorm(mu - y, sigma)
}

#' @noRd
.rspde_cv_CRPS <- function(y, mu, sigma) {
  -.rspde_cv_Exy(mu, sigma, y) + 0.5 * .rspde_cv_Exx(mu, sigma)
}

#' @noRd
.rspde_cv_SCRPS <- function(y, mu, sigma) {
  -.rspde_cv_Exy(mu, sigma, y) / .rspde_cv_Exx(mu, sigma) -
    0.5 * log(.rspde_cv_Exx(mu, sigma))
}

#' @noRd
.rspde_cv_graph_object <- function(object) {
  if (!is.null(object$graph)) {
    return(object$graph)
  }
  if (!is.null(object$latent_model) && !is.null(object$latent_model$graph)) {
    return(object$latent_model$graph)
  }
  NULL
}

#' Build the data structures needed to drive cross-validation predictions.
#'
#' Returns a list with:
#'   - n_obs, repl_vec, y_vec
#'   - newdata_fn(idx): returns a list/data.frame slice of the data
#'   - loc_fn(idx): returns a loc matrix for the indices (NULL when newdata
#'       already encodes locations as for graph fits)
#'   - has_graph, is_spacetime
#'   - predict_extra: list of extra args to forward to predict.rspde_lme
#'   - response_var, formula, model
#'   - graph (R6 metric_graph) for graph fits
#'
#' @noRd
.rspde_cv_build_data <- function(object, user_data) {

  n_obs <- as.integer(object$nobs)
  if (is.na(n_obs) || n_obs <= 0L) {
    stop("Cannot determine number of observations from the fitted object.")
  }

  repl_vec <- object$repl
  if (length(repl_vec) != n_obs) {
    stop("Length of object$repl does not match object$nobs.")
  }

  response_name <- as.character(object$response_var)
  is_spacetime  <- isTRUE(object$spacetime)
  has_graph     <- isTRUE(object$has_graph)
  is_null       <- isTRUE(object$null_model)

  # null (pure OLS) fit 
  if (is_null) {
    mm <- as.matrix(object$model_matrix)
    if (ncol(mm) < 1L) {
      stop("Null-model fit has no model matrix to perform CV on.")
    }
    y_vec <- as.numeric(mm[, 1])
    if (ncol(mm) >= 2L) {
      X_design <- mm[, -1L, drop = FALSE]
    } else {
      X_design <- matrix(1, nrow = n_obs, ncol = 1L,
                         dimnames = list(NULL, "(Intercept)"))
    }
    beta_full <- as.numeric(object$coeff$fixed_effects)
    if (length(beta_full) != ncol(X_design)) {
      stop("Length of fixed-effect coefficients (", length(beta_full),
           ") does not match number of design columns (", ncol(X_design), ").")
    }
    sigma_e_full <- as.numeric(object$coeff$measurement_error[[1]])

    return(list(n_obs = n_obs, repl_vec = repl_vec, y_vec = y_vec,
                X_design = X_design,
                beta_full = beta_full, sigma_e_full = sigma_e_full,
                response_var = response_name,
                is_null = TRUE,
                has_graph = FALSE, is_spacetime = FALSE))
  }

  # graph fits 
  if (has_graph) {
    graph <- .rspde_cv_graph_object(object)
    if (is.null(graph)) {
      stop("Graph object not found on the fit (object$graph / object$latent_model$graph).")
    }
    graph_data <- graph$.__enclos_env__$private$data
    if (is.null(graph_data)) {
      stop("The metric graph does not store any private data.")
    }

    keep <- rep(TRUE, length(graph_data[[".edge_number"]]))
    if (response_name %in% names(graph_data)) {
      keep <- keep & !is.na(graph_data[[response_name]])
    }
    if (!is.null(object$which_repl) && ".group" %in% names(graph_data)) {
      keep <- keep & (as.character(graph_data[[".group"]]) %in%
                        as.character(object$which_repl))
    }
    full_data <- lapply(graph_data, function(col) col[keep])

    n_full <- length(full_data[[".edge_number"]])
    if (n_full != n_obs) {
      stop("Filtered graph data length (", n_full,
           ") does not match object$nobs (", n_obs, ").")
    }

    y_vec <- as.numeric(full_data[[response_name]])
    if (length(y_vec) != n_obs) {
      y_vec <- as.numeric(object$model_matrix[, 1])
    }

    newdata_fn <- function(idx) {
      lapply(full_data, function(col) col[idx])
    }
    loc_fn <- function(idx) NULL
    predict_extra <- list(edge_number = ".edge_number",
                          distance_on_edge = ".distance_on_edge",
                          normalized = TRUE)

    return(list(n_obs = n_obs, repl_vec = repl_vec, y_vec = y_vec,
                full_data = full_data,
                newdata_fn = newdata_fn, loc_fn = loc_fn,
                has_graph = TRUE, is_spacetime = is_spacetime,
                predict_extra = predict_extra,
                response_var = response_name,
                graph = graph,
                formula = object$formula,
                model = object$latent_model,
                use_data_from_graph = TRUE))
  }

  # non-graph fits 
  loc_mat <- as.matrix(object$loc)
  if (nrow(loc_mat) != n_obs) {
    stop("Number of rows of object$loc does not match object$nobs.")
  }

  if (!is.null(user_data)) {
    user_data <- as.data.frame(user_data)
    if (nrow(user_data) != n_obs) {
      stop("Argument 'data' has ", nrow(user_data), " rows but the fit has ",
           n_obs, " observations (after filtering). Please supply data ",
           "aligned with the rows used in fitting.")
    }
    full_data <- as.list(user_data)
  } else {
    mm <- as.data.frame(object$model_matrix)
    if (ncol(mm) >= 1) {
      names(mm)[1] <- response_name
    }
    if ("(Intercept)" %in% names(mm)) {
      mm[["(Intercept)"]] <- NULL
    }
    full_data <- as.list(mm)
  }

  if (!response_name %in% names(full_data)) {
    full_data[[response_name]] <- as.numeric(object$model_matrix[, 1])
  }
  y_vec <- as.numeric(full_data[[response_name]])

  time_vec <- NULL
  if (is_spacetime) {
    time_vec <- as.numeric(object$time)
    if (length(time_vec) != n_obs) {
      stop("Length of object$time does not match object$nobs for spacetime fit.")
    }
  }

  newdata_fn <- function(idx) {
    out <- lapply(full_data, function(col) col[idx])
    out
  }
  loc_fn <- function(idx) loc_mat[idx, , drop = FALSE]
  predict_extra <- list()
  if (is_spacetime) {
    predict_extra$time <- NULL  # we pass a numeric vector directly below
  }

  list(n_obs = n_obs, repl_vec = repl_vec, y_vec = y_vec,
       full_data = full_data,
       newdata_fn = newdata_fn, loc_fn = loc_fn,
       has_graph = FALSE, is_spacetime = is_spacetime,
       predict_extra = predict_extra,
       response_var = response_name,
       loc_mat = loc_mat, time_vec = time_vec,
       formula = object$formula,
       model = object$latent_model,
       use_data_from_graph = FALSE)
}

#' Build the argument list for a predict.rspde_lme call.
#' @noRd
.rspde_cv_predict_args <- function(cv_info, idx, compute_variances,
                                   which_repl) {
  args <- list(newdata = cv_info$newdata_fn(idx),
               compute_variances = compute_variances,
               which_repl = which_repl)
  if (!cv_info$has_graph) {
    args$loc <- cv_info$loc_fn(idx)
    if (cv_info$is_spacetime && !is.null(cv_info$time_vec)) {
      args$time <- cv_info$time_vec[idx]
    }
  }
  for (nm in names(cv_info$predict_extra)) {
    args[[nm]] <- cv_info$predict_extra[[nm]]
  }
  args
}

#' Cross-validation step for a pure OLS (null-model) fit.
#'
#' For pseudo-CV (`true_CV = FALSE`) the full-fit coefficients are reused;
#' for `true_CV = TRUE` the OLS coefficients are refit on `train_indices`.
#' The predictive variance is the standard new-observation variance,
#' `sigma_e^2 * (1 + x' (X'X)^{-1} x)`, computed from whichever design
#' (full or training) was used to estimate beta.
#'
#' @noRd
.rspde_cv_ols_fold <- function(cv_info, train_indices, test_indices, true_CV) {
  X_test <- cv_info$X_design[test_indices, , drop = FALSE]

  if (isTRUE(true_CV)) {
    X_train <- cv_info$X_design[train_indices, , drop = FALSE]
    y_train <- cv_info$y_vec[train_indices]
    fit <- stats::lm.fit(X_train, y_train)
    beta <- fit$coefficients
    n_train <- length(y_train)
    rdf <- fit$df.residual
    if (is.na(rdf) || rdf <= 0L) {
      sigma_e <- cv_info$sigma_e_full
    } else {
      sigma_e <- sqrt(sum(fit$residuals^2) / rdf)
    }
    X_used <- X_train
  } else {
    beta    <- cv_info$beta_full
    sigma_e <- cv_info$sigma_e_full
    X_used  <- cv_info$X_design
  }

  na_beta <- is.na(beta)
  if (any(na_beta)) {
    beta[na_beta] <- 0
  }
  mu <- as.numeric(X_test %*% beta)

  XtX_inv <- tryCatch(chol2inv(chol(crossprod(X_used))),
                      error = function(e) NULL)
  if (is.null(XtX_inv)) {
    XtX_inv <- tryCatch(solve(crossprod(X_used)), error = function(e) NULL)
  }
  if (is.null(XtX_inv)) {
    var_pred <- rep(sigma_e^2, length(test_indices))
  } else {
    leverage <- rowSums((X_test %*% XtX_inv) * X_test)
    var_pred <- sigma_e^2 * (1 + as.numeric(leverage))
  }

  list(mu = mu, var = var_pred)
}


#' Refit an rspde_lme model on a training fold.
#' @noRd
.rspde_cv_refit <- function(object, cv_info, test_indices,
                            parallel_fitting = FALSE,
                            n_cores = 1L,
                            nelder_mead_init = FALSE) {

  refit_args <- list(
    formula = object$formula,
    model = object$latent_model,
    optim_method = object$optim_method,
    optim_controls = if (!is.null(object$optim_controls)) object$optim_controls else list(),
    previous_fit = object,
    parallel = parallel_fitting,
    n_cores = n_cores,
    nelder_mead_init = nelder_mead_init
  )

  if (cv_info$has_graph) {
    graph <- cv_info$graph$clone()
    gd <- graph$.__enclos_env__$private$data
    response_name <- cv_info$response_var
    if (!is.null(gd) && response_name %in% names(gd)) {
      keep <- rep(TRUE, length(gd[[".edge_number"]]))
      if (!is.null(object$which_repl) && ".group" %in% names(gd)) {
        keep <- keep & (as.character(gd[[".group"]]) %in%
                          as.character(object$which_repl))
      }
      keep_idx_in_filtered <- which(keep)
      y_col <- gd[[response_name]]
      y_col[keep_idx_in_filtered[test_indices]] <- NA
      gd[[response_name]] <- y_col
      graph$.__enclos_env__$private$data <- gd
    }
    refit_args$loc <- c(".edge_number", ".distance_on_edge")
    refit_args$model$graph <- graph
    refit_args$use_data_from_graph <- TRUE
    refit_args$repl <- ".group"
    refit_args$which_repl <- object$which_repl
  } else {
    df <- as.data.frame(cv_info$full_data, stringsAsFactors = FALSE)
    df[[cv_info$response_var]][test_indices] <- NA
    refit_args$data <- df
    refit_args$loc <- cv_info$loc_mat
    refit_args$use_data_from_graph <- FALSE
    refit_args$repl <- object$repl
    refit_args$which_repl <- object$which_repl
    if (cv_info$is_spacetime && !is.null(cv_info$time_vec)) {
      refit_args$loc_time <- cv_info$time_vec
    }
  }

  if (!is.null(object$rspde_order)) {
    refit_args$rspde_order <- object$rspde_order
  }
  if (isTRUE(object$mean_correction)) {
    refit_args$mean_correction <- TRUE
  }

  do.call(rspde_lme, refit_args)
}
