#'
#' @title rSPDE inlabru mapper
#' @name bru_get_mapper.inla_rspde
#' @param model An `inla_rspde` object for which to construct or extract a mapper
#' @param \dots Arguments passed on to other methods
#' @rdname bru_get_mapper.inla_rspde
#' @rawNamespace if (getRversion() >= "3.6.0") {
#'   S3method(inlabru::bru_get_mapper, inla_rspde)
#'   S3method(inlabru::ibm_n, bru_mapper_inla_rspde)
#'   S3method(inlabru::ibm_values, bru_mapper_inla_rspde)
#'   S3method(inlabru::ibm_jacobian, bru_mapper_inla_rspde)
#' }
#'
#' @examples
#' \donttest{ #tryCatch version
          #' tryCatch({
#' if (requireNamespace("INLA", quietly = TRUE) &&
#'   requireNamespace("inlabru", quietly = TRUE)) {
#'   library(INLA)
#'   library(inlabru)
#'
#'   set.seed(123)
#'   m <- 100
#'   loc_2d_mesh <- matrix(runif(m * 2), m, 2)
#'   mesh_2d <- inla.mesh.2d(
#'     loc = loc_2d_mesh,
#'     cutoff = 0.05,
#'     max.edge = c(0.1, 0.5)
#'   )
#'   sigma <- 1
#'   range <- 0.2
#'   nu <- 0.8
#'   kappa <- sqrt(8 * nu) / range
#'   op <- matern.operators(
#'     mesh = mesh_2d, nu = nu,
#'     range = range, sigma = sigma, m = 2,
#'     parameterization = "matern"
#'   )
#'   u <- simulate(op)
#'   A <- inla.spde.make.A(
#'     mesh = mesh_2d,
#'     loc = loc_2d_mesh
#'   )
#'   sigma.e <- 0.1
#'   y <- A %*% u + rnorm(m) * sigma.e
#'   y <- as.vector(y)
#'
#'   data_df <- data.frame(
#'     y = y, x1 = loc_2d_mesh[, 1],
#'     x2 = loc_2d_mesh[, 2]
#'   )
#'   rspde_model <- rspde.matern(
#'     mesh = mesh_2d,
#'     nu_upper_bound = 2
#'   )
#'
#'   cmp <- y ~ Intercept(1) +
#'     field(cbind(x1,x2), model = rspde_model)
#'
#'
#'   rspde_fit <- bru(cmp, data = data_df)
#'   summary(rspde_fit)
#' }
#' #stable.tryCatch
          #' }, error = function(e){print("Could not run the example")})
#' }
bru_get_mapper.inla_rspde <- function(model, ...) {
  stopifnot(requireNamespace("inlabru"))
  inlabru_version <- as.character(packageVersion("inlabru"))
  if(inlabru_version >= "2.11.1.9022"){
      n_rep <- model[["rspde.order"]] + 1
      if((model[["est_nu"]] == 0L) && (model[["integer.nu"]])){
          n_rep <- 1
      }
    inlabru::bru_mapper_repeat(inlabru::bru_mapper(model[["mesh"]]), n_rep = n_rep)
  } else{
    mapper <- list(model = model)
    inlabru::bru_mapper_define(mapper, new_class = "bru_mapper_inla_rspde")
  }
}

#' @param mapper A `bru_mapper_inla_rspde` object
#' @rdname bru_get_mapper.inla_rspde
ibm_n.bru_mapper_inla_rspde <- function(mapper, ...) {
  model <- mapper[["model"]]
  integer_nu <- model$integer.nu
  rspde_order <- model$rspde.order
  if (integer_nu) {
    factor_rspde <- 1
  } else {
    factor_rspde <- rspde_order + 1
  }
  factor_rspde * model$n.spde
}
#' @rdname bru_get_mapper.inla_rspde
ibm_values.bru_mapper_inla_rspde <- function(mapper, ...) {
  seq_len(inlabru::ibm_n(mapper))
}
#' @param input The values for which to produce a mapping matrix
#' @rdname bru_get_mapper.inla_rspde
ibm_jacobian.bru_mapper_inla_rspde <- function(mapper, input, ...) {
  model <- mapper[["model"]]
  if (!is.null(input) && !is.matrix(input) && !inherits(input, "Spatial")) {
    input <- as.matrix(input)
  }

  if (model$est_nu) {
    nu <- NULL
  } else {
    nu <- model$nu
  }

  rspde_order <- model$rspde.order
  rSPDE::rspde.make.A(
    mesh = model$mesh, loc = input,
    rspde.order = rspde_order,
    nu = nu
  )
}

#' @noRd

process_formula <- function(bru_result) {
  form <- bru_result$bru_info$model$formula[3]
  form <- as.character(form)
  form <- strsplit(form, "f\\(")
  form <- form[[1]]
  form <- form[-1]
  form_proc <- sub(",.*", "", strsplit(form, "f\\(")[1])
  if (length(form) > 1) {
    for (i in 2:(length(form))) {
      form_proc <- paste(form_proc, " + ", sub(",.*", "", strsplit(form, "f\\(")[i]))
    }
  }
  form_proc <- paste("~", "linkfuninv(", form_proc, ")")
  return(stats::as.formula(form_proc))
}

#' @noRd

process_formula_lhoods <- function(bru_result, like_number) {
      form <- bru_result$bru_info$lhoods[[like_number]]$formula[3]
      form <- as.character(form)
      if(form == "."){
        return(process_formula(bru_result))
      }
      form_proc <- paste("~", "linkfuninv(", form, ")")
      return(stats::as.formula(form_proc))
    }

#' @noRd
# Function to process the link function

process_link <- function(link_name) {
  return_link <- switch(link_name,
    "log" = function(x) {
      INLA::inla.link.log(x, inverse = TRUE)
    },
    "invlog" = function(x) {
      INLA::inla.link.invlog(x, inverse = TRUE)
    },
    "logit" = function(x) {
      INLA::inla.link.logit(x, inverse = TRUE)
    },
    "invlogit" = function(x) {
      INLA::inla.link.invlogit(x, inverse = TRUE)
    },
    "probit" = function(x) {
      INLA::inla.link.probit(x, inverse = TRUE)
    },
    "invprobit" = function(x) {
      INLA::inla.link.invprobit(x, inverse = TRUE)
    },
    "cloglog" = function(x) {
      INLA::inla.link.cloglog(x, inverse = TRUE)
    },
    "invcloglog" = function(x) {
      INLA::inla.link.invcloglog(x, inverse = TRUE)
    },
    "tan" = function(x) {
      INLA::inla.link.tan(x, inverse = TRUE)
    },
    "invtan" = function(x) {
      INLA::inla.link.invtan(x, inverse = TRUE)
    },
    "identity" = function(x) {
      INLA::inla.link.identity(x, inverse = TRUE)
    },
    "invidentity" = function(x) {
      INLA::inla.link.invidentity(x, inverse = TRUE)
    }
  )
  return(return_link)
}

#' @noRd

bru_rerun_with_data <- function(result, idx_data, true_CV, fit_verbose) {
  stopifnot(inherits(result, "bru"))
  if (!true_CV) {
    options <- list(control.mode = list(
      theta = result$mode$theta,
      fixed=TRUE
    ))
  } else {
    options <- list()
  }

  if (fit_verbose) {
    options$verbose <- TRUE
  } else {
    options$verbose <- FALSE
  }

  info <- result[["bru_info"]]
  info[["options"]] <- inlabru::bru_call_options(
    inlabru::bru_options(
      info[["options"]],
      inlabru::as.bru_options(options)
    )
  )

  original_timings <- result[["bru_timings"]]

  for(i_like in seq_along(info[["lhoods"]])){
    info[["lhoods"]][[i_like]]$response_data$BRU_response[-idx_data[[i_like]]] <- NA
  }

  result <- inlabru::iinla(
      model = info[["model"]],
      lhoods = info[["lhoods"]],
      initial = result,
      options = info[["options"]]
    )

  new_timings <- result[["bru_iinla"]][["timings"]]$Iteration >
    max(original_timings$Iteration)
  result$bru_timings <-
    rbind(
      original_timings,
      result[["bru_iinla"]][["timings"]][new_timings, , drop = FALSE]
    )

  # Add bru information to the result
  result$bru_info <- info
  class(result) <- c("bru", class(result))
  return(result)
}


#' @noRd

get_post_var <- function(density_df) {
  min_x <- min(density_df[, "x"])
  max_x <- max(density_df[, "x"])
  denstemp <- function(x) {
    dens <- sapply(x, function(z) {
      if (z < min_x) {
        return(0)
      } else if (z > max_x) {
        return(0)
      } else {
        return(approx(x = density_df[, "x"], y = density_df[, "y"], xout = z)$y)
      }
    })
    return(dens)
  }

  post_var <- stats::integrate(
    f = function(z) {
      denstemp(z) * 1 / z
    }, lower = min_x, upper = max_x,
    subdivisions = nrow(density_df),
    stop.on.error = FALSE
  )$value

  return(post_var)
}

#' @name cross_validation
#' @title Perform cross-validation on a list of fitted models.
#' @description Obtain several scores for a list of fitted models according
#' to a folding scheme.
#' @param models A fitted model obtained from calling the `bru()` function or a list of fitted models. 
#'        All models in the list must have the same number of likelihoods and must be fitted to 
#'        identical datasets. 
#' @param model_names A vector containing the names of the models to appear in the returned `data.frame`. If `NULL`, the names will be of the form `Model 1`, `Model 2`, and so on. By default, it will try to obtain the name from the models list.
#' @param scores A vector containing the scores to be computed. The options are "mse", "crps", "scrps", "dss", "wcrps" and "swcrps". By default, all scores are computed.
#' @param cv_type The type of the folding to be carried out. The options are `k-fold` for `k`-fold cross-validation, in which case the parameter `k` should be provided,
#' `loo`, for leave-one-out and `lpo` for leave-percentage-out, in this case, the parameter `percentage` should be given, and also the `number_folds`
#' with the number of folds to be done. The default is `k-fold`.
#' @param k The number of folds to be used in `k`-fold cross-validation. Will only be used if `cv_type` is `k-fold`.
#' @param percentage The percentage (from 1 to 99) of the data to be used to train the model. Will only be used if `cv_type` is `lpo`.
#' @param number_folds Number of folds to be done if `cv_type` is `lpo`.
#' @param weight_thr When computing "wcrps" or "swcrps", the threshold to be used to compute the weights. Must be supplied if any of these scores are requested. No default value is provided.
#' @param n_samples Number of samples to compute the posterior statistics to be used to compute the scores.
#' @param return_scores_folds If `TRUE`, the scores for each fold will also be returned.
#' @param orientation_results character vector. The options are "negative" and "positive". If "negative", the smaller the scores the better. If "positive", the larger the scores the better.
#' @param include_best Should a row indicating which model was the best for each score be included?
#' @param train_test_indexes A list where each element corresponds to a fold. Each fold contains:
#' - `train`: A list of training index vectors, one for each likelihood.
#' - `test`: A list of test index vectors, one for each likelihood, with the same length as `train`.
#' This list is typically obtained by setting the argument `return_train_test` to `TRUE`. 
#' When supplying `train_test_indexes`, the `cv_type`, `k`, `percentage` and `number_folds` arguments are ignored.
#' @param return_train_test Logical. Should the training and test indexes be returned? If 'TRUE' the train and test indexes will the 'train_test' element of the returned list. 
#' @param return_post_samples If `TRUE` the posterior samples will be included in the returned list.
#' @param return_true_test_values If `TRUE` the true test values will be included in the returned list.
#' @param parallelize_RP Logical. Should the computation of CRPS and SCRPS (and for some cases, DSS) be parallelized?
#' @param n_cores_RP Number of cores to be used if `parallelize_rp` is `TRUE`.
#' @param true_CV Should a `TRUE` cross-validation be performed? If `TRUE` the models will be fitted on the training dataset. If `FALSE`, the parameters will be kept fixed at the ones obtained in the result object.
#' @param save_settings Logical. If `TRUE`, the settings used in the cross-validation will also be returned.
#' @param print Should partial results be printed throughout the computation?
#' @param fit_verbose Should INLA's run during cross-validation be verbose?
#' @return A data.frame with the fitted models and the corresponding scores.
#' @export
cross_validation <- function(models, model_names = NULL, scores = c("mae", "mse", "crps", "scrps", "dss"),
                             cv_type = c("k-fold", "loo", "lpo"),
                             weight_thr=NULL,
                             k = 5, percentage = 20, number_folds = 10,
                             n_samples = 1000, return_scores_folds = FALSE,
                             orientation_results = c("negative", "positive"),
                             include_best = TRUE,
                             train_test_indexes = NULL,
                             return_train_test = FALSE,
                             return_post_samples = FALSE,
                             return_true_test_values = FALSE,
                             parallelize_RP = FALSE, n_cores_RP = parallel::detectCores() - 1,
                             true_CV = TRUE, save_settings = FALSE,
                             print = TRUE,
                             fit_verbose = FALSE) {
  orientation_results <- orientation_results[[1]]
  if (!(orientation_results %in% c("positive", "negative"))) {
    stop("orientation_results must be either 'positive' or 'negative'!")
  }

  if(any(scores %in% c("wcrps", "swcrps")) && is.null(weight_thr)){
    stop("weight_thr must be supplied if 'wcrps' or 'swcrps' are requested!")
  }

  scores <- intersect(scores, c("mae", "mse", "crps", "scrps", "dss", "wcrps", "swcrps"))

  cv_type <- cv_type[[1]]
  if (!(cv_type %in% c("k-fold", "loo", "lpo"))) {
    stop("The possible options for cv_type are 'k-fold', 'loo' or 'lpo'!")
  }

  if (!is.numeric(percentage)) {
    stop("percentage must be numeric!")
  }

  if (percentage %% 1 != 0) {
    warning("Non-integer percentage given, it will be rounded to an integer number.")
    percentage <- round(percentage)
  }

  if (percentage <= 0 || percentage >= 100) {
    stop("percentage must be a number between 1 and 99!")
  }

  if (!is.numeric(number_folds)) {
    stop("number_folds must be numeric!")
  }

  if (number_folds %% 1 != 0) {
    warning("Non-integer number_folds given, it will be rounded to an integer number.")
    number_folds <- round(number_folds)
  }

  if (number_folds <= 0) {
    stop("number_folds must be positive!")
  }

  if (inherits(models, "bru")) {
    models <- list(models)
  } else {
    for (i in 1:length(models)) {
      if (!inherits(models[[i]], "bru")) {
        stop("models must be either a result from a bru call or a list of results from bru() calls!")
      }
    }
  }

  # The number of likelihoods. All models must have the same number of likelihoods.

  n_likelihoods <- length(models[[1]]$bru_info$lhoods)

  for (model_number in seq_along(models)) {
      if (length(models[[model_number]]$bru_info$lhoods) != n_likelihoods) {
          stop(paste("Model", model_number, "does not have the same number of likelihoods as the first model."))
      }
  }

  if (is.null(model_names) && is.list(models)) {
    model_names <- names(models)
  }

  if (!is.null(model_names)) {
    if (!is.character(model_names)) {
      stop("model_names must be a vector of strings!")
    }
    if (length(models) != length(model_names)) {
      stop("model_names must contain one name for each model!")
    }
  } else {
    model_names <- vector(mode = "character", length(models))
    for (i in 1:length(models)) {
      model_names[i] <- paste("Model", i)
    }
  }

  if (!is.numeric(n_samples)) {
    stop("n_samples must be numeric!")
  }

  if (n_samples %% 1 != 0) {
    warning("Non-integer n_samples given, it will be rounded to an integer number.")
    n_samples <- round(n_samples)
  }

  if (n_samples <= 0) {
    stop("n_samples must be positive!")
  }

  if (parallelize_RP) {
    cluster_tmp <- parallel::makeCluster(n_cores_RP)
    doParallel::registerDoParallel(cluster_tmp)
  }

  # Creating lists of train and test datasets

  if (is.null(train_test_indexes)) {
    # Observe that here we are assuming that all models use the same data, which we added in the description as an assumption.

    data_list <- lapply(seq_len(n_likelihoods), function(i) {
        models[[1]]$bru_info$lhoods[[i]]$data
    })
    train_test_indexes <- create_train_test_indices(data_list,
      cv_type = cv_type,
      k = k, percentage = percentage, number_folds = number_folds
    )
  } else {
    if (!is.list(train_test_indexes)) {
        stop("train_test_indexes should be a list!")
    }

    for (i in seq_along(train_test_indexes)) {
        fold_entry <- train_test_indexes[[i]]

        if (!is.list(fold_entry)) {
            stop(paste("train_test_indexes[[", i, "]] should be a list!", sep = ""))
        }

        if (is.null(fold_entry[["train"]])) {
            stop(paste("train_test_indexes[[", i, "]] must contain a 'train' element.", sep = ""))
        }

        if (is.null(fold_entry[["test"]])) {
            stop(paste("train_test_indexes[[", i, "]] must contain a 'test' element.", sep = ""))
        }

        if (!is.list(fold_entry[["train"]])) {
            stop(paste("train_test_indexes[[", i, "]][['train']] must be a list!", sep = ""))
        }

        if (!is.list(fold_entry[["test"]])) {
            stop(paste("train_test_indexes[[", i, "]][['test']] must be a list!", sep = ""))
        }

        if (length(fold_entry[["train"]]) != n_likelihoods) {
            stop(paste("train_test_indexes[[", i, "]][['train']] must have length equal to n_likelihoods:", n_likelihoods))
        }

        if (length(fold_entry[["test"]]) != n_likelihoods) {
            stop(paste("train_test_indexes[[", i, "]][['test']] must have length equal to n_likelihoods:", n_likelihoods))
        }
    }
  }

  post_samples <- list()

  true_test_values = list()

  for (model_number in 1:length(models)) {
    n_likelihoods <- length(models[[model_number]]$bru_info$lhoods)
    post_samples[[model_names[[model_number]]]] <- vector(mode = "list", length = length(train_test_indexes))
    if(return_true_test_values){
      true_test_values[[model_names[[model_number]]]] <- vector(mode = "list", length = length(train_test_indexes))
    }    
    for(j in seq_along(train_test_indexes)){
      post_samples[[model_names[[model_number]]]][[j]] <- vector(mode = "list", length = n_likelihoods)
      if(return_true_test_values){
        true_test_values[[model_names[[model_number]]]][[j]] <- vector(mode = "list", length = n_likelihoods)
      }    
    }
  }
  # Perform the cross-validation

  result_df <- data.frame(Model = model_names)

  n_folds <- length(train_test_indexes)
  n_models <- length(models)
  
  dss <- lapply(1:n_likelihoods, function(i) matrix(numeric(n_folds * n_models), ncol = n_models))
  mse <- lapply(1:n_likelihoods, function(i) matrix(numeric(n_folds * n_models), ncol = n_models))
  mae <- lapply(1:n_likelihoods, function(i) matrix(numeric(n_folds * n_models), ncol = n_models))
  crps <- lapply(1:n_likelihoods, function(i) matrix(numeric(n_folds * n_models), ncol = n_models))
  scrps <- lapply(1:n_likelihoods, function(i) matrix(numeric(n_folds * n_models), ncol = n_models))
  wcrps <- lapply(1:n_likelihoods, function(i) matrix(numeric(n_folds * n_models), ncol = n_models))
  swcrps <- lapply(1:n_likelihoods, function(i) matrix(numeric(n_folds * n_models), ncol = n_models))

  if(any(c("crps", "scrps", "dss", "wcrps", "swcrps") %in% scores)){
    new_n_samples <- 2 * n_samples
  } else {
    new_n_samples <- n_samples
  }

  for (fold in 1:length(train_test_indexes)) {
    for (model_number in 1:length(models)) {
      if (print) {
        cat(paste("Fold:", fold, "/", length(train_test_indexes), "\n"))
        if (!is.null(model_names)) {
          cat(paste("Model:", model_names[[model_number]], "\n"))
        } else {
          cat(paste("Model:", model_number, "\n"))
        }
      }

      train_list <- train_test_indexes[[fold]][["train"]]
      test_list <- train_test_indexes[[fold]][["test"]]

      new_model <- bru_rerun_with_data(models[[model_number]], train_list, true_CV = true_CV, fit_verbose = fit_verbose)

      for(i_lik in 1:n_likelihoods){

          model_family <- models[[model_number]]$.args$family[[i_lik]]

          post_linear_predictors <- sample_posterior_linear_predictor(new_model, i_lik, test_list, new_n_samples, print)

          post_samples[[model_names[[model_number]]]][[fold]][[i_lik]] <- get_posterior_samples(
                      post_linear_predictors = post_linear_predictors, new_model = new_model, 
                      i_lik = i_lik, new_n_samples = new_n_samples, 
                      full_model = models[[model_number]], 
                      true_CV = true_CV, print = print)

          post_samples[[model_names[[model_number]]]][[fold]][[i_lik]] <-  do.call(rbind, post_samples[[model_names[[model_number]]]][[fold]][[i_lik]])
          
          test_data <- models[[model_number]]$bru_info$lhoods[[i_lik]]$response_data$BRU_response[test_list[[i_lik]]]
          
          if(return_true_test_values){
            true_test_values[[model_names[[model_number]]]][[fold]] <- test_data
          }


          if (!(model_family %in% c("stochvol", "stochvolln", "stochvolnig", "stochvolt"))) {
            posterior_mean <- rowMeans(post_samples[[model_names[[model_number]]]][[fold]][[i_lik]])
            if ("mse" %in% scores) {
              mse[[i_lik]][fold, model_number] <- mean((test_data - posterior_mean)^2, na.rm = TRUE)
              if (orientation_results == "positive") {
                mse[[i_lik]][fold, model_number] <- -mse[[i_lik]][fold, model_number]
              }
              if (print) {
                cat(paste0("MSE - Likelihood ",i_lik,": ", mse[[i_lik]][fold, model_number], "\n"))
              }
            }
            if ("mae" %in% scores) {
              mae[[i_lik]][fold, model_number] <- mean(abs(test_data - posterior_mean), na.rm = TRUE)
              if (orientation_results == "positive") {
                mae[[i_lik]][fold, model_number] <- -mae[[i_lik]][fold, model_number]
              }
              if (print) {
                cat(paste0("MAE - Likelihood ",i_lik,": ", mae[[i_lik]][fold, model_number], "\n"))
              }
            }            
          }
        
          if ("dss" %in% scores) {
            post_var <- rowMeans(post_samples[[model_names[[model_number]]]][[fold]][[i_lik]][, 1:n_samples, drop=FALSE]^2) - (rowMeans(post_samples[[model_names[[model_number]]]][[fold]][[i_lik]][, 1:n_samples, drop=FALSE]))^2

            dss[[i_lik]][fold, model_number] <- mean((test_data - rowMeans(post_samples[[model_names[[model_number]]]][[fold]][[i_lik]][, (n_samples + 1):(2 * n_samples), drop=FALSE]))^2 / post_var + log(post_var))
            if (print) {
              cat(paste("DSS - Likelihood ",i_lik,": ", dss[[i_lik]][fold, model_number], "\n"))
            }            
          }

        if(any(c("crps", "scrps", "wcrps", "swcrps") %in% scores)){
          Y1_sample <- post_samples[[model_names[[model_number]]]][[fold]][[i_lik]][, 1:n_samples, drop=FALSE]
          Y2_sample <- post_samples[[model_names[[model_number]]]][[fold]][[i_lik]][, (n_samples + 1):(2 * n_samples), drop=FALSE]
        }

        if(any(c("crps", "scrps") %in% scores)){
          if (parallelize_RP) {
            E1_tmp <- foreach::`%dopar%`(foreach::foreach(i = 1:length(test_data)), {
              mean(abs(Y1_sample[i,] - test_data[i]))
            })
            E2_tmp <- foreach::`%dopar%`(foreach::foreach(i = 1:length(test_data)), {
              mean(abs(Y1_sample[i,] - Y2_sample[i,]))
            })
          } else {
            E1_tmp <- lapply(1:length(test_data), function(i) {
              mean(abs(Y1_sample[i,] - test_data[i]))
            })
            E2_tmp <- lapply(1:length(test_data), function(i) {
              mean(abs(Y1_sample[i,] - Y2_sample[i,]))
            })
          }          
        }

        if(any(c("wcrps", "swcrps") %in% scores)){
          if (parallelize_RP) {
            E1_tmp_thr <- foreach::`%dopar%`(foreach::foreach(i = 1:length(test_data)), {
              mean(abs((Y1_sample[i,]>weight_thr)*(Y1_sample[i,]-weight_thr)-(test_data[i]>weight_thr)*(test_data[i]-weight_thr)))
            })
            E2_tmp_thr <- foreach::`%dopar%`(foreach::foreach(i = 1:length(test_data)), {
              mean(abs((Y1_sample[i,]>weight_thr)*(Y1_sample[i,]-weight_thr)-(Y2_sample[i,]>weight_thr)*(Y2_sample[i,]-weight_thr)))
            })
          } else {
            E1_tmp_thr <- lapply(1:length(test_data), function(i) {
              mean(abs((Y1_sample[i,]>weight_thr)*(Y1_sample[i,]-weight_thr)-(test_data[i]>weight_thr)*(test_data[i]-weight_thr)))
            })
            E2_tmp_thr <- lapply(1:length(test_data), function(i) {
              mean(abs((Y1_sample[i,]>weight_thr)*(Y1_sample[i,]-weight_thr)-(Y2_sample[i,]>weight_thr)*(Y2_sample[i,]-weight_thr)))
            })
          }    
        }

        if ("crps" %in% scores) {
            crps_temp <- lapply(1:length(test_data), function(i) {
              return(-E1_tmp[[i]] + 0.5 * E2_tmp[[i]])
            })

            crps_temp <- unlist(crps_temp)
            crps[[i_lik]][fold, model_number] <- mean(crps_temp)
            if (orientation_results == "negative") {
              crps[[i_lik]][fold, model_number] <- -crps[[i_lik]][fold, model_number]
            }

            if (print) {
              cat(paste0("CRPS - Likelihood ",i_lik,": ", crps[[i_lik]][fold, model_number], "\n"))
            }
          }

          if ("scrps" %in% scores) {
            scrps_temp <- lapply(1:length(test_data), function(i) {
              return(-E1_tmp[[i]] / E2_tmp[[i]] - 0.5 * log(E2_tmp[[i]]))
            })
          scrps_temp <- unlist(scrps_temp)
          scrps[[i_lik]][fold, model_number] <- mean(scrps_temp)
          if (orientation_results == "negative") {
            scrps[[i_lik]][fold, model_number] <- -scrps[[i_lik]][fold, model_number]
          }

          if (print) {
            cat(paste("SCRPS: - Likelihood ",i_lik,": ", scrps[[i_lik]][fold, model_number], "\n"))
          }
        }     

        if("wcrps" %in% scores){
            wcrps_temp <- lapply(1:length(test_data), function(i){
                return(0.5*E2_tmp_thr[[i]]-E1_tmp_thr[[i]])
            })
            wcrps_temp <- unlist(wcrps_temp)
            wcrps[[i_lik]][fold, model_number] <- mean(wcrps_temp)  
            if(orientation_results == "negative"){
                wcrps[[i_lik]][fold, model_number] <- - wcrps[[i_lik]][fold, model_number]
            }    

            if (print) {
              cat(paste0("wCRPS - Likelihood ",i_lik,": ", wcrps[[i_lik]][fold, model_number], "\n"))
            }  
        }

        if("swcrps" %in% scores){
          if(any(E2_tmp_thr == 0)){
            warning("swCRPS cannot be computed. Please, decrease `weight_thr` or increase `n_sample`.")
            swcrps[[i_lik]][fold, model_number] <- NA
          } else{
            swcrps_temp <- lapply(1:length(test_data), function(i){
                return(-E1_tmp_thr[[i]]/E2_tmp_thr[[i]] - 0.5*log(E2_tmp_thr[[i]]))
            })
            swcrps_temp <- unlist(swcrps_temp)
            swcrps[[i_lik]][fold, model_number] <- mean(swcrps_temp)  
          }

            if(orientation_results == "negative"){
                swcrps[[i_lik]][fold, model_number] <- - swcrps[[i_lik]][fold, model_number]
            }   
            
            if (print) {
              cat(paste0("swCRPS - Likelihood ",i_lik,": ", swcrps[[i_lik]][fold, model_number], "\n"))
            }                                               
        }        
      }
     }
    }



  if (n_likelihoods > 1) {
   if ("mse" %in% scores) {
     # Individual likelihood means
     for (i in 1:n_likelihoods) {
       mse_mean <- colMeans(mse[[i]])
       result_df <- data.frame(result_df, mse = mse_mean)
       names(result_df)[names(result_df) == "mse"] <- paste0("mse_lik", i)
     }
     # Total MSE (weighted average across likelihoods)
     mse_weighted <- matrix(0, nrow = nrow(mse[[1]]), ncol = ncol(mse[[1]]))
     total_weights <- matrix(0, nrow = nrow(mse[[1]]), ncol = ncol(mse[[1]]))
     for (i in 1:n_likelihoods) {
       weights <- sapply(1:nrow(mse[[i]]), function(fold) length(train_test_indexes[[fold]][["test"]][[i]]))
       mse_weighted <- mse_weighted + weights * mse[[i]]
       total_weights <- total_weights + weights
     }
     mse_total <- colMeans(mse_weighted / total_weights)
     result_df <- data.frame(result_df, mse_total = mse_total)
   }

   if ("mae" %in% scores) {
     for (i in 1:n_likelihoods) {
       mae_mean <- colMeans(mae[[i]])
       result_df <- data.frame(result_df, mae = mae_mean)
       names(result_df)[names(result_df) == "mae"] <- paste0("mae_lik", i)
     }
     mae_weighted <- matrix(0, nrow = nrow(mae[[1]]), ncol = ncol(mae[[1]]))
     total_weights <- matrix(0, nrow = nrow(mae[[1]]), ncol = ncol(mae[[1]]))
     for (i in 1:n_likelihoods) {
       weights <- sapply(1:nrow(mae[[i]]), function(fold) length(train_test_indexes[[fold]][["test"]][[i]]))
       mae_weighted <- mae_weighted + weights * mae[[i]]
       total_weights <- total_weights + weights
     }
     mae_total <- colMeans(mae_weighted / total_weights)
     result_df <- data.frame(result_df, mae_total = mae_total)
   }

   if ("dss" %in% scores) {
     for (i in 1:n_likelihoods) {
       dss_mean <- colMeans(dss[[i]])
       result_df <- data.frame(result_df, dss = dss_mean)
       names(result_df)[names(result_df) == "dss"] <- paste0("dss_lik", i)
     }
     dss_weighted <- matrix(0, nrow = nrow(dss[[1]]), ncol = ncol(dss[[1]]))
     total_weights <- matrix(0, nrow = nrow(dss[[1]]), ncol = ncol(dss[[1]]))
     for (i in 1:n_likelihoods) {
       weights <- sapply(1:nrow(dss[[i]]), function(fold) length(train_test_indexes[[fold]][["test"]][[i]]))
       dss_weighted <- dss_weighted + weights * dss[[i]]
       total_weights <- total_weights + weights
     }
     dss_total <- colMeans(dss_weighted / total_weights)
     result_df <- data.frame(result_df, dss_total = dss_total)
   }

   if ("crps" %in% scores) {
     for (i in 1:n_likelihoods) {
       crps_mean <- colMeans(crps[[i]])
       result_df <- data.frame(result_df, crps = crps_mean)
       names(result_df)[names(result_df) == "crps"] <- paste0("crps_lik", i)
     }
     crps_weighted <- matrix(0, nrow = nrow(crps[[1]]), ncol = ncol(crps[[1]]))
     total_weights <- matrix(0, nrow = nrow(crps[[1]]), ncol = ncol(crps[[1]]))
     for (i in 1:n_likelihoods) {
       weights <- sapply(1:nrow(crps[[i]]), function(fold) length(train_test_indexes[[fold]][["test"]][[i]]))
       crps_weighted <- crps_weighted + weights * crps[[i]]
       total_weights <- total_weights + weights
     }
     crps_total <- colMeans(crps_weighted / total_weights)
     result_df <- data.frame(result_df, crps_total = crps_total)
   }

   if ("scrps" %in% scores) {
     for (i in 1:n_likelihoods) {
       scrps_mean <- colMeans(scrps[[i]])
       result_df <- data.frame(result_df, scrps = scrps_mean)
       names(result_df)[names(result_df) == "scrps"] <- paste0("scrps_lik", i)
     }
     scrps_weighted <- matrix(0, nrow = nrow(scrps[[1]]), ncol = ncol(scrps[[1]]))
     total_weights <- matrix(0, nrow = nrow(scrps[[1]]), ncol = ncol(scrps[[1]]))
     for (i in 1:n_likelihoods) {
       weights <- sapply(1:nrow(scrps[[i]]), function(fold) length(train_test_indexes[[fold]][["test"]][[i]]))
       scrps_weighted <- scrps_weighted + weights * scrps[[i]]
       total_weights <- total_weights + weights
     }
     scrps_total <- colMeans(scrps_weighted / total_weights)
     result_df <- data.frame(result_df, scrps_total = scrps_total)
   }

   if ("wcrps" %in% scores) {
     for (i in 1:n_likelihoods) {
       wcrps_mean <- colMeans(wcrps[[i]])
       result_df <- data.frame(result_df, wcrps = wcrps_mean)
       names(result_df)[names(result_df) == "wcrps"] <- paste0("wcrps_lik", i)
     }
     wcrps_weighted <- matrix(0, nrow = nrow(wcrps[[1]]), ncol = ncol(wcrps[[1]]))
     total_weights <- matrix(0, nrow = nrow(wcrps[[1]]), ncol = ncol(wcrps[[1]]))
     for (i in 1:n_likelihoods) {
       weights <- sapply(1:nrow(wcrps[[i]]), function(fold) length(train_test_indexes[[fold]][["test"]][[i]]))
       wcrps_weighted <- wcrps_weighted + weights * wcrps[[i]]
       total_weights <- total_weights + weights
     }
     wcrps_total <- colMeans(wcrps_weighted / total_weights)
     result_df <- data.frame(result_df, wcrps_total = wcrps_total)
   }

   if ("swcrps" %in% scores) {
     for (i in 1:n_likelihoods) {
       swcrps_mean <- colMeans(swcrps[[i]])
       result_df <- data.frame(result_df, swcrps = swcrps_mean)
       names(result_df)[names(result_df) == "swcrps"] <- paste0("swcrps_lik", i)
     }
     swcrps_weighted <- matrix(0, nrow = nrow(swcrps[[1]]), ncol = ncol(swcrps[[1]]))
     total_weights <- matrix(0, nrow = nrow(swcrps[[1]]), ncol = ncol(swcrps[[1]]))
     for (i in 1:n_likelihoods) {
       weights <- sapply(1:nrow(swcrps[[i]]), function(fold) length(train_test_indexes[[fold]][["test"]][[i]]))
       swcrps_weighted <- swcrps_weighted + weights * swcrps[[i]]
       total_weights <- total_weights + weights
     }
     swcrps_total <- colMeans(swcrps_weighted / total_weights)
     result_df <- data.frame(result_df, swcrps_total = swcrps_total)
   }

  } else {
   # Original code for n_likelihoods == 1
   if ("mse" %in% scores) {
     mse_mean <- colMeans(mse[[1]])
     result_df <- data.frame(result_df, mse = mse_mean)
   }

   if ("mae" %in% scores) {
     mae_mean <- colMeans(mae[[1]])
     result_df <- data.frame(result_df, mae = mae_mean)
   }

   if ("dss" %in% scores) {
     dss_mean <- colMeans(dss[[1]])
     result_df <- data.frame(result_df, dss = dss_mean)
   }

   if ("crps" %in% scores) {
     crps_mean <- colMeans(crps[[1]])
     result_df <- data.frame(result_df, crps = crps_mean)
   }

   if ("scrps" %in% scores) {
     scrps_mean <- colMeans(scrps[[1]])
     result_df <- data.frame(result_df, scrps = scrps_mean)
   }

   if ("wcrps" %in% scores) {
     wcrps_mean <- colMeans(wcrps[[1]])
     result_df <- data.frame(result_df, wcrps = wcrps_mean)
   }

   if ("swcrps" %in% scores) {
     swcrps_mean <- colMeans(swcrps[[1]])
     result_df <- data.frame(result_df, swcrps = swcrps_mean)
   }   
  }

  if (save_settings) {
    settings_list <- list(
      n_samples = n_samples, cv_type = cv_type, true_CV = true_CV,
      orientation_results = orientation_results
    )
    if (cv_type == "k-fold") {
      settings_list[["k"]] <- k
    } else if (cv_type == "lpo") {
      settings_list[["percentage"]] <- percentage
      settings_list[["number_folds"]] <- number_folds
    }
  }

  # Best model identification section
  if (include_best) {
   n_fit_scores <- ncol(result_df) - 1
   final_row <- c("Best")
  
   for (j in 2:ncol(result_df)) {
     colname <- names(result_df)[j]
     # Skip if it's not a metric column
     if (!any(sapply(c("mse", "mae", "dss", "crps", "scrps","wcrps","swcrps"), function(x) startsWith(colname, x)))) {
       final_row <- c(final_row, "")
       next
     }

      if (orientation_results == "negative") {
        best_tmp <- if(all(is.na(result_df[, j]))) {
          NA
        } else {
          which.min(result_df[, j])
        }
      } else {
        best_tmp <- if(all(is.na(result_df[, j]))) {
          NA
        } else {
          which.max(result_df[, j])
        }
      }
      final_row <- c(final_row, if(is.na(best_tmp)) NA else model_names[best_tmp])
   }
   result_df <- rbind(result_df, final_row)
   row.names(result_df)[nrow(result_df)] <- ""
  }

  # Cluster cleanup
  if (parallelize_RP) {
   parallel::stopCluster(cluster_tmp)
  }

  # Set return flag
  if (return_post_samples) {
   return_scores_folds <- TRUE
  }

  # Prepare output
  if (!return_scores_folds) {
   if (save_settings) {
     out <- list(
       scores_df = result_df,
       settings = settings_list
     )
     if (return_train_test) {
       out[["train_test"]] <- train_test_indexes
     }
     if(return_true_test_values){
       out[["true_test_values"]] <- true_test_values
     }
   } else if (return_train_test) {
     out <- list(scores_df = result_df, train_test = train_test_indexes)
     if(return_true_test_values){
       out[["true_test_values"]] <- true_test_values
     }      
   } else if(return_true_test_values){
     out <- list(scores_df = result_df, true_test_values = true_test_values)
   } else {
     out <- result_df
   }
  } else {
   # Add model names to all matrices in score lists
   add_model_names <- function(score_list) {
     lapply(score_list, function(mat) {
       colnames(mat) <- model_names
       return(mat)
     })
   }
  
   scores_folds <- list()
   if ("dss" %in% scores) scores_folds$dss <- add_model_names(dss)
   if ("mse" %in% scores) scores_folds$mse <- add_model_names(mse)
   if ("mae" %in% scores) scores_folds$mae <- add_model_names(mae)
   if ("crps" %in% scores) scores_folds$crps <- add_model_names(crps)
   if ("scrps" %in% scores) scores_folds$scrps <- add_model_names(scrps)
   if ("wcrps" %in% scores) scores_folds$wcrps <- add_model_names(wcrps)
   if ("swcrps" %in% scores) scores_folds$swcrps <- add_model_names(swcrps)   
  
   out <- list(
     scores_df = result_df,
     scores_folds = scores_folds
   )
  
   if (save_settings) {
     out[["settings"]] <- settings_list
   }
   if (return_train_test) {
     out[["train_test"]] <- train_test_indexes
   }
   if(return_true_test_values){
     out[["true_test_values"]] <- true_test_values
   }
   if (return_post_samples) {
     out[["post_samples"]] <- post_samples
   }
  }
  return(out)
}


#' @noRd 

sample_posterior_linear_predictor <- function(model, i_lik, test_list, n_samples, print){
        link_name <- model$.args$control.family[[i_lik]]$link

        model_family <- model$.args$family[[i_lik]]

        if (link_name == "default") {
          if (model_family %in% c("gaussian", "t")) {
            linkfuninv <- function(x) {
              x
            }
          } else if (model_family %in% c("gamma", "poisson", "stochvol", "stochvolln", "stochvolnig", "stochvolt")) {
            linkfuninv <- function(x) {
                exp(x)
              }
          } else if(model_family  == "binomial"){
              linkfuninv <- function(x) {
                exp(x)/(1 + exp(x))
              }
          } else{
          stop(paste("The family", model_family, "is not supported yet, please, raise an issue in https://github.com/davidbolin/rSPDE/issues requesting the support."))
          }
        } else {
          linkfuninv <- process_link(link_name)
        }

        formula_tmp <- process_formula_lhoods(model, i_lik)
        
        env_tmp <- environment(formula_tmp)
        assign("linkfuninv", linkfuninv, envir = env_tmp)
        if (model_family %in% c("stochvol", "stochvolln", "stochvolnig", "stochvolt")) {
          tmp_n_samples <- new_n_samples
          new_n_samples <- 2 * n_samples
        }
        if (print) {
          cat("Generating samples...\n")
        }

        data <- model$bru_info$lhoods[[i_lik]]$data

        post_samples <- inlabru::generate(model, newdata = data, formula = formula_tmp, n.samples = n_samples)

        return(post_samples[test_list[[i_lik]] , , drop=FALSE])
}


#' @noRd

get_posterior_samples <- function(post_linear_predictors, new_model, i_lik, new_n_samples, full_model, true_CV, print){
    model_family <- new_model$.args$family[[i_lik]]
    
    if(true_CV){
      model_sample <- new_model
    } else{
      model_sample <- full_model
    }

    family_mappings <- map_models_to_strings(full_model)
    if(family_mappings[[i_lik]][[1]] != ".none"){
      meas_err_par <- lapply(family_mappings[[i_lik]], function(param) {
                hyper_sample <- INLA::inla.hyperpar.sample(new_n_samples, model_sample, improve.marginals = TRUE)
                return(hyper_sample[,param])
            })
    }

    if (model_family == "gaussian") {
      sd_sample <- 1 / sqrt(as.vector(meas_err_par[[1]]))
      Y_sample <- lapply(1:nrow(post_linear_predictors), function(i) {
              post_linear_predictors[i, ] + sd_sample * rnorm(new_n_samples)
            }) 
    } else if(model_family == "gamma"){
      phi_sample <- as.vector(meas_err_par[[1]])
      Y_sample <-  lapply(1:nrow(post_linear_predictors), function(i) {
        scale_temp <- post_linear_predictors[i,] / phi_sample
        stats::rgamma(new_n_samples, shape = phi_sample, scale = scale_temp)
      })
    } else if(model_family == "t"){
      sd_sample <- 1 / sqrt(as.vector(meas_err_par[[1]]))
      deg_sample <- as.vector(meas_err_par[[2]])
      Y_sample <-  lapply(1:nrow(post_linear_predictors), function(i) {
        scale_temp <- post_linear_predictors[i,] + 
        sd_sample + rt(new_n_samples, df = deg_sample)
      })
    } else if(model_family == "poisson"){
      Y_sample <- lapply(1:nrow(post_linear_predictors), function(i) {
              stats::rpois(new_n_samples, post_linear_predictors[i,])
            })
    } else if(model_family == "binomial"){
      Y_sample <- lapply(1:nrow(post_linear_predictors), function(i) {
              stats::rbinom(n = new_n_samples, size = 1, prob = post_linear_predictors[i, ])
            })
    } else if(model_family == "stochvol"){
        phi_sample <- as.vector(meas_err_par[[1]])
        Y_sample <- lapply(1:nrow(post_linear_predictors), function(i) {
            sqrt(post_linear_predictors[i, ] + 1 / phi_sample) * rnorm(new_n_samples)
        })
    } else if(model_family == "stochvolln"){
        phi_sample <- as.vector(meas_err_par[[1]])
        mu_sample <- as.vector(meas_err_par[[2]])
        Y_sample <- lapply(1:nrow(post_linear_predictors), function(i) {
           var <- post_linear_predictors[i, ] + 1 / phi_sample
           mean <- mu_sample - 0.5 * var          
           mean + sqrt(var) * rnorm(new_n_samples)
          })
    } else if(model_family == "stochvolnig"){
        shape <- as.vector(meas_err_par[[1]])
        skewness <- as.vector(meas_err_par[[2]])
        gamma <- sqrt(1+skewness^2/shape^2)
        Y_sample <- lapply(1:nrow(post_linear_predictors), function(i) {
            sqrt(post_linear_predictors[i, ]) * GeneralizedHyperbolic::rnig(new_n_samples,
                                   mu = -skewness / gamma,   # mu_new corresponds to delta_old = -skewness_1 / gamma_1
                                   delta = shape / gamma,    # delta_new = sqrt(nu * sigma^2) = shape_1 / gamma_1
                                   alpha = sqrt(skewness^2 * gamma^4 + shape^2 * gamma^2),
                                   # alpha = sqrt(mu^2 / sigma^4 + nu / sigma^2)
                                   #        = sqrt(skewness_1^2 * gamma_1^4 + shape_1^2 * gamma_1^2)

                                   beta = skewness * gamma   # beta = mu_old / sigma^2 = skewness_1 * gamma_1
                              )
                              })        
    } else if(model_family == "stochvolt"){
      degree <- as.vector(meas_err_par[[1]])
      Y_sample <- lapply(1:nrow(post_linear_predictors), function(i) {
        sqrt(post_linear_predictors[i, ]) * rt(new_n_samples, degree)
      })      
    }
    if (print) {
      cat("Samples generated!\n")
    }
    
    return(Y_sample)
}

#' @noRd
map_models_to_strings <- function(models) {
 # Extract families from models
 families <- models$.args$family
 
 # Define base mappings
 mapping <- list(
   "gaussian" = "Precision for the Gaussian observations",
   "gamma" = "Precision-parameter for the Gamma observations", 
   "poisson" = ".none",
   "binomial" = ".none",
   "stochvol" = "Offset precision for stochvol",
   "stochvolln" = c("Offset precision for stochvolln","Mean offset for stochvolln"),
   "stochvolnig" = c("shape parameter for stochvol-nig", "skewness parameter for stochvol-nig"),
   "stochvolt" = "degrees of freedom for stochvol student-t",
   "t" = c("precision for the student-t observations", "degrees of freedom for student-t")
 )
 
 # Initialize result list
 result <- vector("list", length(families))
 
 # Process each family
 for (i in seq_along(families)) {
   family <- families[i]
   base_string <- mapping[[family]]
   
   # If it's first occurrence or .none, use base string
   if (i == 1 || base_string[1] == ".none") {
     result[[i]] <- base_string
   } else {
     # For subsequent occurrences, append [i]
     if (length(base_string) == 1) {
       result[[i]] <- paste0(base_string, "[", i, "]")
     } else {
       result[[i]] <- paste0(base_string, "[", i, "]")
     }
   }
 }
 
 return(result)
}

