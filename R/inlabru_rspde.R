
#'
#' @title rSPDE inlabru mapper
#' @name bru_get_mapper.inla_rspde
#' @param model An `inla_rspde` for which to construct or extract a mapper
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
#' \donttest{ #devel version
#' if (requireNamespace("INLA", quietly = TRUE) && 
#'      requireNamespace("inlabru", quietly = TRUE)){
#' library(INLA)
#' library(inlabru)
#' 
#' set.seed(123)
#' m <- 100
#' loc_2d_mesh <- matrix(runif(m * 2), m, 2)
#' mesh_2d <- inla.mesh.2d(
#'   loc = loc_2d_mesh,
#'   cutoff = 0.05,
#'   max.edge = c(0.1, 0.5)
#' )
#' sigma <- 1
#' range <- 0.2
#' nu <- 0.8
#' kappa <- sqrt(8 * nu) / range
#' op <- matern.operators(
#'   mesh = mesh_2d, nu = nu,
#'   range = range, sigma = sigma, m = 2,
#'   parameterization = "matern"
#' )
#' u <- simulate(op)
#' A <- inla.spde.make.A(
#'   mesh = mesh_2d,
#'   loc = loc_2d_mesh
#' )
#' sigma.e <- 0.1
#' y <- A %*% u + rnorm(m) * sigma.e
#' y <- as.vector(y)
#' 
#' data_df <- data.frame(y=y, x1 = loc_2d_mesh[,1],
#'                        x2 = loc_2d_mesh[,2])
#' coordinates(data_df) <- c("x1", "x2")
#' rspde_model <- rspde.matern(
#'   mesh = mesh_2d,
#'   nu_upper_bound = 2
#' )
#' 
#' cmp <- y ~ Intercept(1) + 
#'            field(coordinates, model = rspde_model)
#' 
#' 
#' rspde_fit <- bru(cmp, data = data_df)
#' summary(rspde_fit)
#' }
#' #devel.tag
#' }
bru_get_mapper.inla_rspde <- function(model,...) {
  mapper <- list(model = model)
  inlabru::bru_mapper_define(mapper, new_class = "bru_mapper_inla_rspde")
}

#' @param mapper A `bru_mapper_inla_rspde` object
#' @rdname bru_get_mapper.inla_rspde
ibm_n.bru_mapper_inla_rspde <- function(mapper, ...) {
  model <- mapper[["model"]]
  integer_nu <- model$integer.nu
  rspde_order <- model$rspde.order
  if(integer_nu){
            factor_rspde <- 1
  } else{
            factor_rspde <- rspde_order + 1
  }
  factor_rspde*model$n.spde
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

  if(model$est_nu){
    nu <- NULL
  } else{
   nu <- model$nu
  }

  rspde_order <- model$rspde.order
  rSPDE::rspde.make.A(mesh = model$mesh, loc=input,
                                rspde.order = rspde_order,
                                nu=nu)
}

#' @noRd 
# Function to process bru's formula

process_formula <- function(bru_result){
    form <- bru_result$bru_info$model$formula[3]
    form <- as.character(form)
    form <- strsplit(form, "f\\(")
    form <- form[[1]]
    form <- form[-1]
    form_proc <- sub(",.*", "",strsplit(form, "f\\(")[1])
    if(length(form)>1){
            for(i in 2:(length(form))){
                form_proc <- paste(form_proc, " + ",  sub(",.*", "",strsplit(form, "f\\(")[i]))
            }
    }
    form_proc <- paste("~", "linkfuninv(", form_proc, ")")
    return(stats::as.formula(form_proc))
}

#' @noRd 
# Function to process the link function

process_link <- function(link_name){
  return_link <- switch(link_name,
  "log" = function(x){INLA::inla.link.log(x, inverse=TRUE)},
  "invlog" = function(x){INLA::inla.link.invlog(x, inverse=TRUE)},
  "logit" = function(x){INLA::inla.link.logit(x, inverse=TRUE)},
   "invlogit" = function(x){INLA::inla.link.invlogit(x, inverse=TRUE)},
  "probit" = function(x){INLA::inla.link.probit(x, inverse=TRUE)},
 "invprobit" = function(x){INLA::inla.link.invprobit(x, inverse=TRUE)},
 "cloglog" = function(x){INLA::inla.link.cloglog(x, inverse=TRUE)},
 "invcloglog" = function(x){INLA::inla.link.invcloglog(x, inverse=TRUE)},
 "tan" = function(x){INLA::inla.link.tan(x, inverse=TRUE)},
 "invtan" = function(x){INLA::inla.link.invtan(x, inverse=TRUE)},
 "identity" = function(x){INLA::inla.link.identity(x, inverse=TRUE)},
 "invidentity" = function(x){INLA::inla.link.invidentity(x, inverse=TRUE)}
  )
  return(return_link)
}

#' @noRd 

bru_rerun_with_data <- function(result, idx_data, true_CV, fit_verbose) {
stopifnot(inherits(result, "bru"))
  if(!true_CV){
      options <- list(control.mode = list(restart = FALSE,
              theta=result$mode$theta, fixed = TRUE))
  } else{
    options <- list()
  }

  if(fit_verbose){
    options$verbose <- TRUE
  } else{
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

  lhoods_tmp <- info[["lhoods"]]
  lhoods_tmp[[1]]$response_data$BRU_response <- lhoods_tmp[[1]]$response_data$BRU_response[idx_data]
  
    # lhoods_tmp[[1]]$data <- lhoods_tmp[[1]]$data[idx_data,]

    lhoods_tmp[[1]]$data <- select_indexes(lhoods_tmp[[1]]$data, idx_data)

  if(length(lhoods_tmp[[1]]$E)>1){
    lhoods_tmp[[1]]$E <- lhoods_tmp[[1]]$E[idx_data]
  }

  if(length(lhoods_tmp[[1]]$Ntrials)>1){
    lhoods_tmp[[1]]$Ntrials <- lhoods_tmp[[1]]$Ntrials[idx_data]
  }

  if(length(lhoods_tmp[[1]]$weights)>1){
    lhoods_tmp[[1]]$weights <- lhoods_tmp[[1]]$weights[idx_data]
  }

  if(isS4(lhoods_tmp[[1]]$data)){
    lhoods_tmp[[1]]$drange <- lapply(lhoods_tmp[[1]]$data@data, function(i){range(i)})
  } else{
    lhoods_tmp[[1]]$drange <- lapply(lhoods_tmp[[1]]$data, function(i){range(i)})
  }
  
  # Get the components list

  list_of_components <- names(info[["model"]][["effects"]])

  backup_list <- list()

  total_length <- NULL
  small_length <- length(idx_data)

  for(comp in list_of_components){
    name_input_group <- info[["model"]][["effects"]][[comp]][["group"]][["input"]][["input"]]
    if(!is.null(name_input_group)){
      name_input_group <- as.character(name_input_group)
      comp_group_tmp <-  info[["model"]][["effects"]][[comp]][["env"]][[name_input_group]]
      if(is.null(total_length) && !is.null(comp_group_tmp)){
        total_length <- length(comp_group_tmp)

        if(length(comp_group_tmp) == total_length){
          backup_list[[comp]][["group_val"]] <- info[["model"]][["effects"]][[comp]][["env"]][[name_input_group]]
          comp_group_tmp <- comp_group_tmp[idx_data]
          if(!is.null(comp_group_tmp)){
            assign(name_input_group, comp_group_tmp, envir = info[["model"]][["effects"]][[comp]][["env"]])
          }
        }
      }
    }
    name_input_repl <- info[["model"]][["effects"]][[comp]][["replicate"]][["input"]][["input"]]
    if(!is.null(name_input_repl)){
      name_input_repl <- as.character(name_input_repl)
      comp_repl_tmp <-  info[["model"]][["effects"]][[comp]][["env"]][[name_input_repl]]
      if(is.null(total_length) && !is.null(comp_repl_tmp)){
        total_length <- length(comp_repl_tmp)
        if(length(comp_repl_tmp) == total_length){
          backup_list[[comp]][["repl_val"]] <- info[["model"]][["effects"]][[comp]][["env"]][[name_input_repl]]
          comp_repl_tmp <- comp_repl_tmp[idx_data]
          if(!is.null(comp_repl_tmp)){
            assign(name_input_repl, comp_repl_tmp, envir = info[["model"]][["effects"]][[comp]][["env"]])
          }
        }
      }
    }
  }


  if(!true_CV){
  result <- inlabru::iinla(
    model = info[["model"]],
    lhoods = lhoods_tmp,
    options = info[["options"]]
  )
  } else{
      result <- inlabru::iinla(
    model = info[["model"]],
    lhoods = lhoods_tmp,
        initial = result,
    options = info[["options"]]
  )
  }

  # Assigning back:

  for(comp in list_of_components){
    name_input_group <- info[["model"]][["effects"]][[comp]][["group"]][["input"]][["input"]]
    if(!is.null(name_input_group)){
      name_input_group <- as.character(name_input_group)
      if(!is.null(backup_list[[comp]][["group_val"]])){
        assign(name_input_group, backup_list[[comp]][["group_val"]], envir = info[["model"]][["effects"]][[comp]][["env"]])
      }
    }
    name_input_repl <- info[["model"]][["effects"]][[comp]][["replicate"]][["input"]][["input"]]
    if(!is.null(name_input_repl)){
      name_input_repl <- as.character(name_input_repl)
      if(!is.null(backup_list[[comp]][["repl_val"]])){
        assign(name_input_repl, backup_list[[comp]][["repl_val"]], envir = info[["model"]][["effects"]][[comp]][["env"]])
      }
    }
  }


  timing_end <- Sys.time()
  result$bru_timings <-
    rbind(
      original_timings[1, , drop = FALSE],
      result[["bru_iinla"]][["timings"]]
    )

  # Add bru information to the result
  result$bru_info <- info
  class(result) <- c("bru", class(result))
  return(result)
}


#' @noRd 

get_post_var <- function(density_df){
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
        denstemp(z) * 1/z
      }, lower = min_x, upper = max_x,
      subdivisions = nrow(density_df),
                  stop.on.error = FALSE
    )$value

    return(post_var)
}

#' @noRd 

prepare_df_pred <- function(df_pred, result, idx_test){
  info <- result[["bru_info"]]
  list_of_components <- names(info[["model"]][["effects"]])
  lhoods_tmp <- info[["lhoods"]]

  for(comp in list_of_components){
    name_input_group <- info[["model"]][["effects"]][[comp]][["group"]][["input"]][["input"]]
    if(!is.null(name_input_group)){
      name_input_group <- as.character(name_input_group)
      comp_group_tmp <-  info[["model"]][["effects"]][[comp]][["env"]][[name_input_group]]
      if(!is.null(comp_group_tmp)){
        if(!is.null(dim(comp_group_tmp))){
          comp_group_tmp <- comp_group_tmp[idx_test, , drop=FALSE]
        } else{
          comp_group_tmp <- comp_group_tmp[idx_test]
        }
      } else{
        if(!is.null(dim(lhoods_tmp[[1]]$data[[name_input_group]]))){
          comp_group_tmp <- lhoods_tmp[[1]]$data[[name_input_group]][idx_test, , drop=FALSE]
        } else{
          comp_group_tmp <- lhoods_tmp[[1]]$data[[name_input_group]][idx_test]
        }
      }
      df_pred[[name_input_group]] <- comp_group_tmp
    }
    name_input_repl <- info[["model"]][["effects"]][[comp]][["replicate"]][["input"]][["input"]]
    if(!is.null(name_input_repl)){
      name_input_repl <- as.character(name_input_repl)
      comp_repl_tmp <-  info[["model"]][["effects"]][[comp]][["env"]][[name_input_repl]]
      if(!is.null(comp_repl_tmp)){
        if(!is.null(dim(comp_repl_tmp))){
          comp_repl_tmp <- comp_repl_tmp[idx_test, ,drop=FALSE]
        } else{
          comp_repl_tmp <- comp_repl_tmp[idx_test]
        }
      } else{
        if(!is.null(dim(lhoods_tmp[[1]]$data[[name_input_repl]]))){
          comp_repl_tmp <- lhoods_tmp[[1]]$data[[name_input_repl]][idx_test, , drop=FALSE]
        } else{
          comp_repl_tmp <- lhoods_tmp[[1]]$data[[name_input_repl]][idx_test]
        }
      }
        df_pred[[name_input_repl]] <- comp_repl_tmp
    }
  }
  return(df_pred)
}


#' @name cross_validation
#' @title Perform cross-validation on a list of fitted models.
#' @description Obtain several scores for a list of fitted models according 
#' to a folding scheme.
#' @param models A fitted model obtained from calling the `bru()` function or a list of models fitted with the `bru()` function.
#' @param model_names A vector containing the names of the models to appear in the returned `data.frame`. If `NULL`, the names will be of the form `Model 1`, `Model 2`, and so on. By default, it will try to obtain the name from the models list.
#' @param scores A vector containing the scores to be computed. The options are "mse", "crps", "scrps" and "dss". By default, all scores are computed.
#' @param cv_type The type of the folding to be carried out. The options are `k-fold` for `k`-fold cross-validation, in which case the parameter `k` should be provided, 
#' `loo`, for leave-one-out and `lpo` for leave-percentage-out, in this case, the parameter `percentage` should be given, and also the `number_folds` 
#' with the number of folds to be done. The default is `k-fold`.
#' @param k The number of folds to be used in `k`-fold cross-validation. Will only be used if `cv_type` is `k-fold`.
#' @param percentage The percentage (from 1 to 99) of the data to be used to train the model. Will only be used if `cv_type` is `lpo`.
#' @param number_folds Number of folds to be done if `cv_type` is `lpo`.
#' @param n_samples Number of samples to compute the posterior statistics to be used to compute the scores.
#' @param return_scores_folds If `TRUE`, the scores for each fold will also be returned.
#' @param orientation_results character vector. The options are "negative" and "positive". If "negative", the smaller the scores the better. If "positive", the larger the scores the better.
#' @param include_best Should a row indicating which model was the best for each score be included?
#' @param train_test_indexes A list containing two entries `train`, which is a list whose elements are vectors of indexes of the training data, and `test`, which is a list whose elements are vectors of indexes of the test data.
#' Typically this will be returned list obtained by setting the argument `return_train_test` to `TRUE`.
#' @param return_train_test Logical. Should the training and test indexes be returned? If 'TRUE' the train and test indexes will the 'train_test' element of the returned list.
#' @param return_post_samples If `TRUE` the posterior samples will be included in the returned list.
#' @param parallelize_RP Logical. Should the computation of CRPS and SCRPS be parallelized?
#' @param n_cores_RP Number of cores to be used if `parallelize_rp` is `TRUE`.
#' @param true_CV Should a `TRUE` cross-validation be performed? If `TRUE` the models will be fitted on the training dataset. If `FALSE`, the parameters will be kept fixed at the ones obtained in the result object.
#' @param save_settings Logical. If `TRUE`, the settings used in the cross-validation will also be returned.
#' @param print Should partial results be printed throughout the computation?
#' @param fit_verbose Should INLA's run during cross-validation be verbose?
#' @return A data.frame with the fitted models and the corresponding scores.
#' @export

cross_validation <- function(models, model_names = NULL, scores = c("mse", "crps", "scrps", "dss"),
                              cv_type = c("k-fold", "loo", "lpo"),
                              k = 5, percentage = 20, number_folds = 10,
                              n_samples = 1000, return_scores_folds = FALSE,
                              orientation_results = c("negative", "positive"),
                              include_best = TRUE,
                              train_test_indexes = NULL,
                              return_train_test = FALSE,
                              return_post_samples = FALSE,                              
                              parallelize_RP = FALSE, n_cores_RP = parallel::detectCores()-1,
                              true_CV = FALSE, save_settings = FALSE, 
                              print = TRUE,
                              fit_verbose = FALSE){

                                orientation_results <- orientation_results[[1]]
                                if(!(orientation_results %in% c("positive", "negative"))){
                                  stop("orientation_results must be either 'positive' or 'negative'!")
                                }

                                scores <- intersect(scores, c("mse", "crps", "scrps", "dss"))

                                cv_type <- cv_type[[1]]
                                if(!(cv_type %in% c("k-fold", "loo", "lpo"))){
                                  stop("The possible options for cv_type are 'k-fold', 'loo' or 'lpo'!")
                                }

                                if(!is.numeric(percentage)){
                                  stop("percentage must be numeric!")
                                }

                                if(percentage %%1 != 0){
                                  warning("Non-integer percentage given, it will be rounded to an integer number.")
                                  percentage <- round(percentage)
                                }

                                if(percentage <= 0 || percentage >= 100){
                                  stop("percentage must be a number between 1 and 99!")
                                }

                                if(!is.numeric(number_folds)){
                                  stop("number_folds must be numeric!")
                                }

                                if(number_folds %% 1 != 0){
                                  warning("Non-integer number_folds given, it will be rounded to an integer number.")
                                  number_folds <- round(number_folds)
                                }

                                if(number_folds <= 0){
                                  stop("number_folds must be positive!")
                                }

                                if(inherits(models, "bru")){
                                  models <- list(models)
                                } else{
                                  for(i in 1:length(models)){
                                    if(!inherits(models[[i]],"bru")){
                                      stop("models must be either a result from a bru call or a list of results from bru() calls!")
                                    }
                                  }
                                }                                


                                if(!is.numeric(n_samples)){
                                  stop("n_samples must be numeric!")
                                }

                                if(n_samples %% 1 != 0){
                                  warning("Non-integer n_samples given, it will be rounded to an integer number.")
                                  n_samples <- round(n_samples)
                                }

                                if(n_samples <= 0){
                                  stop("n_samples must be positive!")
                                }

                                if(parallelize_RP){
                                      cluster_tmp <-  parallel::makeCluster(n_cores_RP)
                                      doParallel::registerDoParallel(cluster_tmp)
                                }                                

                                # Getting the data if NULL
                                data <- models[[1]]$bru_info$lhoods[[1]]$data

                                if(is.vector(data)){
                                  data <- as.data.frame(data)
                                }

                                # Creating lists of train and test datasets

                                if(is.null(train_test_indexes)){
                                  train_test_indexes <- create_train_test_indices(data, cv_type = cv_type, 
                                            k = k, percentage = percentage, number_folds = number_folds)
                                  train_list <- train_test_indexes[["train"]]
                                  test_list <- train_test_indexes[["test"]]
                                } else{
                                  if(!is.list(train_test_indexes)){
                                    stop("train_test_indexes should be a list!")
                                  }
                                  if(is.null(train_test_indexes[["train"]])){
                                    stop("train_test_indexes must contain a 'train' element.")
                                  }
                                  if(is.null(train_test_indexes[["test"]])){
                                    stop("train_test_indexes must contain a 'test' element.")
                                  }
                                  if(!is.list(train_test_indexes[["train"]])){
                                    stop("train_test_indexes$train must be a list!")
                                  }
                                  if(!is.list(train_test_indexes[["test"]])){
                                    stop("train_test_indexes$test must be a list!")
                                  }
                                  train_list <- train_test_indexes[["train"]]
                                  test_list <- train_test_indexes[["test"]]
                                }

                                post_samples <- list()
                                hyper_samples <- list()                                

                                for(model_number in 1:length(models)){
                                    post_samples[[model_names[[model_number]]]] <- vector(mode = "list", length = length(train_list))
                                    hyper_samples[[model_names[[model_number]]]] <- vector(mode = "list", length = 2)                                    
                                    hyper_samples[[model_names[[model_number]]]][[1]] <- vector(mode = "list", length = length(train_list))
                                    hyper_samples[[model_names[[model_number]]]][[2]] <- vector(mode = "list", length = length(train_list))
                                }
                                # Perform the cross-validation
                                
                                result_df <- data.frame(Model = model_names)

                                dss <- matrix(numeric(length(train_list)*length(models)), ncol = length(models))
                                mse <- matrix(numeric(length(train_list)*length(models)), ncol = length(models))
                                crps <- matrix(numeric(length(train_list)*length(models)), ncol = length(models))
                                scrps <- matrix(numeric(length(train_list)*length(models)), ncol = length(models))

                                # Get the formulas for the models

                                formula_list <- lapply(models, function(model){process_formula(model)})

                                if(("crps" %in% scores) || ("scrps" %in% scores)){
                                  new_n_samples <- 2 * n_samples
                                } else{
                                  new_n_samples <- n_samples
                                }

                                for(fold in 1:length(train_list)){                                 
                                  for(model_number in 1:length(models)){

                                    test_data <- models[[model_number]]$bru_info$lhoods[[1]]$response_data$BRU_response[test_list[[fold]]]

                                    link_name <- models[[model_number]]$.args$control.family[[1]]$link

                                        if(link_name == "default"){
                                          if(models[[model_number]]$.args$family == "gaussian"){     
                                            linkfuninv <- function(x){x}
                                          } else if(models[[model_number]]$.args$family == "gamma"){
                                            linkfuninv <- function(x){exp(x)}
                                          } else if (models[[model_number]]$.args$family == "poisson"){
                                            linkfuninv <- function(x){exp(x)}
                                          }
                                        } else {
                                          linkfuninv <- process_link(link_name)
                                        } 

                                        formula_tmp <- formula_list[[model_number]]
                                        env_tmp <- environment(formula_tmp)
                                        assign("linkfuninv", linkfuninv, envir = env_tmp)

                                        post_predict <- group_predict(models = models[[model_number]], model_names = model_names[[model_number]],
                                                              formula = formula_tmp, train_indices = train_list[[fold]],
                                                              test_indices = test_list[[fold]], n_samples = new_n_samples,
                                                              pseudo_predict = !true_CV, return_samples = TRUE, return_hyper_samples = TRUE,
                                                              n_hyper_samples = 2, compute_posterior_means = TRUE, print = print, fit_verbose = fit_verbose)

                                        hyper_marginals <- post_predict[["hyper_marginals"]][[model_names[[model_number]]]][[1]]
                                        hyper_summary <- post_predict[["hyper_summary"]][[model_names[[model_number]]]][[1]]
                                        hyper_samples_1 <- post_predict[["hyper_samples"]][[model_names[[model_number]]]][[1]][[1]]
                                        hyper_samples_2 <- post_predict[["hyper_samples"]][[model_names[[model_number]]]][[2]][[1]]
                                        posterior_samples <- post_predict[["post_samples"]][[model_names[[model_number]]]][[1]]
                                        posterior_mean <- post_predict[["post_means"]][[model_names[[model_number]]]][[1]]

                                        if(return_post_samples){
                                          post_samples[[model_names[[model_number]]]][[fold]] <- posterior_samples
                                          hyper_samples[[model_names[[model_number]]]][[1]][[fold]] <- hyper_samples_1
                                          hyper_samples[[model_names[[model_number]]]][[2]][[fold]] <- hyper_samples_2
                                        }                           


                                        if("mse" %in% scores){
                                          mse[fold, model_number] <- mean((test_data - posterior_mean)^2)          
                                          if(orientation_results == "positive"){
                                            mse[fold, model_number] <- - mse[fold, model_number] 
                                          }               
                                         if(print){
                                            cat(paste("MSE:",mse[fold, model_number],"\n"))
                                          }                                                        
                                        }



                                      if(models[[model_number]]$.args$family == "gaussian"){                                        
                                        if("dss" %in% scores){
                                          density_df <- hyper_marginals$`Precision for the Gaussian observations`
                                          Expect_post_var <- tryCatch(get_post_var(density_df), error = function(e) NA)
                                          if(is.na(Expect_post_var)){
                                            Expect_post_var <- 1/hyper_summary["Precision for the Gaussian observations","mean"]
                                          }    

                                          posterior_variance_of_mean <- rowMeans(posterior_samples^2) - posterior_mean^2
                                          post_var <- Expect_post_var + posterior_variance_of_mean                                          

                                          dss[fold, model_number] <- mean((test_data - posterior_mean)^2/post_var + log(post_var))
                                          if(orientation_results == "positive"){
                                            dss[fold, model_number] <- - dss[fold, model_number]
                                          }
                                          if(print){
                                            cat(paste("DSS:", dss[fold, model_number],"\n"))
                                          }
                                        }

                                        if(("crps" %in% scores) || ("scrps" %in% scores)){
                                          phi_sample_1 <- as.vector(hyper_samples_1[,"Precision for the Gaussian observations"])
                                          sd_sample_1 <- 1/sqrt(phi_sample_1)  

                                          phi_sample_2 <- as.vector(hyper_samples_2[,"Precision for the Gaussian observations"])
                                          sd_sample_2 <- 1/sqrt(phi_sample_2)  

                                          if(parallelize_RP){
                                            Y1_sample <- foreach::`%dopar%`(foreach::foreach(i = 1:length(test_data)), {
                                                            posterior_samples[i,1:n_samples] + sd_sample_1 * rnorm(n_samples)
                                                          })
                                            Y2_sample <- foreach::`%dopar%`(foreach::foreach(i = 1:length(test_data)), {
                                                            posterior_samples[i,(n_samples+1):(2*n_samples)] + sd_sample_2 * rnorm(n_samples)
                                                          })
                                            E1_tmp <- foreach::`%dopar%`(foreach::foreach(i = 1:length(test_data)), {
                                                            mean(abs(Y1_sample[[i]]-test_data[i]))
                                                          })
                                            E2_tmp <- foreach::`%dopar%`(foreach::foreach(i = 1:length(test_data)), {
                                                            mean(abs(Y1_sample[[i]]-Y2_sample[[i]]))
                                                          })
                                            
                                          } else{
                                            Y1_sample <- lapply(1:length(test_data), function(i){
                                                      posterior_samples[i,1:n_samples] + sd_sample_1 * rnorm(n_samples)
                                                    })
                                            Y2_sample <- lapply(1:length(test_data), function(i){
                                                      posterior_samples[i,(n_samples+1):(2*n_samples)] + sd_sample_2 * rnorm(n_samples)
                                                    })
                                            E1_tmp <- lapply(1:length(test_data), function(i){
                                                      mean(abs(Y1_sample[[i]]-test_data[i]))
                                                    })
                                            E2_tmp <- lapply(1:length(test_data), function(i){
                                                      mean(abs(Y1_sample[[i]]-Y2_sample[[i]]))
                                                    })              
                                          }
                                    
                                          if("crps" %in% scores){
                                              
                                              crps_temp <- lapply(1:length(test_data), function(i){
                                                return(E1_tmp[[i]] - 0.5*E2_tmp[[i]])
                                              })

                                              crps_temp <- unlist(crps_temp)
                                              crps[fold, model_number] <- mean(crps_temp)  
                                              if(orientation_results == "negative"){
                                                crps[fold, model_number] <- - crps[fold, model_number]
                                              }    

                                            if(print){
                                              cat(paste("CRPS:",crps[fold, model_number],"\n"))
                                            }      
                                          }
                                          if("scrps" %in% scores){
                                              scrps_temp <- lapply(1:length(test_data), function(i){
                                                return(-E1_tmp[[i]]/E2_tmp[[i]] - 0.5*log(E2_tmp[[i]]))
                                              })
                                              scrps_temp <- unlist(scrps_temp)
                                              scrps[fold, model_number] <- mean(scrps_temp)  
                                              if(orientation_results == "negative"){
                                                scrps[fold, model_number] <- - scrps[fold, model_number]
                                              }   

                                            if(print){
                                              cat(paste("SCRPS:",scrps[fold, model_number],"\n"))
                                            }                                                  
                                          }

                                        }

                                      } else if (models[[model_number]]$.args$family == "gamma"){
                                        if("dss" %in% scores){
                                          Expected_post_var <- hyper_marginals["Precision parameter for the Gamma observations","mean"]/(posterior_mean^2)
                                          posterior_variance_of_mean <- rowMeans(posterior_samples^2) - posterior_mean^2

                                          post_var <- Expected_post_var + posterior_variance_of_mean                                          
                                          dss[fold, model_number] <- mean((test_data - posterior_mean)^2/post_var + log(post_var))
                                          if(orientation_results == "positive"){
                                            dss[fold, model_number] <- - dss[fold, model_number]
                                          }                                          
                                          if(print){
                                            cat(paste("DSS:", dss[fold, model_number],"\n"))
                                          }
                                        }

                                        if(("crps" %in% scores) || ("scrps" %in% scores)){
                                          phi_sample_1 <- as.vector(hyper_samples_1[,"Precision parameter for the Gamma observations"])

                                          phi_sample_2 <- as.vector(hyper_samples_2[,"Precision parameter for the Gamma observations"])

                                          if(parallelize_RP){
                                            Y1_sample <- foreach::`%dopar%`(foreach::foreach(i = 1:length(test_data)), {
                                                            scale_temp <- posterior_samples[i,1:n_samples] / phi_sample_1
                                                            stats::rgamma(n_samples, shape = phi_sample_1, scale = scale_temp)
                                                          })
                                            Y2_sample <- foreach::`%dopar%`(foreach::foreach(i = 1:length(test_data)), {
                                                            scale_temp <- posterior_samples[i,(n_samples+1):(2*n_samples)] / phi_sample_2
                                                            stats::rgamma(n_samples, shape = phi_sample_2, scale = scale_temp)
                                                          })
                                            E1_tmp <- foreach::`%dopar%`(foreach::foreach(i = 1:length(test_data)), {
                                                            mean(abs(Y1_sample[[i]]-test_data[i]))
                                                          })
                                            E2_tmp <- foreach::`%dopar%`(foreach::foreach(i = 1:length(test_data)), {
                                                            mean(abs(Y1_sample[[i]]-Y2_sample[[i]]))
                                                          })
                                            
                                          } else{
                                            Y1_sample <- lapply(1:length(test_data), function(i){
                                                            scale_temp <- posterior_samples[i,1:n_samples] / phi_sample_1
                                                            stats::rgamma(n_samples, shape = phi_sample_1, scale = scale_temp)
                                                    })
                                            Y2_sample <- lapply(1:length(test_data), function(i){
                                                            scale_temp <- posterior_samples[i,(n_samples+1):(2*n_samples)] / phi_sample_2
                                                            stats::rgamma(n_samples, shape = phi_sample_2, scale = scale_temp)
                                                    })
                                            E1_tmp <- lapply(1:length(test_data), function(i){
                                                      mean(abs(Y1_sample[[i]]-test_data[i]))
                                                    })
                                            E2_tmp <- lapply(1:length(test_data), function(i){
                                                      mean(abs(Y1_sample[[i]]-Y2_sample[[i]]))
                                                    })              
                                          }
                                    
                                          if("crps" %in% scores){
                                              crps_temp <- lapply(1:length(test_data), function(i){
                                                return(E1_tmp[[i]] - 0.5*E2_tmp[[i]])
                                              })

                                              crps_temp <- unlist(crps_temp)
                                              crps[fold, model_number] <- mean(crps_temp)  
                                              if(orientation_results == "negative"){
                                                crps[fold, model_number] <- crps[fold, model_number]
                                              }    

                                            if(print){
                                              cat(paste("CRPS:",crps[fold, model_number],"\n"))
                                            }      
                                          }
                                          if("scrps" %in% scores){
                                              scrps_temp <- lapply(1:length(test_data), function(i){
                                                return(-E1_tmp[[i]]/E2_tmp[[i]] - 0.5*log(E2_tmp[[i]]))
                                              })
                                              scrps_temp <- unlist(scrps_temp)
                                              scrps[fold, model_number] <- mean(scrps_temp)  
                                              if(orientation_results == "negative"){
                                                scrps[fold, model_number] <- - scrps[fold, model_number]
                                              }   

                                        }
                                        }
                                     } else if (models[[model_number]]$.args$family == "poisson"){
                                      
                                        if("dss" %in% scores){
                                          posterior_variance_of_mean <- rowMeans(posterior_samples^2) - posterior_mean^2
                                          post_var <- posterior_mean + posterior_variance_of_mean

                                          dss[fold, model_number] <- mean((test_data - posterior_mean)^2/post_var + log(post_var))
                                          if(orientation_results == "positive"){
                                            dss[fold, model_number] <- - dss[fold, model_number]
                                          }                                                 
                                          if(print){
                                            cat(paste("DSS:", dss[fold, model_number],"\n"))
                                          }
                                        }

                                     if(("crps" %in% scores) || ("scrps" %in% scores)){

                                          if(parallelize_RP){
                                            Y1_sample <- foreach::`%dopar%`(foreach::foreach(i = 1:length(test_data)), {
                                                            stats::rpois(n_samples, posterior_samples[i,1:n_samples])
                                                          })
                                            Y2_sample <- foreach::`%dopar%`(foreach::foreach(i = 1:length(test_data)), {
                                                            stats::rpois(n_samples, posterior_samples[i,(n_samples+1):(2*n_samples)])
                                                          })
                                            E1_tmp <- foreach::`%dopar%`(foreach::foreach(i = 1:length(test_data)), {
                                                            mean(abs(Y1_sample[[i]]-test_data[i]))
                                                          })
                                            E2_tmp <- foreach::`%dopar%`(foreach::foreach(i = 1:length(test_data)), {
                                                            mean(abs(Y1_sample[[i]]-Y2_sample[[i]]))
                                                          })
                                            
                                          } else{
                                            Y1_sample <- lapply(1:length(test_data), function(i){
                                                            stats::rpois(n_samples, posterior_samples[i,1:n_samples])
                                                    })
                                            Y2_sample <- lapply(1:length(test_data), function(i){
                                                            stats::rpois(n_samples, posterior_samples[i,(n_samples+1):(2*n_samples)])
                                                    })
                                            E1_tmp <- lapply(1:length(test_data), function(i){
                                                      mean(abs(Y1_sample[[i]]-test_data[i]))
                                                    })
                                            E2_tmp <- lapply(1:length(test_data), function(i){
                                                      mean(abs(Y1_sample[[i]]-Y2_sample[[i]]))
                                                    })              
                                          }
                                    
                                          if("crps" %in% scores){
                                              
                                              crps_temp <- lapply(1:length(test_data), function(i){
                                                return(E1_tmp[[i]] - 0.5*E2_tmp[[i]])
                                              })

                                              crps_temp <- unlist(crps_temp)
                                              crps[fold, model_number] <- mean(crps_temp)  
                                              if(orientation_results == "negative"){
                                                crps[fold, model_number] <- crps[fold, model_number]
                                              }    

                                            if(print){
                                              cat(paste("CRPS:",crps[fold, model_number],"\n"))
                                            }      
                                          }
                                          if("scrps" %in% scores){
                                              scrps_temp <- lapply(1:length(test_data), function(i){
                                                return(-E1_tmp[[i]]/E2_tmp[[i]] - 0.5*log(E2_tmp[[i]]))
                                              })
                                              scrps_temp <- unlist(scrps_temp)
                                              scrps[fold, model_number] <- mean(scrps_temp)  
                                              if(orientation_results == "negative"){
                                                scrps[fold, model_number] <- - scrps[fold, model_number]
                                              }
                                          }   

                                        }
                                     }

                                  } 

                                }
                                

                                if("dss" %in% scores){
                                  dss_mean <- colMeans(dss)
                                  result_df <- data.frame(result_df, dss = dss_mean)
                                } 
                                if("mse" %in% scores){
                                  mse_mean <- colMeans(mse)
                                  result_df <- data.frame(result_df, mse = mse_mean)
                                }
                                if("crps" %in% scores){
                                  crps_mean <- colMeans(crps)
                                  result_df <- data.frame(result_df, crps = crps_mean)
                                }

                                if("scrps" %in% scores){
                                  scrps_mean <- colMeans(scrps)
                                  result_df <- data.frame(result_df, scrps = scrps_mean)
                                }

        if(save_settings){
          settings_list <- list(n_samples = n_samples, cv_type = cv_type, true_CV = true_CV,
                                  orientation_results = orientation_results)
          if(cv_type == "k-fold"){
            settings_list[["k"]] <- k
          } else if(cv_type == "lpo"){
            settings_list[["percentage"]] <- percentage
            settings_list[["number_folds"]] <- number_folds
          }
        }

        if(include_best){
          n_fit_scores <- ncol(result_df)-1
          final_row <- c("Best")
          for(j in 2:ncol(result_df)){
            if(orientation_results == "negative"){
              best_tmp <- which.min(result_df[,j])
              final_row <- c(final_row, model_names[best_tmp])
            } else{
              best_tmp <- which.max(result_df[,j])
              final_row <- c(final_row, model_names[best_tmp])
            }
          }
          result_df <- rbind(result_df, final_row)
          row.names(result_df)[nrow(result_df)] <- ""
        }


        if(parallelize_RP){
              parallel::stopCluster(cluster_tmp)
        }
        
        if(!return_scores_folds){
            if(save_settings){
              out <- list(scores_df = result_df,
                      settings = settings_list)
              if(return_train_test){
                out[["train_test"]] <- list(train = train_list, test = test_list)
              }
            } else if(return_train_test){
              out <- list(scores_df = result_df, train_test = list(train = train_list, test = test_list))
            } else{
              out <- result_df
            }
        } else{
          colnames(dss) <- model_names
          colnames(mse) <- model_names
          colnames(crps) <- model_names
          colnames(scrps) <- model_names
          out <- list(scores_df = result_df,
                      scores_folds = list(dss = dss, mse = mse, crps = crps, scrps = scrps))
          if(save_settings){
             out[["settings"]] <- settings_list
            }
          if(return_train_test){
            out[["train_test"]] <- list(train = train_list, test = test_list)
          }
        }

   if(return_post_samples){
    out[["post_samples"]] <- post_samples
    out[["hyper_samples"]] <- hyper_samples
   }

  return(out)
}




#' @name group_predict
#' @title Perform prediction on a testing set based on a training set
#' @description Compute prediction of a formula-based expression on a testing set based on a training set.
#' @param models A fitted model obtained from calling the `bru()` function or a list of models fitted with the `bru()` function.
#' @param model_names A vector containing the names of the models to appear in the returned `data.frame`. If `NULL`, the names will be of the form `Model 1`, `Model 2`, and so on. By default, it will try to obtain the name from the models list.
#' @param formula A formula where the right hand side defines an R expression to evaluate for each generated sample. If `NULL``, the latent and hyperparameter states are returned as named list elements. See the manual for the `predict` method in the `inlabru` package.
#' @param train_indices A list containing the indices of the observations for the model to be trained, or a numerical vector containing the indices.
#' @param test_indices A list containing the indices of the test data, where the prediction will be done, or a numerical vector containing the indices.
#' @param n_samples Number of samples to compute the posterior statistics to be used to compute the scores.
#' @param pseudo_predict If `TRUE`, the models will NOT be refitted on the training data, and the parameters obtained on the entire dataset will be used. If `FALSE`, the models will be refitted on the training data.
#' @param return_samples Should the posterior samples be returned?
#' @param return_hyper_samples Should samples for the hyperparameters be obtained?
#' @param n_hyper_samples Number of independent samples of the hyper parameters of size `n_samples`.
#' @param compute_posterior_means Should the posterior means be computed from the posterior samples?
#' @param print Should partial results be printed throughout the computation?
#' @param fit_verbose Should INLA's run during the prediction be verbose?
#' @return A data.frame with the fitted models and the corresponding scores.
#' @export

group_predict <- function(models, model_names = NULL, formula = NULL,
        train_indices, test_indices, n_samples = 1000, 
        pseudo_predict = TRUE, 
        return_samples = FALSE, return_hyper_samples = FALSE,
        n_hyper_samples = 1,
        compute_posterior_means = TRUE,
        print = TRUE, fit_verbose = FALSE){

                                if(length(train_indices) != length(test_indices)){
                                  if(!is.numeric(train_indices) || !is.numeric(test_indices)){
                                    stop("train_indices and test_indices must be lists of the same length or must be numerical vectors containing the indices!")
                                  }
                                }

                                if(is.numeric(train_indices) && is.numeric(test_indices)){
                                  train_indices <- list(train_indices)
                                  test_indices <- list(test_indices)
                                }
          
                                if(!is.numeric(n_samples)){
                                  stop("n_samples must be numeric!")
                                }

                                if(n_samples %% 1 != 0){
                                  warning("Non-integer n_samples given, it will be rounded to an integer number.")
                                  n_samples <- round(n_samples)
                                }

                                if(n_samples <= 0){
                                  stop("n_samples must be positive!")
                                }

                                if(!is.list(models)){
                                  stop("models should either be a result from a bru() call or a list of results from bru() calls!")
                                }
                                if(inherits(models, "bru")){
                                  models <- list(models)
                                } else{
                                  for(i in 1:length(models)){
                                    if(!inherits(models[[i]],"bru")){
                                      stop("models must be either a result from a bru call or a list of results from bru() calls!")
                                    }
                                  }
                                }

                                if(is.null(model_names) && is.list(models)){
                                  model_names <- names(models)
                                }

                                if(!is.null(model_names)){
                                  if(!is.character(model_names)){
                                    stop("model_names must be a vector of strings!")
                                  }
                                  if(length(models)!= length(model_names)){
                                    stop("model_names must contain one name for each model!")
                                  }
                                } else{
                                  model_names <- vector(mode = "character", length(models))
                                  for(i in 1:length(models)){
                                    model_names[i] <- paste("Model",i)
                                  }
                                }

                                # Getting the data if NULL
                                data <- models[[1]]$bru_info$lhoods[[1]]$data

                                if(is.vector(data)){
                                  data <- as.data.frame(data)
                                }

                            post_samples <- list()
                            post_means <- list()
                            hyper_samples <- list()
                            hyper_marginals <- list()
                            hyper_summary <- list()
                            
                            for(model_number in 1:length(models)){
                                  post_samples[[model_names[[model_number]]]] <- vector(mode = "list", length = length(train_indices))
                                  post_means[[model_names[[model_number]]]] <- vector(mode = "list", length = length(train_indices))                                  
                                  hyper_samples[[model_names[[model_number]]]] <- vector(mode = "list", length = n_hyper_samples)
                                  hyper_marginals[[model_names[[model_number]]]] <- vector(mode = "list", length = length(train_indices))
                                  hyper_summary[[model_names[[model_number]]]] <- vector(mode = "list", length = length(train_indices))
                                  for(j in 1:n_hyper_samples){
                                      hyper_samples[[model_names[[model_number]]]][[j]] <- vector(mode = "list", length = length(train_indices))
                                  }
                            }

                                
                            for(fold in 1:length(train_indices)){                                 
                                  for(model_number in 1:length(models)){
                                        if(print){
                                            cat(paste("Fold:",fold,"/",length(train_indices),"\n"))
                                            if(!is.null(model_names)){
                                              cat(paste("Model:",model_names[[model_number]],"\n"))                                                                                          
                                            } else{
                                              cat(paste("Model:", model_number,"\n"))
                                            }
                                          }

                                      # Generate posterior samples of the mean
                                      if(is.null(models[[model_number]]$.args)){
                                        stop("There was a problem with INLA's fit. Please, check your model specifications carefully and re-fit the model.")
                                      }


                                      df_train <- select_indexes(data, train_indices[[fold]])
                                      df_pred <- select_indexes(data, test_indices[[fold]])

                                      df_pred <- prepare_df_pred(df_pred, models[[model_number]], test_indices[[fold]])
                                      new_model <- bru_rerun_with_data(models[[model_number]], train_indices[[fold]], true_CV = !pseudo_predict, fit_verbose = fit_verbose)

                                      if(print){
                                          cat("Generating samples...\n")
                                      }

                                      post_samples[[model_names[[model_number]]]][[fold]] <- inlabru::generate(new_model, newdata = df_pred, formula = formula, n.samples = n_samples)

                                        if(print){ 
                                          cat("Samples generated!\n")
                                        }

                                        if(nrow(post_samples[[model_names[[model_number]]]][[fold]]) == 1){
                                          post_samples[[model_names[[model_number]]]][[fold]] <- matrix(rep(post_samples[[model_names[[model_number]]]][[fold]], length(test_indices[[fold]])),ncol=ncol(post_samples[[model_names[[model_number]]]][[fold]]), byrow = TRUE)
                                        }                                        

                                        if(compute_posterior_means){
                                          post_means[[model_names[[model_number]]]][[fold]] <- rowMeans(post_samples[[model_names[[model_number]]]][[fold]])
                                        }

                                        hyper_marginals[[model_names[[model_number]]]][[fold]] <- new_model$marginals.hyperpar

                                        hyper_summary[[model_names[[model_number]]]][[fold]] <- new_model$summary.hyperpar

                                        if(return_hyper_samples){
                                          for(j in 1:n_hyper_samples){
                                              hyper_samples[[model_names[[model_number]]]][[j]][[fold]] <- INLA::inla.hyperpar.sample(n_samples, new_model, improve.marginals=TRUE) 
                                          }
                                        }
                                  }
                            }          


              out <- list()
              if(return_samples){
                out[["post_samples"]] <- post_samples
              }

              if(return_hyper_samples){
                out[["hyper_samples"]] <- hyper_samples
              }

              out[["hyper_marginals"]] <- hyper_marginals
              out[["hyper_summary"]] <- hyper_summary
              if(compute_posterior_means){
                out[["post_means"]] <- post_means
              }

              return(out)
          }