#' @name rspde.matern1d
#' @title Matern rSPDE model object for INLA
#' @description Creates an INLA object for a stationary Matern model with
#' general smoothness parameter.
#' @param loc A vector of spatial locations.
#' @param nu.upper.bound Upper bound for the smoothness parameter. If `NULL`, it will be set to 2.
#' @param rspde.order The order of the covariance-based rational SPDE approach. The default order is 1.
#' @param nu If nu is set to a parameter, nu will be kept fixed and will not
#' be estimated. If nu is `NULL`, it will be estimated.
#' @param parameterization Which parameterization to use? `matern` uses range, std. deviation and nu (smoothness). `spde` uses kappa, tau and nu (smoothness). `matern2` uses range-like (1/kappa), variance and nu (smoothness). The default is `spde`.
#' @param prior.kappa a `list` containing the elements `meanlog` and
#' `sdlog`, that is, the mean and standard deviation on the log scale.
#' @param prior.nu a list containing the elements `mean` and `prec`
#' for beta distribution, or `loglocation` and `logscale` for a
#' truncated lognormal distribution. `loglocation` stands for
#' the location parameter of the truncated lognormal distribution in the log
#' scale. `prec` stands for the precision of a beta distribution.
#' `logscale` stands for the scale of the truncated lognormal
#' distribution on the log scale. Check details below.
#' @param prior.tau a list containing the elements `meanlog` and
#' `sdlog`, that is, the mean and standard deviation on the log scale.
#' @param prior.range a `list` containing the elements `meanlog` and
#' `sdlog`, that is, the mean and standard deviation on the log scale. Will not be used if prior.kappa is non-null.
#' @param prior.std.dev a `list` containing the elements `meanlog` and
#' `sdlog`, that is, the mean and standard deviation on the log scale. Will not be used if prior.tau is non-null.
#' @param start.lkappa Starting value for log of kappa.
#' @param start.nu Starting value for nu.
#' @param start.theta Starting values for the model parameters. In the stationary case, if `parameterization='matern'`, then `theta[1]` is the std.dev and `theta[2]` is the range parameter.
#' If `parameterization = 'spde'`, then `theta[1]` is `tau` and `theta[2]` is `kappa`.
#' @param theta.prior.mean A vector for the mean priors of `theta`.
#' @param theta.prior.prec A precision matrix for the prior of `theta`.
#' @param prior.std.dev.nominal Prior std. deviation to be used for the priors and for the starting values.
#' @param prior.range.nominal Prior range to be used for the priors and for the starting values.
#' @param prior.kappa.mean Prior kappa to be used for the priors and for the starting values.
#' @param prior.tau.mean Prior tau to be used for the priors and for the starting values.
#' @param start.lstd.dev Starting value for log of std. deviation. Will not be used if start.ltau is non-null. Will be only used in the stationary case and if `parameterization = 'matern'`.
#' @param start.lrange Starting value for log of range. Will not be used if start.lkappa is non-null. Will be only used in the stationary case and if `parameterization = 'matern'`.
#' @param start.ltau Starting value for log of tau. Will be only used in the stationary case and if `parameterization = 'spde'`.
#' @param start.lkappa Starting value for log of kappa. Will be only used in the stationary case and if `parameterization = 'spde'`.
#' @param prior.theta.param Should the lognormal prior be on `theta` or on the SPDE parameters (`tau` and `kappa` on the stationary case)?
#' @param prior.nu.dist The distribution of the smoothness parameter.
#' The current options are "beta" or "lognormal". The default is "lognormal".
#' @param nu.prec.inc Amount to increase the precision in the beta prior
#' distribution.
#' @param type.rational.approx Which type of rational approximation
#' should be used? The current types are "chebfun", "brasil" or "chebfunLB".
#' @param debug INLA debug argument
#' @param shared_lib Which shared lib to use for the cgeneric implementation?
#' If "detect", it will check if the shared lib exists locally, in which case it will
#' use it. Otherwise it will use INLA's shared library.
#' If "INLA", it will use the shared lib from INLA's installation. If 'rSPDE', then
#' it will use the local installation (does not work if your installation is from CRAN).
#' Otherwise, you can directly supply the path of the .so (or .dll) file.
#' @param ... Only being used internally.
#'
#' @return An INLA model.
#' @export
rspde.matern1d <- function(loc,
                         nu.upper.bound = NULL, 
                         rspde.order = 1,
                         nu = NULL,
                         parameterization = c("spde", "matern", "matern2"),
                         start.nu = NULL,
                         start.theta = NULL,
                         prior.nu = NULL,
                         theta.prior.mean = NULL,
                         theta.prior.prec = 0.1,
                         prior.std.dev.nominal = 1,
                         prior.range.nominal = NULL,
                         prior.kappa.mean = NULL,
                         prior.tau.mean = NULL,
                         start.lstd.dev = NULL,
                         start.lrange = NULL,
                         start.ltau = NULL,
                         start.lkappa = NULL,
                         prior.theta.param = c("theta", "spde"),
                         prior.nu.dist = c("beta", "lognormal"),
                         nu.prec.inc = 1,
                         type.rational.approx = c(
                             "chebfun",
                             "brasil", "chebfunLB"
                         ),
                         debug = FALSE,
                         shared_lib = "detect",
                         ...) {
    type.rational.approx <- type.rational.approx[[1]]
    
    prior.theta.param <- prior.theta.param[[1]]
    
    if (!(prior.theta.param %in% c("theta", "spde"))) {
        stop("theta.theta.param should be either 'theta' or 'spde'!")
    }
    
    parameterization <- parameterization[[1]]
    
    prior.nu.dist <- prior.nu.dist[[1]]
    if (!prior.nu.dist %in% c("beta", "lognormal")) {
        stop("prior.nu.dist should be either 'beta' or 'lognormal'!")
    }
    
    if (!parameterization %in% c("matern", "spde", "matern2")) {
        stop("parameterization should be either 'matern', 'spde' or 'matern2'!")
    }
    
    if (!type.rational.approx %in% c("chebfun", "brasil", "chebfunLB")) {
        stop("type.rational.approx should be either 'chebfun', 'brasil' or 'chebfunLB'!")
    }
    
    if(length(unique(diff(loc))) == 1) {
        equally_spaced <- TRUE
    } else {
        equally_spaced <- FALSE
    }
    integer.nu <- FALSE
    
    stationary <- FALSE
    
    if(is.null(nu.upper.bound)){
        nu.upper.bound <- 2
    }

    if(nu.upper.bound + 0.5 - floor(nu.upper.bound + 0.5) == 0){
        nu.upper.bound <- nu.upper.bound - 1e-5
    }
    
    fixed_nu <- !is.null(nu)
    if (fixed_nu) {
        nu_order <- nu
        start.nu <- nu
    } else {
        nu_order <- nu.upper.bound
    }
    d = 1
    beta <- nu_order / 2 + d / 4
    
    m_alpha <- floor(2 * beta)
    
    if (!is.null(nu)) {
        if (!is.numeric(nu)) {
            stop("nu must be numeric!")
        }
    }
    
    if (fixed_nu) {
        alpha <- nu + d / 2
        integer_alpha <- (alpha %% 1 == 0)
        if (!integer_alpha) {
            if (rspde.order > 0) {
                rational_table <- get_rational_coefficients(rspde.order, type.rational.approx)
            } 
        } else {
            rational_table <- get_rational_coefficients(1, type.rational.approx)
        }
    } else {
        integer_alpha <- FALSE
        if (rspde.order > 0) {
            rational_table <- get_rational_coefficients(rspde.order, type.rational.approx)
        }
    }
    
    ### Location of object files
    
    rspde_lib <- get_shared_library(shared_lib)
    
    ### PRIORS AND STARTING VALUES
    
    # Prior nu
    
    if (is.null(prior.nu$loglocation)) {
        prior.nu$loglocation <- log(min(1, nu.upper.bound / 2))
    }
    
    if (is.null(prior.nu[["mean"]])) {
        prior.nu[["mean"]] <- min(1, nu.upper.bound / 2)
    }
    
    if (is.null(prior.nu$prec)) {
        mu_temp <- prior.nu[["mean"]] / nu.upper.bound
        prior.nu$prec <- max(1 / mu_temp, 1 / (1 - mu_temp)) + nu.prec.inc
    }
    
    if (is.null(prior.nu[["logscale"]])) {
        prior.nu[["logscale"]] <- 1
    }
    
    # Start nu
    
    if (is.null(start.nu)) {
        if (prior.nu.dist == "beta") {
            start.nu <- prior.nu[["mean"]]
        } else if (prior.nu.dist == "lognormal") {
            start.nu <- exp(prior.nu[["loglocation"]])
        } else {
            stop("prior.nu.dist should be either beta or lognormal!")
        }
    } else if (start.nu > nu.upper.bound || start.nu < 0) {
        stop("start.nu should be a number between 0 and nu.upper.bound!")
    }
    
    
    # Prior kappa and prior range
    param <- get_parameters_rSPDE(
        NULL, alpha,
        matrix(c(0, 1, 0), 1, 3),
        matrix(c(0, 0, 1), 1, 3),
        matrix(c(0, 1, 0), 1, 3),
        matrix(c(0, 0, 1), 1, 3),
        start.nu,
        start.nu + 1 / 2,
        parameterization,
        prior.std.dev.nominal,
        prior.range.nominal,
        prior.tau.mean,
        prior.kappa.mean,
        theta.prior.mean,
        theta.prior.prec,
        mesh.range = diff(range(loc)),
        d = 1,
        n.spde = length(loc)
    )
    
    if (is.null(start.theta)) {
        start.theta <- param$theta.prior.mean
    }
    
    theta.prior.mean <- param$theta.prior.mean
    theta.prior.prec <- param$theta.prior.prec
    

    # Starting values
    if (parameterization == "spde") {
        if (!is.null(start.lkappa)) {
            start.theta[2] <- start.lkappa
        }
        if (!is.null(start.ltau)) {
            start.theta[1] <- start.ltau
        }
    } else if (parameterization == "matern") {
        if (!is.null(start.lrange)) {
            start.theta[2] <- start.lrange
        }
        if (!is.null(start.lstd.dev)) {
            start.theta[1] <- start.lstd.dev
        }
    } else if (parameterization == "matern2") {
        if (!is.null(start.lrange)) {
            start.theta[2] <- start.lrange
        }
        if (!is.null(start.lstd.dev)) {
            start.theta[1] <- 2 * start.lstd.dev
        }
    }

    
    if (!fixed_nu) {
        tmp <- matern.rational.precision(loc = loc,
                                         order = rspde.order,
                                         nu = nu.upper.bound,
                                         kappa = 1,
                                         sigma = 1)
        graph_opt <- tmp$Q
        A <- tmp$A
        n_cgeneric <- dim(graph_opt)[1]        
                
        graph_opt <- transpose_cgeneric(graph_opt)
                
        model <- do.call(
                    eval(parse(text = "INLA::inla.cgeneric.define")),
                    list(
                        model = "inla_cgeneric_rspde_1d_general_model",
                        shlib = rspde_lib,
                        n = as.integer(n_cgeneric), 
                        debug = debug,
                        nu_upper_bound = nu.upper.bound,
                        rational_table = as.matrix(rational_table),
                        graph_opt_i = graph_opt@i,
                        graph_opt_j = graph_opt@j,
                        start_theta = start.theta,
                        theta_prior_mean = theta.prior.mean,
                        theta_prior_prec = theta.prior.prec,
                        prior_nu_loglocation = prior.nu$loglocation,
                        prior_nu_mean = prior.nu$mean,
                        prior_nu_prec = prior.nu$prec,
                        prior_nu_logscale = prior.nu$logscale,
                        start_nu = start.nu,
                        rspde_order = as.integer(rspde.order),
                        prior_nu_dist = prior.nu.dist,
                        parameterization = parameterization,
                        prior_theta_param = prior.theta.param,
                        loc = loc,
                        es = as.integer(equally_spaced),
                        nu_fixed = as.integer(0)
                    )
                )
        
            
    model$cgeneric_type <- "general"
    } else if (!integer_alpha) {
        tmp <- matern.rational.precision(loc = loc,
                                         order = rspde.order,
                                         nu = nu,
                                         kappa = 1,
                                         sigma = 1)
        graph_opt <- tmp$Q
        A <- tmp$A
        n_cgeneric <- dim(graph_opt)[1]
        
        graph_opt <- transpose_cgeneric(graph_opt)
        
        model <- do.call(
            eval(parse(text = "INLA::inla.cgeneric.define")),
            list(
                model = "inla_cgeneric_rspde_1d_general_model",
                shlib = rspde_lib,
                n = as.integer(n_cgeneric), 
                debug = debug,
                nu_upper_bound = nu,
                rational_table = as.matrix(rational_table),
                graph_opt_i = graph_opt@i,
                graph_opt_j = graph_opt@j,
                start_theta = start.theta,
                theta_prior_mean = theta.prior.mean,
                theta_prior_prec = theta.prior.prec,
                prior_nu_loglocation = prior.nu$loglocation,
                prior_nu_mean = prior.nu$mean,
                prior_nu_prec = prior.nu$prec,
                prior_nu_logscale = prior.nu$logscale,
                start_nu = nu,
                rspde_order = as.integer(rspde.order),
                prior_nu_dist = prior.nu.dist,
                parameterization = parameterization,
                prior_theta_param = prior.theta.param,
                loc = loc,
                es = as.integer(equally_spaced),
                nu_fixed = as.integer(1)
            )
        )
        
            
    
            
        model$cgeneric_type <- "frac_alpha"
    } else {
       
        tmp <- matern.rational.precision(loc = loc,
                                         order = rspde.order,
                                         nu = nu,
                                         kappa = 1,
                                         sigma = 1)
        graph_opt <- tmp$Q
        A <- tmp$A
        n_cgeneric <- dim(graph_opt)[1]
        
        graph_opt <- transpose_cgeneric(graph_opt)
        
        model <- do.call(
            eval(parse(text = "INLA::inla.cgeneric.define")),
            list(
                model = "inla_cgeneric_rspde_1d_general_model",
                shlib = rspde_lib,
                n = as.integer(n_cgeneric), 
                debug = debug,
                nu_upper_bound = nu,
                rational_table = as.matrix(rational_table),
                graph_opt_i = graph_opt@i,
                graph_opt_j = graph_opt@j,
                start_theta = start.theta,
                theta_prior_mean = theta.prior.mean,
                theta_prior_prec = theta.prior.prec,
                prior_nu_loglocation = prior.nu$loglocation,
                prior_nu_mean = prior.nu$mean,
                prior_nu_prec = prior.nu$prec,
                prior_nu_logscale = prior.nu$logscale,
                start_nu = nu,
                rspde_order = as.integer(0),
                prior_nu_dist = prior.nu.dist,
                parameterization = parameterization,
                prior_theta_param = prior.theta.param,
                loc = loc,
                es = as.integer(equally_spaced),
                nu_fixed = as.integer(1)
            )
        )
        model$cgeneric_type <- "int_alpha"
    }
        
    model$nu <- nu
    model$theta.prior.mean <- theta.prior.mean
    model$prior.nu <- prior.nu
    model$theta.prior.prec <- theta.prior.prec
    model$start.nu <- start.nu
    model$integer.nu <- ifelse(fixed_nu, integer_alpha, FALSE)
    model$start.theta <- start.theta
    if (integer.nu) {
        rspde.order <- 0
    }
    model$rspde.order <- rspde.order

    rspde_check_cgeneric_symbol(model)
    
    class(model) <- c("inla_rspde_matern1d", class(model))
    model$dim <- d
    model$est_nu <- !fixed_nu
    model$nu.upper.bound <- nu.upper.bound
    model$prior.nu.dist <- prior.nu.dist
    model$debug <- debug
    model$type.rational.approx <- type.rational.approx
    model$parameterization <- parameterization
    model$A <- A
    model$index <- INLA::inla.spde.make.index(name = "field", n.spde = dim(A)[2])
    model$rspde_version <- as.character(packageVersion("rSPDE"))
    model$stationary = TRUE
    model$loc <- loc
    return(model)
}



#' @title rSPDE stationary inlabru mapper
#' @name bru_get_mapper.inla_rspde_matern1d
#' @param model An `inla_rspde_matern1d` object for which to construct or extract a mapper
#' @param \dots Arguments passed on to other methods
#' @rdname bru_get_mapper.inla_rspde_matern1d
#' @rawNamespace if (getRversion() >= "3.6.0") {
#'   S3method(inlabru::bru_get_mapper, inla_rspde_matern1d)
#'   S3method(inlabru::ibm_n, bru_mapper_inla_rspde_matern1d)
#'   S3method(inlabru::ibm_values, bru_mapper_inla_rspde_matern1d)
#'   S3method(inlabru::ibm_jacobian, bru_mapper_inla_rspde_matern1d)
#' }

bru_get_mapper.inla_rspde_matern1d <- function(model, ...) {
  stopifnot(requireNamespace("inlabru"))
  mapper <- list(model = model)
  inlabru::bru_mapper_define(mapper, new_class = "bru_mapper_inla_rspde_matern1d")
}


#' @param mapper A `bru_mapper_inla_rspde_matern1d` object
#' @rdname bru_get_mapper.inla_rspde_matern1d
ibm_n.bru_mapper_inla_rspde_matern1d <- function(mapper, ...) {
  model <- mapper[["model"]]
  return(model$f$n)
}
#' @rdname bru_get_mapper.inla_rspde_matern1d
ibm_values.bru_mapper_inla_rspde_matern1d <- function(mapper, ...) {
  seq_len(inlabru::ibm_n(mapper))
}
#' @param input The values for which to produce a mapping matrix
#' @rdname bru_get_mapper.inla_rspde_matern1d
ibm_jacobian.bru_mapper_inla_rspde_matern1d <- function(mapper, input, ...) {
  model <- mapper[["model"]]
  loc <- model[["loc"]]
  indices <- match_with_tolerance(input, loc)
  A <- model[["A"]]
  return(A[indices,, drop=FALSE])
}

#' @name predict.inla_rspde_matern1d
#' @title Predict method for 'inlabru' stationary Matern 1d models
#' @description Auxiliar function to obtain predictions of the stationary Matern 1d models
#' using 'inlabru'.
#' @param object An `inla_rspde_matern1d` object built with the `rspde.matern1d()`
#' function.
#' @param cmp The 'inlabru' component used to fit the model.
#' @param bru_fit A fitted model using 'inlabru' or 'INLA'.
#' @param newdata A data.frame of covariates needed for the prediction. 
#' @param formula A formula where the right hand side defines an R expression to
#' evaluate for each generated sample. If NULL, the latent and hyperparameter
#' states are returned as named list elements. See Details for more information.
#' @param n.samples Integer setting the number of samples to draw in order to
#' calculate the posterior statistics. The default is rather low but provides a
#' quick approximate result.
#' @param seed Random number generator seed passed on to `inla.posterior.sample()`
#' @param probs	A numeric vector of probabilities with values in the standard
#' unit interval to be passed to stats::quantile
#' @param return_original_order Should the predictions be returned in the
#' original order?
#' @param num.threads	Specification of desired number of threads for parallel
#' computations. Default NULL, leaves it up to 'INLA'. When seed != 0, overridden to "1:1"
#' @param include	Character vector of component labels that are needed by the
#' predictor expression; Default: NULL (include all components that are not
#' explicitly excluded)
#' @param exclude	Character vector of component labels that are not used by the
#' predictor expression. The exclusion list is applied to the list as determined
#' by the include parameter; Default: NULL (do not remove any components from
#' the inclusion list)
#' @param drop logical; If keep=FALSE, data is a SpatialDataFrame, and the
#' prediciton summary has the same number of rows as data, then the output is a
#' SpatialDataFrame object. Default FALSE.
#' @param tolerance Tolerance for merging locations.
#' @param... Additional arguments passed on to `inla.posterior.sample()`.
#' @return A list with predictions.
#' @export

predict.inla_rspde_matern1d <- function(object,
                                           cmp,
                                           bru_fit,
                                           newdata = NULL,
                                           formula = NULL,
                                           n.samples = 100,
                                           seed = 0L,
                                           probs = c(0.025, 0.5, 0.975),
                                           return_original_order = TRUE,
                                           num.threads = NULL,
                                           include = NULL,
                                           exclude = NULL,
                                           drop = FALSE,
                                           tolerance = 1e-4,
                                           ...){
  if(length(bru_fit$bru_info$lhoods) > 1){
    stop("Only models with one likelihood implemented.")
  }

  name_locations <- bru_fit$bru_info$model$effects$field$main$input$input
  
  original_data <- bru_fit$bru_info$lhoods[[1]]$data

  new_data <- newdata
  new_data[["__new"]] <- TRUE
  n_locations <- nrow(newdata[[name_locations]])
  names_columns <- names(original_data)

  new_data <- merge_with_tolerance(original_data, new_data, by = as.character(name_locations), tolerance = tolerance)

  spde____model <- rspde.matern1d(loc = new_data[[name_locations]], 
                                    rspde.order = object[["rspde.order"]],
                                    nu.upper.bound = object[["nu.upper.bound"]],
                                    nu = object[["nu"]],
                                    parameterization = object[["parameterization"]],
                                    start.nu = object[["start.nu"]],
                                    start.theta = object[["start.theta"]],
                                    prior.nu = object[["prior.nu"]],
                                    theta.prior.mean = object[["theta.prior.mean"]],
                                    theta.prior.prec = object[["theta.prior.prec"]],
                                    prior.nu.dist = object[["prior.nu.dist"]],
                                    type.rational.approx = object[["type.rational.approx"]])

  cmp_c <- as.character(cmp)
  name_model <- deparse(substitute(object))
  cmp_c[3] <- sub(name_model, "spde____model", cmp_c[3])
  cmp_new <- as.formula(paste(cmp_c[2], cmp_c[1], cmp_c[3]))

  info <- bru_fit[["bru_info"]]
  info[["options"]] <- inlabru::bru_call_options(inlabru::bru_options(info[["options"]]))

  bru_fit_new <- inlabru::bru(cmp_new,
          data = new_data, options = info[["options"]])
  
  pred <- predict(object = bru_fit_new,
                    newdata = newdata,
                    formula = formula,
                    n.samples = n.samples,
                    seed = seed,
                    probs = probs,
                    num.threads = num.threads,
                    include = include,
                    exclude = exclude,
                    drop = drop,
                    ...)

  return(pred)
}
