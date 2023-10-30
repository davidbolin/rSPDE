
#' @name rspde.matern
#' @title Matern rSPDE model object for INLA
#' @description Creates an INLA object for a stationary Matern model with
#' general smoothness parameter.
#' @param mesh The mesh to build the model. It can be an `inla.mesh` or
#' an `inla.mesh.1d` object. Otherwise, should be a list containing elements d, the dimension, C, the mass matrix,
#' and G, the stiffness matrix.
#' @param nu.upper.bound Upper bound for the smoothness parameter.
#' @param rspde.order The order of the covariance-based rational SPDE approach.
#' @param nu If nu is set to a parameter, nu will be kept fixed and will not
#' be estimated. If nu is `NULL`, it will be estimated.
#' @param B.sigma Matrix with specification of log-linear model for \eqn{\sigma} (for 'matern' parameterization) or for \eqn{\sigma^2} (for 'matern2' parameterization). Will be used if `parameterization = 'matern'` or `parameterization = 'matern2'`. 
#' @param B.range Matrix with specification of log-linear model for \eqn{\rho}, which is a range-like parameter (it is exactly the range parameter in the stationary case). Will be used if `parameterization = 'matern'` or `parameterization = 'matern2'`.
#' @param parameterization Which parameterization to use? `matern` uses range, std. deviation and nu (smoothness). `spde` uses kappa, tau and nu (smoothness). `matern2` uses range-like (1/kappa), variance and nu (smoothness). The default is `spde`.
#' @param B.tau Matrix with specification of log-linear model for \eqn{\tau}. Will be used if `parameterization = 'spde'`.
#' @param B.kappa Matrix with specification of log-linear model for \eqn{\kappa}. Will be used if `parameterization = 'spde'`.
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
#' distribution. Check details below.
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

rspde.matern <- function(mesh,
                         nu.upper.bound = 4, rspde.order = 2,
                         nu = NULL, 
                         B.sigma = matrix(c(0, 1, 0), 1, 3), 
                         B.range = matrix(c(0, 0, 1), 1, 3), 
                         parameterization = c("spde", "matern", "matern2"),
                         B.tau = matrix(c(0, 1, 0), 1, 3), 
                         B.kappa = matrix(c(0, 0, 1), 1, 3), 
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
                         type.rational.approx = c("chebfun",
                         "brasil", "chebfunLB"),
                         debug = FALSE,
                         shared_lib = "detect",
                         ...) {
  type.rational.approx <- type.rational.approx[[1]]

  prior.theta.param <- prior.theta.param[[1]]

  if(!(prior.theta.param %in% c("theta", "spde"))){
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

  if (parameterization == "spde"){
    if(!missing(B.range)){
      warning("B.range was passed, but will not be used since the parameterization is 'spde'.")
    }
    if(!missing(B.sigma)){
      warning("B.sigma was passed, but will not be used since the parameterization is 'spde'.")
    }    

    if(ncol(B.kappa)!=ncol(B.tau)){
      stop("B.kappa and B.tau must have the same number of columns.")
    }

    tmp_B <- rbind(B.kappa, B.tau)
    tmp_B <- colSums(tmp_B^2)
    tmp_B <- tmp_B[-1]
    if(any(tmp_B == 0)){
      stop("The only column that is allowed to be zero simultaneously on B.kappa and B.tau is the first column.")
    }
  }

  if(parameterization %in% c("matern", "matern2")){
    if(!missing(B.kappa)){
      warning("B.kappa was passed, but will not be used since the parameterization is NOT 'spde'.")      
    }
    if(!missing(B.tau)){
      warning("B.tau was passed, but will not be used since the parameterization is NOT 'spde'.")      
    }

    if(ncol(B.sigma)!=ncol(B.range)){
      stop("B.sigma and B.range must have the same number of columns.")
    }

    tmp_B <- rbind(B.sigma, B.range)
    tmp_B <- colSums(tmp_B^2)
    tmp_B <- tmp_B[-1]
    if(any(tmp_B == 0)){
      stop("The only column that is allowed to be zero simultaneously on B.sigma and B.range is the first column.")
    }
  }

  integer.nu <- FALSE

  stationary <- FALSE

  if(parameterization == "spde"){
    if(nrow(B.tau) == 1 && nrow(B.kappa) == 1){
      if(B.tau[1,1] == 0 && B.tau[1,2] == 1 && B.tau[1,3] == 0 &&
          B.kappa[1,1] == 0 && B.kappa[1,2] == 0 && B.kappa[1,3] == 1){
            stationary <- TRUE
          }
    }
  } else{
    if(nrow(B.sigma) == 1 && nrow(B.range) == 1){
      if(B.sigma[1,1] == 0 && B.sigma[1,2] == 1 && B.sigma[1,3] == 0 &&
          B.range[1,1] == 0 && B.range[1,2] == 0 && B.range[1,3] == 1){
            stationary <- TRUE
          }
    }
  }

  if(inherits(mesh, c("inla.mesh", "inla.mesh.1d"))){
    d <- get_inla_mesh_dimension(mesh)
  } else if (!is.null(mesh$d)){
    d <- mesh$d
  } else{
    stop("The mesh object should either be an INLA mesh object or contain d, the dimension!")
  }

  if (nu.upper.bound - floor(nu.upper.bound) == 0) {
    nu.upper.bound <- nu.upper.bound - 1e-5
  }
  fixed_nu <- !is.null(nu)
  if (fixed_nu) {
    nu_order <- nu
    start.nu <- nu
  } else {
    nu_order <- nu.upper.bound
  }

  beta <- nu_order / 2 + d / 4

  m_alpha <- floor(2 * beta)

  if (!is.null(nu)) {
    if (!is.numeric(nu)) {
      stop("nu must be numeric!")
    }
  }

  if (d == 1) {
    if (nu.upper.bound > 2) {
      warning("In dimension 1 you can have unstable results
      for nu.upper.bound > 2. Consider changing
      nu.upper.bound to 2 or 1.")
    }
  }


  if (fixed_nu) {
    alpha <- nu + d / 2
    integer_alpha <- (alpha %% 1 == 0)
    if(!integer_alpha){
      if(rspde.order > 0){
      n_m <- rspde.order
            mt <- get_rational_coefficients(rspde.order, type.rational.approx)
            r <- sapply(1:(n_m), function(i) {
            approx(mt$alpha, mt[[paste0("r", i)]], cut_decimals(2 * beta))$y
            })
            p <- sapply(1:(n_m), function(i) {
            approx(mt$alpha, mt[[paste0("p", i)]], cut_decimals(2 * beta))$y
            })
            k <- approx(mt$alpha, mt$k, cut_decimals(2 * beta))$y
      }
    }
  } else {
    integer_alpha <- FALSE
    if(rspde.order > 0){
      rational_table <- get_rational_coefficients(rspde.order, type.rational.approx)
    }
  }

  ### Location of object files

  rspde_lib <- shared_lib

  if(shared_lib == "INLA"){
    rspde_lib <- INLA::inla.external.lib('rSPDE')
  } else if(shared_lib == "rSPDE"){
    rspde_lib <- system.file('shared', package='rSPDE')
    if(Sys.info()['sysname']=='Windows') {
		rspde_lib <- paste0(rspde_lib, "/rspde_cgeneric_models.dll")
            } else {
		rspde_lib <- paste0(rspde_lib, "/rspde_cgeneric_models.so")
            }
  } else if(shared_lib == "detect"){
    rspde_lib_local <- system.file('shared', package='rSPDE')
    if(Sys.info()['sysname']=='Windows') {
		rspde_lib_local <- paste0(rspde_lib_local, "/rspde_cgeneric_models.dll")
            } else {
		rspde_lib_local <- paste0(rspde_lib_local, "/rspde_cgeneric_models.so")
            }
    if(file.exists(rspde_lib_local)){
      rspde_lib <- rspde_lib_local
    } else{
      rspde_lib <- INLA::inla.external.lib('rSPDE')
    }
  }

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

  if(!inherits(mesh, "metric_graph")){
    param <- get_parameters_rSPDE(mesh, 2 * beta, 
                      B.tau, 
                      B.kappa, 
                      B.sigma,
                      B.range, 
                      start.nu,
                      start.nu + d/2,
                      parameterization,
                      prior.std.dev.nominal, 
                      prior.range.nominal, 
                      prior.tau.mean, 
                      prior.kappa.mean, 
                      theta.prior.mean, 
                      theta.prior.prec) 
  } else{
    tmp_function <- function(vec_param){
      vec_param
    }
    param <- tmp_function(...)
  }



    if(is.null(start.theta)){
      start.theta <- param$theta.prior.mean
    }

    theta.prior.mean <- param$theta.prior.mean
    theta.prior.prec <- param$theta.prior.prec

    B.tau <- param$B.tau
    B.kappa <- param$B.kappa

  # Starting values
  if(stationary){
      if(parameterization == "spde"){
          if (!is.null(start.lkappa)) {
            start.theta[2] <- start.lkappa
          }
          if (!is.null(start.ltau)) {
            start.theta[1] <- start.ltau
          }
      } else if(parameterization == "matern"){
          if(!is.null(start.lrange)){
            start.theta[2] <- start.lrange
          }
          if(!is.null(start.lstd.dev)){
            start.theta[1] <- start.lstd.dev
          }
      } else if(parameterization == "matern2"){
          if(!is.null(start.lrange)){
            start.theta[2] <- start.lrange
          }
          if(!is.null(start.lstd.dev)){
            start.theta[1] <- 2*start.lstd.dev
          }
      }
  }



  ### STATIONARY PART
  if(stationary){

    if(inherits(mesh, c("inla.mesh", "inla.mesh.1d"))){
      if (integer_alpha) {
        integer.nu <- TRUE
        if (d == 1) {
          fem_mesh <- fem_mesh_order_1d(mesh, m_order = m_alpha + 1)
        } else {
          # fem_mesh <- INLA::inla.mesh.fem(mesh, order = m_alpha)
          # fem_mesh <- fmesher::fm_fem(mesh, order = m_alpha)
          fem_mesh <- fm_fem(mesh, order = m_alpha)
        }
      } else {
        if (d == 1) {
          fem_mesh <- fem_mesh_order_1d(mesh, m_order = m_alpha + 2)
        } else {
          # fem_mesh <- INLA::inla.mesh.fem(mesh, order = m_alpha + 1)
          # fem_mesh <- fmesher::fm_fem(mesh, order = m_alpha + 1)
          fem_mesh <- fm_fem(mesh, order = m_alpha + 1)
        }
      }
    } else{
      if(is.null(mesh$C) || is.null(mesh$G)){
        stop("If mesh is not an inla.mesh object, you should manually supply a list with elements c0, g1, g2...")
      }
      fem_mesh <- generic_fem_mesh_order(mesh, m_order = m_alpha + 2)
    }


    n_cgeneric <- ncol(fem_mesh[["c0"]])

    fem_mesh_orig <- fem_mesh

    fem_mesh <- fem_mesh[setdiff(names(fem_mesh),c("ta","va"))]

    fem_mesh <- lapply(fem_mesh, transpose_cgeneric)

    if (integer_alpha) {
      result_sparsity <- analyze_sparsity_rspde(
        nu.upper.bound = nu_order, dim = d,
        rspde.order = rspde.order,
        fem_mesh_matrices = fem_mesh,
        include_higher_order = FALSE
      )
    } else if(rspde.order > 0) {
        result_sparsity <- analyze_sparsity_rspde(
          nu.upper.bound = nu_order, dim = d,
          rspde.order = rspde.order,
          fem_mesh_matrices = fem_mesh
        )
      positions_matrices <- result_sparsity$positions_matrices
    } else{
            result_sparsity <- analyze_sparsity_rspde(
          nu.upper.bound = nu_order, dim = d,
          rspde.order = rspde.order,
          fem_mesh_matrices = fem_mesh, 
          include_lower_order = FALSE
        )
      positions_matrices <- result_sparsity$positions_matrices
    }

    idx_matrices <- result_sparsity$idx_matrices

    if(rspde.order > 0 || integer_alpha){
      positions_matrices_less <- result_sparsity$positions_matrices_less
    }

    # if (integer_alpha) {
    if(rspde.order > 0 || integer_alpha){  
    n_tmp <- length(
        fem_mesh[[paste0("g", m_alpha)]]@x[idx_matrices[[m_alpha + 1]]]
        )

    tmp <-
    rep(0, n_tmp)
    tmp[positions_matrices_less[[1]]] <-
    fem_mesh$c0@x[idx_matrices[[1]]]
    matrices_less <- tmp

    if(m_alpha > 1){
      tmp <-
      rep(0, n_tmp)
      tmp[positions_matrices_less[[2]]] <-
      fem_mesh$g1@x[idx_matrices[[2]]]
      matrices_less <- c(matrices_less, tmp)
    }



    if (m_alpha > 2) {
    for (j in 2:(m_alpha - 1)) {
        tmp <-
        rep(0, n_tmp)
        tmp[positions_matrices_less[[j + 1]]] <-
        fem_mesh[[paste0("g", j)]]@x[idx_matrices[[j + 1]]]
        matrices_less <- c(matrices_less, tmp)
    }
    }

    tmp <- fem_mesh[[paste0("g", m_alpha)]]@x[idx_matrices[[m_alpha + 1]]]

    matrices_less <- c(matrices_less, tmp)
    }
    # }


    if (!integer_alpha) {
      if (m_alpha == 0) {
        n_tmp <- length(fem_mesh$g1@x)
        tmp <-
        rep(0, n_tmp)
        tmp[positions_matrices[[1]]] <-
        fem_mesh$c0@x[idx_matrices[[1]]]
        matrices_full <- tmp
        matrices_full <- c(matrices_full, fem_mesh$g1@x)
      } else if (m_alpha > 0) {

        n_tmp <- length(
            fem_mesh[[paste0("g", m_alpha + 1)]]@x[idx_matrices[[m_alpha + 2]]]
        )

        tmp <- 
        rep(0, n_tmp)
        tmp[positions_matrices[[1]]] <-
        fem_mesh$c0@x[idx_matrices[[1]]]
        
        matrices_full <- tmp
        
        tmp <-
        rep(0, n_tmp)
        tmp[positions_matrices[[2]]] <-
        fem_mesh$g1@x[idx_matrices[[2]]]

        matrices_full <- c(matrices_full, tmp)

      if (m_alpha > 1) {
        for (j in 2:(m_alpha)) {
          tmp <-
          rep(0, n_tmp)
          tmp[positions_matrices[[j + 1]]] <-
          fem_mesh[[paste0("g", j)]]@x[idx_matrices[[j + 1]]]
          matrices_full <- c(matrices_full, tmp)
        }
      }

        tmp <-
        fem_mesh[[paste0("g", m_alpha + 1)]]@x[idx_matrices[[m_alpha + 2]]]
        matrices_full <- c(matrices_full, tmp)
      }
    }


  if (!fixed_nu) {

    if(rspde.order == 0){

      # fem_mesh has already been transposed
        graph_opt <-  fem_mesh[[paste0("g", m_alpha + 1)]]

  model <- do.call(eval(parse(text='INLA::inla.cgeneric.define')),
        list(model="inla_cgeneric_rspde_stat_parsim_gen_model",
            shlib=rspde_lib,
            n=as.integer(n_cgeneric), debug=debug,
            d = as.double(d),
            nu.upper.bound = nu.upper.bound,
            matrices_full = as.double(matrices_full),
            graph_opt_i = graph_opt@i,
            graph_opt_j = graph_opt@j,
            theta.prior.mean = theta.prior.mean,
            theta.prior.prec = theta.prior.prec,
            prior.nu.loglocation = prior.nu$loglocation,
            prior.nu.mean = prior.nu$mean,
            prior.nu.prec = prior.nu$prec,
            prior.nu.logscale = prior.nu$logscale,
            start.theta = start.theta,
            start.nu = start.nu,
            prior.nu.dist = prior.nu.dist,
            parameterization = parameterization,
            prior.theta.param = prior.theta.param))


    } else{
          graph_opt <- get.sparsity.graph.rspde(
        fem_mesh_matrices = fem_mesh_orig, dim = d,
        nu = nu.upper.bound,
        rspde.order = rspde.order,
        force_non_integer = TRUE
      )


    graph_opt <- transpose_cgeneric(graph_opt) 


    # matrices_less <- restructure_matrices_less(matrices_less, m_alpha)
    # matrices_full <- restructure_matrices_full(matrices_full, m_alpha)

    model <- do.call(eval(parse(text='INLA::inla.cgeneric.define')),
        list(model="inla_cgeneric_rspde_stat_general_model",
            shlib=rspde_lib,
            n=as.integer(n_cgeneric)*(rspde.order+1), debug=debug,
            d = as.double(d),
            nu.upper.bound = nu.upper.bound,
            matrices_less = as.double(matrices_less),
            matrices_full = as.double(matrices_full),
            rational_table = as.matrix(rational_table),
            graph_opt_i = graph_opt@i,
            graph_opt_j = graph_opt@j,
            start.theta = start.theta,
            theta.prior.mean = theta.prior.mean,
            theta.prior.prec = theta.prior.prec,
            prior.nu.loglocation = prior.nu$loglocation,
            prior.nu.mean = prior.nu$mean,
            prior.nu.prec = prior.nu$prec,
            prior.nu.logscale = prior.nu$logscale,
            start.nu = start.nu,
            rspde.order = as.integer(rspde.order),
            prior.nu.dist = prior.nu.dist,
            parameterization = parameterization,
            prior.theta.param = prior.theta.param))
    }
    
    model$cgeneric_type <- "general"
  } else if (!integer_alpha) {

      if(rspde.order == 0){
        graph_opt <- fem_mesh[[paste0("g", m_alpha + 1)]]

        model <- do.call(eval(parse(text='INLA::inla.cgeneric.define')),
        list(model="inla_cgeneric_rspde_stat_parsim_fixed_model",
            shlib=rspde_lib,
            n=as.integer(n_cgeneric), debug=debug,
            d = as.double(d),
            nu = nu,
            matrices_full = as.double(matrices_full),
            graph_opt_i = graph_opt@i,
            graph_opt_j = graph_opt@j,
            theta.prior.mean = theta.prior.mean,
            theta.prior.prec = theta.prior.prec,
            start.theta = start.theta,
            parameterization = parameterization,
            prior.theta.param = prior.theta.param))

      } else{

      graph_opt <- get.sparsity.graph.rspde(
        fem_mesh_matrices = fem_mesh_orig, dim = d,
        nu = nu,
        rspde.order = rspde.order,
        force_non_integer = TRUE
      )


    graph_opt <- transpose_cgeneric(graph_opt) 

    # matrices_less <- restructure_matrices_less(matrices_less, m_alpha)
    # matrices_full <- restructure_matrices_full(matrices_full, m_alpha)

    model <- do.call(eval(parse(text='INLA::inla.cgeneric.define')),
        list(model="inla_cgeneric_rspde_stat_frac_model",
            shlib=rspde_lib,
            n=as.integer(n_cgeneric)*(rspde.order+1), debug=debug,
            nu = nu,
            matrices_less = as.double(matrices_less),
            matrices_full = as.double(matrices_full),
            r_ratapprox = as.vector(r),
            p_ratapprox = as.vector(p),
            k_ratapprox = k,
            graph_opt_i = graph_opt@i,
            graph_opt_j = graph_opt@j,
            theta.prior.mean = theta.prior.mean,
            theta.prior.prec = theta.prior.prec,
            start.theta = start.theta,
            rspde.order = as.integer(rspde.order),
            parameterization = parameterization,
            d = as.integer(d),
            prior.theta.param = prior.theta.param))
      }
    
    model$cgeneric_type <- "frac_alpha"
  } else {
      graph_opt <- get.sparsity.graph.rspde(
        fem_mesh_matrices = fem_mesh_orig, dim = d,
        nu = nu,
        rspde.order = rspde.order
      )
      graph_opt <- transpose_cgeneric(graph_opt) 

    # matrices_less <- restructure_matrices_less(matrices_less, m_alpha)

    model <- do.call(eval(parse(text='INLA::inla.cgeneric.define')),
        list(model="inla_cgeneric_rspde_stat_int_model",
            shlib=rspde_lib,
            n=as.integer(n_cgeneric), debug=debug,
            matrices_less = as.double(matrices_less),
            m_alpha = as.integer(m_alpha),
            graph_opt_i = graph_opt@i,
            graph_opt_j = graph_opt@j,
            theta.prior.mean = theta.prior.mean,
            theta.prior.prec = theta.prior.prec,
            start.theta = start.theta,
            nu = nu,
            parameterization = parameterization,
            prior.theta.param = prior.theta.param
            ))
    model$cgeneric_type <- "int_alpha"
  }

### END OF STATIONARY PART
  } else {
### NONSTATIONARY PART

    if(inherits(mesh, c("inla.mesh", "inla.mesh.1d"))){
      if (integer_alpha) {
        integer.nu <- TRUE
        if (d == 1) {
          fem_mesh <- fem_mesh_order_1d(mesh, m_order = m_alpha + 1)
        } else {
          # fem_mesh <- INLA::inla.mesh.fem(mesh, order = m_alpha)
          # fem_mesh <- fmesher::fm_fem(mesh, order = m_alpha)
          fem_mesh <- fm_fem(mesh, order = m_alpha)
        }
      } else {
        if (d == 1) {
          fem_mesh <- fem_mesh_order_1d(mesh, m_order = m_alpha + 2)
        } else {
          # fem_mesh <- INLA::inla.mesh.fem(mesh, order = m_alpha + 1)
          # fem_mesh <- fmesher::fm_fem(mesh, order = m_alpha + 1)
          fem_mesh <- fm_fem(mesh, order = m_alpha + 1)
        }
      }
    } else{
      if(is.null(mesh$C) || is.null(mesh$G)){
        stop("If mesh is not an inla.mesh object, you should manually supply a list with elements c0, g1, g2...")
      }
      fem_mesh <- generic_fem_mesh_order(mesh, m_order = m_alpha + 2)
    }

    fem_mesh_orig <- fem_mesh

    C <- fem_mesh[["c0"]]
    G <- fem_mesh[["g1"]]

    n_cgeneric <- ncol(fem_mesh[["c0"]])

    if (!fixed_nu) {

    if(rspde.order == 0){
      warning("The order cannot be zero for nonstationary models! The order was changed to 1.")
      rspde.order <- 1
    }


        graph_opt <- get.sparsity.graph.rspde(
        fem_mesh_matrices = fem_mesh, dim = d,
        nu = nu.upper.bound,
        rspde.order = rspde.order,
        force_non_integer = TRUE)

        matern_par_tmp <- as.integer(!(parameterization == "spde"))
        matern_par_tmp <- matern_par_tmp + as.integer(parameterization == "matern2")


    graph_opt <- transpose_cgeneric(graph_opt) 

        model <- do.call(eval(parse(text='INLA::inla.cgeneric.define')),
        list(model="inla_cgeneric_rspde_nonstat_general_model",
            shlib=rspde_lib,
            n=as.integer(n_cgeneric)*(rspde.order+1), debug=debug,
            d = as.double(d),
            nu_upper_bound = nu.upper.bound,
            rational_table = as.matrix(rational_table),
            graph_opt_i = graph_opt@i,
            graph_opt_j = graph_opt@j,
            C = C,
            G = G,
            B_tau = B.tau,
            B_kappa = B.kappa,
            prior.nu.loglocation = prior.nu$loglocation,
            prior.nu.logscale = prior.nu$logscale,
            prior.nu.mean = prior.nu$mean,
            prior.nu.prec = prior.nu$prec,
            start.nu = start.nu,
            rspde_order = as.integer(rspde.order),
            prior.nu.dist = prior.nu.dist,
            start.theta = start.theta,
            theta.prior.mean = param$theta.prior.mean,
            theta.prior.prec = param$theta.prior.prec,
            matern_par = matern_par_tmp,
            prior.theta.param = prior.theta.param
            ))
    
    model$cgeneric_type <- "general"
    } else if (!integer_alpha) {
      
      if(rspde.order == 0){
          warning("The order cannot be zero for nonstationary models! The order was changed to 1.")
          rspde.order <- 1
        }

        graph_opt <- get.sparsity.graph.rspde(
        fem_mesh_matrices = fem_mesh, dim = d,
        nu = nu,
        rspde.order = rspde.order,
        force_non_integer = TRUE
      )

    graph_opt <- transpose_cgeneric(graph_opt) 


    model <- do.call(eval(parse(text='INLA::inla.cgeneric.define')),
        list(model="inla_cgeneric_rspde_nonstat_fixed_model",
            shlib=rspde_lib,
            n=as.integer(n_cgeneric)*(rspde.order+1), debug=debug,
            d = as.double(d),
            r_ratapprox = as.vector(r),
            p_ratapprox = as.vector(p),
            k_ratapprox = k,
            nu = nu,
            graph_opt_i = graph_opt@i,
            graph_opt_j = graph_opt@j,
            C = C,
            G = G,
            B_tau = B.tau,
            B_kappa = B.kappa,
            rspde_order = as.integer(rspde.order),
            start.theta = start.theta,
            theta.prior.mean = param$theta.prior.mean,
            theta.prior.prec = param$theta.prior.prec,
            prior.theta.param = prior.theta.param
            ))
    
    model$cgeneric_type <- "frac_alpha"

    } else{

        graph_opt <- get.sparsity.graph.rspde(
        fem_mesh_matrices = fem_mesh, dim = d,
        nu = nu,
        rspde.order = rspde.order
      )
      graph_opt <- transpose_cgeneric(graph_opt) 

    model <- do.call(eval(parse(text='INLA::inla.cgeneric.define')),
        list(model="inla_cgeneric_rspde_nonstat_int_model",
            shlib=rspde_lib,
            n=as.integer(n_cgeneric), debug=debug,
            graph_opt_i = graph_opt@i,
            graph_opt_j = graph_opt@j,
            alpha = as.integer(alpha),
            C = C,
            G = G,
            B_tau = B.tau,
            B_kappa = B.kappa,
            start.theta = start.theta,
            theta.prior.mean = param$theta.prior.mean,
            theta.prior.prec = param$theta.prior.prec,
            prior.theta.param = prior.theta.param
            ))
    
    model$cgeneric_type <- "int_alpha"

    }



### END OF NONSTATIONARY PART
  }


  model$nu <- nu
  model$theta.prior.mean <- theta.prior.mean
  model$prior.nu <- prior.nu
  model$theta.prior.prec <- theta.prior.prec
  model$start.nu <- start.nu
  model$integer.nu <- ifelse(fixed_nu, integer_alpha, FALSE)
  model$start.theta <- start.theta
  model$stationary <- stationary
  if (integer.nu) {
    rspde.order <- 0
  }
  model$rspde.order <- rspde.order
  class(model) <- c("inla_rspde", class(model))
  model$dim <- d
  model$est_nu <- !fixed_nu
  model$n.spde <- mesh$n
  model$nu.upper.bound <- nu.upper.bound
  model$prior.nu.dist <- prior.nu.dist
  model$debug <- debug
  model$type.rational.approx <- type.rational.approx
  model$mesh <- mesh
  model$fem_mesh <- fem_mesh_orig
  model$parameterization <- parameterization
  return(model)
}

#' @noRd 

transpose_cgeneric <- function(Cmatrix){
    Cmatrix <- INLA::inla.as.sparse(Cmatrix)
    # Cmatrix <- as(Cmatrix, "TsparseMatrix")
    ii <- Cmatrix@i
    Cmatrix@i <- Cmatrix@j
    Cmatrix@j <- ii
    idx <- which(Cmatrix@i <= Cmatrix@j)
    Cmatrix@i <- Cmatrix@i[idx]
    Cmatrix@j <- Cmatrix@j[idx]
    Cmatrix@x <- Cmatrix@x[idx]
    return(Cmatrix)
}


#' @noRd 

restructure_matrices_less <- function(matrices_less, m_alpha){
  n_temp <- length(matrices_less)
  temp_vec <- numeric(n_temp)
  N_temp <- n_temp/(m_alpha+1)
  for(i in 1:N_temp){
    for(j in 1:(m_alpha+1)){
      temp_vec[(m_alpha+1)*(i-1) + j] <- matrices_less[(j-1) * N_temp + i]
    }
  }
  return(temp_vec)
}


#' @name spde.make.A
#' @title Observation/prediction matrices for rSPDE models with integer smoothness.
#' @description Constructs observation/prediction weight matrices
#' for rSPDE models with integer smoothness based on `inla.mesh` or
#' `inla.mesh.1d` objects.
#' @param mesh An `inla.mesh`,
#' an `inla.mesh.1d` object or a `metric_graph` object.
#' @param loc Locations, needed if an INLA mesh is provided
#' @param A The A matrix from the standard SPDE approach, such as the matrix
#' returned by `inla.spde.make.A`. Should only be provided if
#' `mesh` is not provided.
#' @param index For each observation/prediction value, an index into loc. Default is seq_len(nrow(A.loc)).
#' @param group For each observation/prediction value, an index into
#' the group model.
#' @param repl For each observation/prediction value, the replicate index.
#' @param n.group The size of the group model.
#' @param n.repl The total number of replicates.
#' @return The \eqn{A} matrix for rSPDE models.
#' @export
#' @examples
#' \donttest{ #devel version
#' if (requireNamespace("fmesher", quietly = TRUE)){
#' library(fmesher)
#' 
#' set.seed(123)
#' loc <- matrix(runif(100 * 2) * 100, 100, 2)
#' mesh <- fm_mesh_2d(
#'   loc = loc,
#'   cutoff = 50,
#'   max.edge = c(50, 500)
#' )
#' A <- spde.make.A(mesh, loc = loc)
#' }
#' #devel.tag
#' }
spde.make.A <- function(mesh = NULL,
                         loc = NULL,
                         A = NULL,
                         index = NULL,
                         group = NULL,
                         repl = 1L,
                         n.group = NULL,
                         n.repl = NULL) {
  if (!is.null(mesh)) {
    cond1 <- inherits(mesh, "inla.mesh.1d")
    cond2 <- inherits(mesh, "inla.mesh")
    cond3 <- inherits(mesh, "metric_graph")
    stopifnot(cond1 || cond2 || cond3)
    if(cond1 || cond2){
      dim <- get_inla_mesh_dimension(mesh)  
    } else if(cond3){
      dim <- 1
    }
  } else if (is.null(dim)) {
    stop("If mesh is not provided, then you should provide the dimension d!")
  }
  if (!is.null(mesh)) {
    if (is.null(loc) && !inherits(mesh, "metric_graph")) {
      stop("If you provided mesh, you should also provide the locations, loc.")
    }
  }

  if (!is.null(mesh)) {
    if(cond1 || cond2){
      # A <- fmesher::fm_basis(
      #   x = mesh, loc = loc, repl = repl)

      A <- fm_basis(
        x = mesh, loc = loc, repl = repl)

        if(!is.null(index)){
          A <- A[index,]
        }

        if(!is.null(n.group)){
          A <- kronecker(Matrix::Diagonal(n.group), A)
        } else{
          if(is.null(group)){
            group <- 1L
          }
          # blk_grp <- fmesher::fm_block(group)
          # A <- fmesher::fm_row_kron(Matrix::t(blk_grp), A)
          blk_grp <- fm_block(group)
          A <- fm_row_kron(Matrix::t(blk_grp), A)
        }

        if(!is.null(n.repl)){
          A <- kronecker(Matrix::Diagonal(n.repl), A)
        } else{
          # blk_rep <- fmesher::fm_block(repl)
          # A <- fmesher::fm_row_kron(Matrix::t(blk_rep), A)
          blk_rep <- fm_block(repl)
          A <- fm_row_kron(Matrix::t(blk_rep), A)
        }

    } else if(cond3){
      if(is.null(mesh$mesh)){
        stop("The graph object should contain a mesh!")
      }
      if(!is.null(group) || !is.null(n.group)){
        stop("Groups are still not implemented for metric graphs.")
      }
      if(!is.null(n.repl)){
         if(is.null(loc)){
          A <- kronecker(Matrix::Diagonal(n.repl), mesh$fem_basis(mesh$get_PtE()))
        } else{
          A <- kronecker(Matrix::Diagonal(n.repl), mesh$fem_basis(loc))
        }
      } else if(!is.null(index)){

          if(min(repl)!= 1){
            stop("The indexes of the replicates should begin at 1!")
          }
          if(any(!is.integer(repl))){
            stop("The indexes of the replicates should be integers!")
          }

        if(is.null(loc)){
          loc_PtE <- mesh$get_PtE()
        } else{
          loc_PtE <- loc
        }


        if(max(repl) == 1){
          A <- mesh$fem_basis(loc_PtE[index,])
        } else{
          stopifnot(length(index) == length(repl))

          if(max(abs(diff(repl)))>1){
            stop("The indexes of the replicates should increase by steps of size 1!")
          }

          total_repl <- max(repl)
          index_tmp <- index[repl==1]
          A <- mesh$fem_basis(loc_PtE[index_tmp,])
          for(i in 2:total_repl){
            index_tmp <- index[repl==i]
            A_tmp <- mesh$fem_basis(loc_PtE[index_tmp,])
            A <- bdiag(A, A_tmp)
          }
          if(any(diff(repl)) < 0){
            col_indexes <- 1:ncol(A)
            new_col_indexes <- col_indexes[repl==1]
            for(i in 2:total_repl){
              new_col_indexes <- c(new_col_indexes, col_indexes[repl==i])
            }
            A <- A[,new_col_indexes]
          }
        }
      } else if(length(repl)>1){
        stop("When using replicates, you should provide index!")
      } else{
        if(is.null(loc)){
          A <- mesh$fem_basis(mesh$get_PtE())
        } else{
          A <- mesh$fem_basis(loc)
        }
      }
    }
  } else if (is.null(A)) {
    stop("If mesh is not provided, then you should provide the A matrix from
         the standard SPDE approach!")
  }



  attr(A, "inla_rspde_Amatrix") <- TRUE
  attr(A, "integer_nu") <- TRUE
  return(A)
}


#' @name rspde.make.A
#' @title Observation/prediction matrices for rSPDE models.
#' @description Constructs observation/prediction weight matrices
#' for rSPDE models based on `inla.mesh` or
#' `inla.mesh.1d` objects.
#' @param mesh An `inla.mesh`,
#' an `inla.mesh.1d` object or a `metric_graph` object.
#' @param loc Locations, needed if an INLA mesh is provided
#' @param A The A matrix from the standard SPDE approach, such as the matrix
#' returned by `inla.spde.make.A`. Should only be provided if
#' `mesh` is not provided.
#' @param dim the dimension. Should only be provided if an
#' `mesh` is not provided.
#' @param rspde.order The order of the covariance-based rational SPDE approach.
#' @param nu If `NULL`, then the model will assume that nu will
#' be estimated. If nu is fixed, you should provide the value of nu.
#' @param index For each observation/prediction value, an index into loc. Default is seq_len(nrow(A.loc)).
#' @param group For each observation/prediction value, an index into
#' the group model.
#' @param repl For each observation/prediction value, the replicate index.
#' @param n.group The size of the group model.
#' @param n.repl The total number of replicates.
#' @return The \eqn{A} matrix for rSPDE models.
#' @export
#' @examples
#' \donttest{ #devel version
#' if (requireNamespace("INLA", quietly = TRUE)){
#' library(INLA)
#' 
#' set.seed(123)
#' loc <- matrix(runif(100 * 2) * 100, 100, 2)
#' mesh <- inla.mesh.2d(
#'   loc = loc,
#'   cutoff = 50,
#'   max.edge = c(50, 500)
#' )
#' A <- rspde.make.A(mesh, loc = loc, rspde.order = 3)
#' }
#' #devel.tag
#' }
rspde.make.A <- function(mesh = NULL,
                         loc = NULL,
                         A = NULL,
                         dim = NULL,
                         rspde.order = 2, nu = NULL,
                         index = NULL,
                         group = NULL,
                         repl = 1L,
                         n.group = NULL,
                         n.repl = NULL) {
  if (!is.null(mesh)) {
    cond1 <- inherits(mesh, "inla.mesh.1d")
    cond2 <- inherits(mesh, "inla.mesh")
    cond3 <- inherits(mesh, "metric_graph")
    stopifnot(cond1 || cond2 || cond3)
    if(cond1 || cond2){
      dim <- get_inla_mesh_dimension(mesh)  
    } else if(cond3){
      dim <- 1
    }
  } else if (is.null(dim)) {
    stop("If mesh is not provided, then you should provide the dimension d!")
  }
  if (!is.null(mesh)) {
    if (is.null(loc) && !inherits(mesh, "metric_graph")) {
      stop("If you provided mesh, you should also provide the locations, loc.")
    }
  }

  if (!is.null(mesh)) {
    if(cond1 || cond2){
      # A <- fmesher::fm_basis(
      #   x = mesh, loc = loc, repl = repl)
      A <- fm_basis(
        x = mesh, loc = loc)

        if(!is.null(index)){
          A <- A[index,]
        }

        if(!is.null(n.group)){
          A <- kronecker(Matrix::Diagonal(n.group), A)
        } else{
          if(is.null(group)){
            group <- 1L
          }
          # blk_grp <- fmesher::fm_block(group)
          # A <- fmesher::fm_row_kron(Matrix::t(blk_grp), A)
          blk_grp <- fm_block(group)
          A <- fm_row_kron(Matrix::t(blk_grp), A)
        }

        if(!is.null(n.repl)){
          A <- kronecker(Matrix::Diagonal(n.repl), A)
        } else{
          # blk_rep <- fmesher::fm_block(repl)
          # A <- fmesher::fm_row_kron(Matrix::t(blk_rep), A)
          blk_rep <- fm_block(repl)
          A <- fm_row_kron(Matrix::t(blk_rep), A)
        }

    } else if(cond3){
      if(is.null(mesh$mesh)){
        stop("The graph object should contain a mesh!")
      }

        if(is.null(loc)){
          A <- mesh$fem_basis(mesh$get_PtE())
        } else{
          A <- mesh$fem_basis(loc)
        }

        if(!is.null(index)){
          A <- A[index,]
        }

        if(!is.null(n.group)){
          A <- kronecker(Matrix::Diagonal(n.group), A)
        } else{
          if(is.null(group)){
            group <- 1L
          }
          blk_grp <- fm_block(group)
          A <- fm_row_kron(Matrix::t(blk_grp), A)
        }

        if(!is.null(n.repl)){
          A <- kronecker(Matrix::Diagonal(n.repl), A)
        } else{
          blk_rep <- fm_block(repl)
          A <- fm_row_kron(Matrix::t(blk_rep), A)
        }        
    }
  } else if (is.null(A)) {
    stop("If mesh is not provided, then you should provide the A matrix from
         the standard SPDE approach!")
  }


  if (!is.null(nu)) {
    if (!is.numeric(nu)) {
      stop("nu must be numeric!")
    }
  }

  fixed_nu <- !is.null(nu)
  if (fixed_nu) {
    alpha <- nu + dim / 2
    integer_alpha <- (alpha %% 1 == 0)
    if (integer_alpha) {
      Abar <- A
      integer_nu <- TRUE
    } else {
    if(rspde.order > 0){
      Abar <- kronecker(matrix(1, 1, rspde.order + 1), A)
    } else{
      Abar <- A
    }
      integer_nu <- FALSE
    }
  } else {
    if(rspde.order > 0){
      Abar <- kronecker(matrix(1, 1, rspde.order + 1), A)
    } else{
      Abar <- A
    }
    integer_nu <- FALSE
  }

  if (integer_nu) {
    rspde.order <- 0
  }


  attr(Abar, "inla_rspde_Amatrix") <- TRUE
  attr(Abar, "rspde.order") <- rspde.order
  attr(Abar, "integer_nu") <- integer_nu
  return(Abar)
}


#' @name rspde.make.index
#' @title rSPDE model index vector generation
#' @description Generates a list of named index vectors for an rSPDE model.
#' @param name A character string with the base name of the effect.
#' @param mesh An `inla.mesh`,
#' an `inla.mesh.1d` object or a `metric_graph` object.
#' @param rspde.order The order of the rational approximation
#' @param nu If `NULL`, then the model will assume that nu will
#' be estimated. If nu is fixed, you should provide the value of nu.
#' @param n.spde The number of basis functions in the mesh model.
#' @param n.group The size of the group model.
#' @param n.repl The total number of replicates.
#' @param dim the dimension of the domain. Should only be provided if
#' `mesh` is not provided.
#' @return A list of named index vectors.
#' \item{name}{Indices into the vector of latent variables}
#' \item{name.group}{'group' indices}
#' \item{name.repl}{Indices for replicates}
#' @export
#' @examples
#' \donttest{ #devel version
#' if (requireNamespace("INLA", quietly = TRUE)){
#' library(INLA)
#' 
#' set.seed(123)
#'
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
#' Abar <- rspde.make.A(mesh = mesh_2d, loc = loc_2d_mesh)
#' mesh.index <- rspde.make.index(name = "field", mesh = mesh_2d)
#' st.dat <- inla.stack(
#'   data = list(y = as.vector(y)),
#'   A = Abar,
#'   effects = mesh.index
#' )
#' rspde_model <- rspde.matern(
#'   mesh = mesh_2d,
#'   nu.upper.bound = 2
#' )
#' f <- y ~ -1 + f(field, model = rspde_model)
#' rspde_fit <- inla(f,
#'   data = inla.stack.data(st.dat),
#'   family = "gaussian",
#'   control.predictor =
#'     list(A = inla.stack.A(st.dat))
#' )
#' result <- rspde.result(rspde_fit, "field", rspde_model)
#' summary(result)
#' }
#' #devel.tag
#' }
rspde.make.index <- function(name, n.spde = NULL, n.group = 1,
                             n.repl = 1, mesh = NULL,
                             rspde.order = 2, nu = NULL, dim = NULL) {
  if (is.null(n.spde) && is.null(mesh)) {
    stop("You should provide either n.spde or mesh!")
  }

  if (!is.null(mesh)) {
    cond1 <- inherits(mesh, "inla.mesh.1d")
    cond2 <- inherits(mesh, "inla.mesh")
    cond3 <- inherits(mesh, "metric_graph")
    stopifnot(cond1 || cond2 || cond3)
    if(cond1 || cond2){
      n_mesh <- mesh$n
      dim <- get_inla_mesh_dimension(mesh)
    } else if(cond3){
      dim <- 1
      # n_mesh <- nrow(mesh$mesh$VtE)
      n_mesh <- nrow(mesh$mesh$VtE)
    }
  } else {
    n_mesh <- n.spde
    if (is.null(dim)) {
      stop("You should provide the dimension d!")
    }
  }

  name.group <- paste(name, ".group", sep = "")
  name.repl <- paste(name, ".repl", sep = "")

  if (!is.null(nu)) {
    if (!is.numeric(nu)) {
      stop("nu must be numeric!")
    }
  }

  fixed_nu <- !is.null(nu)

  if (fixed_nu) {
    alpha <- nu + dim / 2
    integer_alpha <- (alpha %% 1 == 0)

    if (integer_alpha) {
      factor_rspde <- 1
      integer_nu <- TRUE
    } else {
      if(rspde.order > 0){
        factor_rspde <- rspde.order + 1
      } else{
        factor_rspde <- 1
      }
      integer_nu <- FALSE
    }
  } else {
      if(rspde.order > 0){
        factor_rspde <- rspde.order + 1
      } else{
        factor_rspde <- 1
      }
    integer_nu <- FALSE
  }

  out <- list()
  out[[name]] <- as.vector(sapply(1:factor_rspde, function(i) {
    rep(rep(((i - 1) * n_mesh + 1):(i * n_mesh), times = n.group),
    times = n.repl)
  }))
  out[[name.group]] <- rep(rep(rep(1:n.group, each = n_mesh),
  times = n.repl), times = factor_rspde)
  out[[name.repl]] <- rep(rep(1:n.repl, each = n_mesh * n.group),
  times = factor_rspde)
  class(out) <- c("inla_rspde_index", class(out))
  if (integer_nu) {
    rspde.order <- 0
  }
  attr(out, "rspde.order") <- rspde.order
  attr(out, "integer_nu") <- integer_nu
  attr(out, "n.mesh") <- n_mesh
  attr(out, "name") <- name
  attr(out, "n.group") <- n.group
  attr(out, "n.repl") <- n.repl
  return(out)
}



#' Data extraction from metric graphs for 'rSPDE' models
#'
#' Extracts data from metric graphs to be used by 'INLA' and 'inlabru'.
#'
#' @param graph_rspde An `inla_metric_graph_spde` object built with the
#' `rspde.metric_graph()` function.
#' @param name A character string with the base name of the effect.
#' @param repl Which replicates? If there is no replicates, one
#' can set `repl` to `NULL`. If one wants all replicates,
#' then one sets to `repl` to `.all`.
#' @param group Which groups? If there is no groups, one
#' can set `group` to `NULL`. If one wants all groups,
#' then one sets to `group` to `.all`.
#' @param group_col Which "column" of the data contains the group variable?
#' @param only_pred Should only return the `data.frame` to the prediction data?
#' @param loc  Locations. If not given, they will be chosen as the available locations on the metric graph internal dataset.
#' @param loc_name Character with the name of the location variable to be used in
#' 'inlabru' prediction.
#' @param tibble Should the data be returned as a `tidyr::tibble`?
#' @param drop_na Should the rows with at least one NA for one of the columns be removed? DEFAULT is `FALSE`. This option is turned to `FALSE` if `only_pred` is `TRUE`.
#' @param drop_all_na Should the rows with all variables being NA be removed? DEFAULT is `TRUE`. This option is turned to `FALSE` if `only_pred` is `TRUE`.
#' @return An 'INLA' and 'inlabru' friendly list with the data.
#' @export

graph_data_rspde <- function (graph_rspde, name = "field", repl = NULL, group = NULL, 
                                group_col = NULL,
                                only_pred = FALSE,
                                loc = NULL,
                                loc_name = NULL,
                                tibble = FALSE,
                                drop_na = FALSE, drop_all_na = TRUE){

  ret <- list()

  rspde.order <- graph_rspde$rspde.order

  nu <- graph_rspde$nu

  graph_tmp <- graph_rspde$mesh$clone()
  if(is.null((graph_tmp$.__enclos_env__$private$data))){
    stop("The graph has no data!")
  }
  if(only_pred){
    idx_anyNA <- !idx_not_any_NA(graph_tmp$.__enclos_env__$private$data)
    graph_tmp$.__enclos_env__$private$data <- lapply(graph_tmp$.__enclos_env__$private$data, function(dat){return(dat[idx_anyNA])})
    drop_na <- FALSE
    drop_all_na <- FALSE
  }

  if(is.null(repl)){
    groups <- graph_tmp$.__enclos_env__$private$data[[".group"]]
    repl <- groups[1]
  } else if(repl[1] == ".all") {
    groups <- graph_tmp$.__enclos_env__$private$data[[".group"]]
    repl <- unique(groups)
  } 

   ret[["data"]] <- select_repl_group(graph_tmp$.__enclos_env__$private$data, repl = repl, group = group, group_col = group_col)   

   repl_vec <- ret[["data"]][[".group"]]
   if(!is.null(group_col)){
    group_vec <- ret[["data"]][[group_col]]
    group <- unique(group_vec)
   } else{
    group_vec <- rep(1, length(ret[["data"]][[".group"]]))
    group <- 1
   }


  n.repl <- length(unique(repl))

  if(is.null(group)){
    n.group <- 1
  } else if (group[1] == ".all"){
    n.group <- length(unique(graph_tmp$.__enclos_env__$private$data[[group_col]]))
  } else{
    n.group <- length(unique(group))
  }

  if(tibble){
    ret[["data"]] <-tidyr::as_tibble(ret[["data"]])
  }

  if(drop_all_na){
    is_tbl <- inherits(ret, "tbl_df")
      idx_temp <- idx_not_all_NA(ret[["data"]])
      ret[["data"]] <- lapply(ret[["data"]], function(dat){dat[idx_temp]}) 
      if(is_tbl){
        ret[["data"]] <- tidyr::as_tibble(ret[["data"]])
      }
  }    
  if(drop_na){
    if(!inherits(ret[["data"]], "tbl_df")){
      idx_temp <- idx_not_any_NA(ret[["data"]])
      ret[["data"]] <- lapply(ret[["data"]], function(dat){dat[idx_temp]})
    } else{
      ret[["data"]] <- tidyr::drop_na(ret[["data"]])
    }
  }
  
  if(!is.null(loc_name)){
      ret[["data"]][[loc_name]] <- cbind(ret[["data"]][[".edge_number"]],
                          ret[["data"]][[".distance_on_edge"]])
  }


  if(!inherits(ret[["data"]], "metric_graph_data")){
    class(ret[["data"]]) <- c("metric_graph_data", class(ret))
  }



  ret[["index"]] <- rspde.make.index(mesh = graph_tmp, n.group = n.group, n.repl = n.repl, nu = nu, dim = 1, rspde.order = rspde.order, name = name)

  ret[["repl"]] <- ret[["data"]][[".group"]]

   if(!is.null(group_col)){
    group_vec <- ret[["data"]][[group_col]]
    group <- unique(group_vec)
   } else{
    group_vec <- rep(1, length(ret[["data"]][[".group"]]))
    group <- 1
   }

  repl_vec <- ret[["repl"]]   

  n_obs <- sum(ret[["data"]][[".group"]] == ret[["data"]][[".group"]][1])

  # index_basis <- rep(rep(1:n_obs, times = n.group), 
  #           times = n.repl)

  ret[["basis"]] <- Matrix::Matrix(nrow=0,ncol=0)       

  loc_basis <- cbind(ret[["data"]][[".edge_number"]], ret[["data"]][[".distance_on_edge"]])     
  
  # We assume the data is ordered by group, then repl. This will be handled by the advanced grouping we are implementing

  for(repl_ in repl){
    idx_rep <- (repl_vec == repl_)
    for(group_ in group){
      idx_grp <- (group_vec == group_)
      idx_grp_rep <- as.logical(idx_grp * idx_rep)
      ret[["basis"]] <- Matrix::bdiag(ret[["basis"]], graph_tmp$fem_basis(loc_basis[idx_grp_rep,]))
    }
  }

  if(!graph_rspde$integer.nu){
    ret[["basis"]] <- kronecker(matrix(1, 1, rspde.order + 1), 
                ret[["basis"]])
  }
  
  return(ret)
}



#' Extraction of vector of replicates for 'INLA'
#'
#' Extracts the vector of replicates from an 'rSPDE'
#' model object for 'INLA'
#'
#' @param graph_spde An `rspde_metric_graph` object built with the
#' `rspde.metric_graph()` function from the 'rSPDE' package.
#' @param repl Which replicates? If there is no replicates, one
#' can set `repl` to `NULL`. If one wants all replicates,
#' then one sets to `repl` to `.all`.
#' @param group Which groups? If there is no groups, one
#' can set `group` to `NULL`. If one wants all groups,
#' then one sets to `group` to `.all`.
#' @param group_col Which "column" of the data contains the group variable?
#' @return The vector of replicates paired with groups.
#' @noRd

graph_repl_rspde <- function (graph_spde, repl = NULL, group = NULL, group_col = NULL){
  # graph_tmp <- graph_spde$graph_spde$clone()
  if(is.null(repl)){
    # groups <- graph_tmp$data[[".group"]]
    groups <- graph_spde$graph_spde$.__enclos_env__$private$data[[".group"]]
    repl <- groups[1] 
  } else if(repl[1] == ".all") {
    # ret <- graph_tmp$data
    groups <- graph_spde$graph_spde$.__enclos_env__$private$data[[".group"]]
    repl <- unique(groups)
  } 

  ret <- select_repl_group(graph_spde$graph_spde$.__enclos_env__$private$data, repl = repl, group = group, group_col = group_col)    
  return(ret[[".group"]])
}

#' Select replicate and group
#' @noRd
#'
select_repl_group <- function(data_list, repl, group, group_col){
    if(!is.null(group) && is.null(group_col)){
      stop("If you specify group, you need to specify group_col!")
    }
    if(!is.null(group)){
      grp <- data_list[[group_col]]
      grp <- which(grp %in% group)
      data_result <- lapply(data_list, function(dat){dat[grp]})
      replicates <- data_result[[".group"]]
      replicates <- which(replicates %in% repl)
      data_result <- lapply(data_result, function(dat){dat[replicates]})
      return(data_result)
    } else{
      replicates <- data_list[[".group"]]
      replicates <- which(replicates %in% repl)
      data_result <- lapply(data_list, function(dat){dat[replicates]})
      return(data_result)
    }
}

#' Creation of index vector for 'INLA'
#'
#' Creates the vector of indexes from an 'rSPDE'
#' model object for 'INLA'
#'
#' @param graph_spde An `rspde_metric_graph` object built with the
#' `rspde.metric_graph()` function from the 'rSPDE' package.
#' @param n.repl The total number of replicates.
#' @param n.group The size of the group model.
#' @return The vector of indexes.
#' @noRd

graph_index_rspde <- function (graph_spde, n.repl = 1, n.group = 1){
        n_obs <- nrow(graph_spde$mesh$get_PtE())
        rep(rep(1:n_obs, times = n.group), 
            times = n.repl)
}


#' @name rspde.result
#' @title rSPDE result extraction from INLA estimation results
#' @description Extract field and parameter values and distributions
#' for an rspde effect from an inla result object.
#' @param inla An `inla` object obtained from a call to
#' `inla()`.
#' @param name A character string with the name of the rSPDE effect
#' in the inla formula.
#' @param rspde The `inla_rspde` object used for the effect in
#' the inla formula.
#' @param compute.summary Should the summary be computed?
#' @param parameterization If 'detect', the parameterization from the model will be used. Otherwise, the options are 'spde', 'matern' and 'matern2'.
#' @param n_samples The number of samples to be used if parameterization is different from the one used to fit the model.
#' @param n_density The number of equally spaced points to estimate the density.
#' @return If the model was fitted with `matern` parameterization (the default), it returns a list containing:
#' \item{marginals.range}{Marginal densities for the range parameter}
#' \item{marginals.log.range}{Marginal densities for log(range)}
#' \item{marginals.std.dev}{Marginal densities for std. deviation}
#' \item{marginals.log.std.dev}{Marginal densities for log(std. deviation)}
#' \item{marginals.values}{Marginal densities for the field values}
#' \item{summary.log.range}{Summary statistics for log(range)}
#' \item{summary.log.std.dev}{Summary statistics for log(std. deviation)}
#' \item{summary.values}{Summary statistics for the field values}
#' If `compute.summary` is `TRUE`, then the list will also contain
#' \item{summary.kappa}{Summary statistics for kappa}
#' \item{summary.tau}{Summary statistics for tau}
#' If the model was fitted with the `spde` parameterization, it returns a list containing:
#' \item{marginals.kappa}{Marginal densities for kappa}
#' \item{marginals.log.kappa}{Marginal densities for log(kappa)}
#' \item{marginals.log.tau}{Marginal densities for log(tau)}
#' \item{marginals.tau}{Marginal densities for tau}
#' \item{marginals.values}{Marginal densities for the field values}
#' \item{summary.log.kappa}{Summary statistics for log(kappa)}
#' \item{summary.log.tau}{Summary statistics for log(tau)}
#' \item{summary.values}{Summary statistics for the field values}
#' If `compute.summary` is `TRUE`, then the list will also contain
#' \item{summary.kappa}{Summary statistics for kappa}
#' \item{summary.tau}{Summary statistics for tau}
#' 
#' For both cases, if nu was estimated, then the list will also contain
#' \item{marginals.nu}{Marginal densities for nu}
#' If nu was estimated and a beta prior was used, then the list will
#' also contain
#' \item{marginals.logit.nu}{Marginal densities for logit(nu)}
#' \item{summary.logit.nu}{Marginal densities for logit(nu)}
#' If nu was estimated and a truncated lognormal prior was used,
#' then the list will also contain
#' \item{marginals.log.nu}{Marginal densities for log(nu)}
#' \item{summary.log.nu}{Marginal densities for log(nu)}
#' If nu was estimated and `compute.summary` is `TRUE`,
#' then the list will also contain
#' \item{summary.nu}{Summary statistics for nu}
#' @export
#' @examples
#' \donttest{ #devel version
#' if (requireNamespace("INLA", quietly = TRUE)){
#' library(INLA)
#' 
#' set.seed(123)
#'
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
#' Abar <- rspde.make.A(mesh = mesh_2d, loc = loc_2d_mesh)
#' mesh.index <- rspde.make.index(name = "field", mesh = mesh_2d)
#' st.dat <- inla.stack(
#'   data = list(y = as.vector(y)),
#'   A = Abar,
#'   effects = mesh.index
#' )
#' rspde_model <- rspde.matern(
#'   mesh = mesh_2d,
#'   nu.upper.bound = 2
#' )
#' f <- y ~ -1 + f(field, model = rspde_model)
#' rspde_fit <- inla(f,
#'   data = inla.stack.data(st.dat),
#'   family = "gaussian",
#'   control.predictor =
#'     list(A = inla.stack.A(st.dat))
#' )
#' result <- rspde.result(rspde_fit, "field", rspde_model)
#' summary(result)
#' }
#' #devel.tag
#' }
rspde.result <- function(inla, name, rspde, compute.summary = TRUE, parameterization = "detect", n_samples = 5000, n_density = 1024) {
  check_class_inla_rspde(rspde)

  stationary <- rspde$stationary

  parameterization <- parameterization[[1]]
  parameterization <- tolower(parameterization)

  if(!(parameterization %in% c("detect", "spde", "matern", "matern2"))){
    stop("The possible options for parameterization are 'detect', 'spde', 'matern' and 'matern2'.")
  }

  nu.upper.bound <- rspde$nu.upper.bound
  result <- list()

  par_model <- rspde$parameterization

  if(parameterization == "detect"){
    parameterization <- rspde$parameterization
  } 

  if(stationary){
          if (!rspde$est_nu) {
              if(parameterization == "spde"){
                row_names <- c("tau", "kappa")
              } else if (parameterization == "matern") {
                row_names <- c("std.dev", "range")
              } else if (parameterization == "matern2") {
                row_names <- c("var", "r")
              }
            } else {
              if(parameterization == "spde"){
                row_names <- c("tau", "kappa", "nu")
              } else if (parameterization == "matern") {
                row_names <- c("std.dev", "range", "nu")
              } else if (parameterization == "matern2") {
                row_names <- c("var", "r", "nu")
              }
            }


            result$summary.values <- inla$summary.random[[name]]

            if (!is.null(inla$marginals.random[[name]])) {
              result$marginals.values <- inla$marginals.random[[name]]
            }

            if(parameterization == "spde"){
              name_theta1 <- "tau"
              name_theta2 <- "kappa"
              } else if (parameterization == "matern") {
              name_theta1 <- "std.dev"
              name_theta2 <- "range"
              } else if (parameterization == "matern2") {
              name_theta1 <- "var"
              name_theta2 <- "r"
            }

            if(par_model == "spde"){
              name_theta1_model <- "tau"
              name_theta2_model <- "kappa"
              } else if (par_model == "matern") {
              name_theta1_model <- "std.dev"
              name_theta2_model <- "range"
              } else if (par_model == "matern2") {
              name_theta1_model <- "var"
              name_theta2_model <- "r"
            }            

            result[[paste0("summary.log.",name_theta1_model)]] <- INLA::inla.extract.el(
              inla$summary.hyperpar,
              paste("Theta1 for ", name, "$", sep = "")
            )
            rownames(  result[[paste0("summary.log.",name_theta1_model)]]) <- paste0("log(",name_theta1_model,")")

            result[[paste0("summary.log.",name_theta2_model)]] <- INLA::inla.extract.el(
              inla$summary.hyperpar,
              paste("Theta2 for ", name, "$", sep = "")
            )
            rownames(result[[paste0("summary.log.",name_theta2_model)]]) <- paste0("log(", name_theta2_model,")")
            if (rspde$est_nu) {
              result$summary.logit.nu <- INLA::inla.extract.el(
                inla$summary.hyperpar,
                paste("Theta3 for ", name, "$", sep = "")
              )
              rownames(result$summary.logit.nu) <- "logit(nu)"
            }

            if (!is.null(inla$marginals.hyperpar[[paste0("Theta1 for ", name)]])) {
              result[[paste0("marginals.log.",name_theta1_model)]] <- INLA::inla.extract.el(
                inla$marginals.hyperpar,
                paste("Theta1 for ", name, "$", sep = "")
              )
              names(result[[paste0("marginals.log.",name_theta1_model)]]) <- name_theta1_model
              result[[paste0("marginals.log.",name_theta2_model)]] <- INLA::inla.extract.el(
                inla$marginals.hyperpar,
                paste("Theta2 for ", name, "$", sep = "")
              )
              names(result[[paste0("marginals.log.",name_theta2_model)]]) <- name_theta2_model

              if (rspde$est_nu) {
                result$marginals.logit.nu <- INLA::inla.extract.el(
                  inla$marginals.hyperpar,
                  paste("Theta3 for ", name, "$", sep = "")
                )
                names(result$marginals.logit.nu) <- "nu"
              }


              if(par_model == parameterization){
                        result[[paste0("marginals.",name_theta1)]] <- lapply(
                          result[[paste0("marginals.log.",name_theta1)]],
                          function(x) {
                            INLA::inla.tmarginal(
                              function(y) exp(y),
                              x
                            )
                          }
                        )
                        result[[paste0("marginals.",name_theta2)]] <- lapply(
                          result[[paste0("marginals.log.",name_theta2)]],
                          function(x) {
                            INLA::inla.tmarginal(
                              function(y) exp(y),
                              x
                            )
                          }
                        )
                        if (rspde$est_nu) {
                          result$marginals.nu <- lapply(
                            result$marginals.logit.nu,
                            function(x) {
                              INLA::inla.tmarginal(
                                function(y) {
                                  nu.upper.bound * exp(y) / (1 + exp(y))
                                },
                                x
                              )
                            }
                          )
                        }
              } else{
                if(par_model == "spde"){
                        dim <- rspde$dim
                        hyperpar_sample <- INLA::inla.hyperpar.sample(n_samples, inla)
                        if (rspde$est_nu) {
                          nu_est <- rspde$nu.upper.bound * exp(hyperpar_sample[, paste0('Theta3 for ',name)])/(1+exp(hyperpar_sample[, paste0('Theta3 for ',name)]))
                        } else{
                          nu_est <- rspde[["nu"]]
                        }
                        tau_est <- exp(hyperpar_sample[, paste0('Theta1 for ',name)])
                        kappa_est <- exp(hyperpar_sample[, paste0('Theta2 for ',name)])

                        sigma_est <- sqrt(gamma(0.5) / (tau_est^2 * kappa_est^(2 * nu_est) *
                                (4 * pi)^(dim / 2) * gamma(nu_est + dim / 2)))

                        if(parameterization == "matern"){
                              range_est <- sqrt(8 * nu_est)/kappa_est
                              density_theta1 <- stats::density(sigma_est, n = n_density)
                              density_theta2 <- stats::density(range_est, n = n_density)            
                        } else if (parameterization == "matern2"){
                              var_est <- sigma_est^2
                              r_est <- 1/kappa_est
                              density_theta1 <- stats::density(var_est, n = n_density)
                              density_theta2 <- stats::density(r_est, n = n_density)
                        }

                              result[[paste0("marginals.",name_theta1)]] <- list()
                              result[[paste0("marginals.",name_theta1)]][[name_theta1]] <- cbind(density_theta1$x, density_theta1$y)
                              colnames(result[[paste0("marginals.",name_theta1)]][[name_theta1]]) <- c("x","y")

                              result[[paste0("marginals.",name_theta2)]] <- list()
                              result[[paste0("marginals.",name_theta2)]][[name_theta2]] <- cbind(density_theta2$x, density_theta2$y)
                              colnames(result[[paste0("marginals.",name_theta2)]][[name_theta2]]) <- c("x","y")    
                } else if(par_model == "matern"){
                        hyperpar_sample <- INLA::inla.hyperpar.sample(n_samples, inla)
                       dim <- rspde$dim
                        if (rspde$est_nu) {
                          nu_est <- rspde$nu.upper.bound * exp(hyperpar_sample[, paste0('Theta3 for ',name)])/(1+exp(hyperpar_sample[, paste0('Theta3 for ',name)]))
                        } else{
                          nu_est <- rspde[["nu"]]
                        }
                        sigma_est <- exp(hyperpar_sample[, paste0('Theta1 for ',name)])
                        range_est <- exp(hyperpar_sample[, paste0('Theta2 for ',name)])

                        kappa_est <- sqrt(8 * nu_est)/range_est

                        if(parameterization == "spde"){
                              tau_est <- sqrt(gamma(0.5) / (sigma_est^2 * kappa_est^(2 * nu_est) *
                                (4 * pi)^(dim / 2) * gamma(nu_est + dim / 2)))
                              density_theta1 <- stats::density(tau_est, n = n_density)
                              density_theta2 <- stats::density(kappa_est, n = n_density)              
                        } else if (parameterization == "matern2"){
                              var_est <- sigma_est^2
                              r_est <- 1/kappa_est
                              density_theta1 <- stats::density(var_est, n = n_density)
                              density_theta2 <- stats::density(r_est, n = n_density)
                        }

                              result[[paste0("marginals.",name_theta1)]] <- list()
                              result[[paste0("marginals.",name_theta1)]][[name_theta1]] <- cbind(density_theta1$x, density_theta1$y)
                              colnames(result[[paste0("marginals.",name_theta1)]][[name_theta1]]) <- c("x","y")

                              result[[paste0("marginals.",name_theta2)]] <- list()
                              result[[paste0("marginals.",name_theta2)]][[name_theta2]] <- cbind(density_theta2$x, density_theta2$y)
                              colnames(result[[paste0("marginals.",name_theta2)]][[name_theta2]]) <- c("x","y")    
                } else if(par_model == "matern2"){
                        hyperpar_sample <- INLA::inla.hyperpar.sample(n_samples, inla)
                       dim <- rspde$dim
                        if (rspde$est_nu) {
                          nu_est <- rspde$nu.upper.bound * exp(hyperpar_sample[, paste0('Theta3 for ',name)])/(1+exp(hyperpar_sample[, paste0('Theta3 for ',name)]))
                        } else{
                          nu_est <- rspde[["nu"]]
                        }
                        var_est <- exp(hyperpar_sample[, paste0('Theta1 for ',name)])
                        r_est <- exp(hyperpar_sample[, paste0('Theta2 for ',name)])

                        kappa_est <- 1/r_est
                        sigma_est <- sqrt(var_est)

                        if(parameterization == "spde"){
                              tau_est <- sqrt(gamma(0.5) / (sigma_est^2 * kappa_est^(2 * nu_est) *
                                (4 * pi)^(dim / 2) * gamma(nu_est + dim / 2)))
                              density_theta1 <- stats::density(tau_est, n = n_density)
                              density_theta2 <- stats::density(kappa_est, n = n_density)              
                        } else if (parameterization == "matern"){
                              range_est <- sqrt(8 * nu_est)/kappa_est
                              density_theta1 <- stats::density(sigma_est, n = n_density)
                              density_theta2 <- stats::density(range_est, n = n_density)
                        }

                              result[[paste0("marginals.",name_theta1)]] <- list()
                              result[[paste0("marginals.",name_theta1)]][[name_theta1]] <- cbind(density_theta1$x, density_theta1$y)
                              colnames(result[[paste0("marginals.",name_theta1)]][[name_theta1]]) <- c("x","y")

                              result[[paste0("marginals.",name_theta2)]] <- list()
                              result[[paste0("marginals.",name_theta2)]][[name_theta2]] <- cbind(density_theta2$x, density_theta2$y)
                              colnames(result[[paste0("marginals.",name_theta2)]][[name_theta2]]) <- c("x","y")   
                }

                        if (rspde$est_nu) {
                          result$marginals.nu <- lapply(
                            result$marginals.logit.nu,
                            function(x) {
                              INLA::inla.tmarginal(
                                function(y) {
                                  nu.upper.bound * exp(y) / (1 + exp(y))
                                },
                                x
                              )
                            }
                          )
                        }

              }

            }





            if (compute.summary) {
              norm_const <- function(density_df) {
                min_x <- min(density_df[, "x"])
                max_x <- max(density_df[, "x"])
                denstemp <- function(x) {
                  dens <- sapply(x, function(z) {
                    if (z < min_x) {
                      return(0)
                    } else if (z > max_x) {
                      return(0)
                    } else {
                      return(approx(x = density_df[, "x"],
                      y = density_df[, "y"], xout = z)$y)
                    }
                  })
                  return(dens)
                }
                norm_const <- stats::integrate(
                  f = function(z) {
                    denstemp(z)
                  }, lower = min_x, upper = max_x,
                  subdivisions = nrow(density_df),
                  stop.on.error = FALSE
                )$value
                return(norm_const)
              }

              norm_const_theta1 <- norm_const(result[[paste0("marginals.",name_theta1)]][[name_theta1]])
              result[[paste0("marginals.",name_theta1)]][[name_theta1]][, "y"] <-
              result[[paste0("marginals.",name_theta1)]][[name_theta1]][, "y"] / norm_const_theta1

              norm_const_theta2 <- norm_const(result[[paste0("marginals.",name_theta2)]][[name_theta2]])
              result[[paste0("marginals.",name_theta2)]][[name_theta2]][, "y"] <-
              result[[paste0("marginals.",name_theta2)]][[name_theta2]][, "y"] / norm_const_theta2

              result[[paste0("summary.",name_theta1)]] <- create_summary_from_density(result[[paste0("marginals.",name_theta1)]][[name_theta1]],
              name = name_theta1)

              result[[paste0("summary.",name_theta2)]] <-
              create_summary_from_density(result[[paste0("marginals.",name_theta2)]][[name_theta2]], name = name_theta2)

              if (rspde$est_nu) {
                norm_const_nu <- norm_const(result$marginals.nu$nu)
                result$marginals.nu$nu[, "y"] <-
                result$marginals.nu$nu[, "y"] / norm_const_nu

                result$summary.nu <- create_summary_from_density(result$marginals.nu$nu,
                name = "nu")
              }
            }
    } else{
      n_par <- length(rspde$start.theta)

      if (!rspde$est_nu) {
        if(parameterization == "spde"){
                row_names <- sapply(1:n_par, function(i){paste0("Theta",i,".spde")})
         } else if(parameterization == "matern"){
                row_names <- sapply(1:n_par, function(i){paste0("Theta",i,".matern")})
         } else if(parameterization == "matern2"){
                row_names <- sapply(1:n_par, function(i){paste0("Theta",i,".matern2")})
         }
      } else {
       if(parameterization == "spde"){
              row_names <- sapply(1:n_par, function(i){paste0("Theta",i,".spde")})
              row_names <- c(row_names, "nu")
         } else if(parameterization == "matern"){
               row_names <- sapply(1:n_par, function(i){paste0("Theta",i,".matern")})
               row_names <- c(row_names, "nu")
         } else if(parameterization == "matern2"){
               row_names <- sapply(1:n_par, function(i){paste0("Theta",i,".matern2")})
               row_names <- c(row_names, "nu")
         }
       }

       for(i in 1:n_par){
          result[[paste0("summary.",row_names[i])]] <- INLA::inla.extract.el(
              inla$summary.hyperpar,
              paste0("Theta",i," for ", name, "$", sep = "")
            )
            rownames(  result[[paste0("summary.",row_names[i])]]) <- row_names[i]
       }

            if (rspde$est_nu) {
              result$summary.logit.nu <- INLA::inla.extract.el(
                inla$summary.hyperpar,
                paste0("Theta",n_par+1," for ", name, "$", sep = "")
              )
              rownames(result$summary.logit.nu) <- "logit(nu)"
            }


       for(i in 1:n_par){
            if (!is.null(inla$marginals.hyperpar[[paste0("Theta",i," for ", name)]])) {
              result[[paste0("marginals.",row_names[i])]] <- INLA::inla.extract.el(
                inla$marginals.hyperpar,
                paste0("Theta",i," for ", name, "$", sep = "")
              )
              names(result[[paste0("marginals.",row_names[i])]]) <- row_names[i]
            }
       }


              if (rspde$est_nu) {
                result$marginals.logit.nu <- INLA::inla.extract.el(
                  inla$marginals.hyperpar,
                  paste0("Theta",n_par+1," for ", name, "$", sep = "")
                )
                names(result$marginals.logit.nu) <- "nu"

                result$marginals.nu <- lapply(
                  result$marginals.logit.nu,
                  function(x) {
                    INLA::inla.tmarginal(
                      function(y) {
                        nu.upper.bound * exp(y) / (1 + exp(y))
                      },
                      x
                    )
                  }
                )
              }

            if (compute.summary) {
              norm_const <- function(density_df) {
                min_x <- min(density_df[, "x"])
                max_x <- max(density_df[, "x"])
                denstemp <- function(x) {
                  dens <- sapply(x, function(z) {
                    if (z < min_x) {
                      return(0)
                    } else if (z > max_x) {
                      return(0)
                    } else {
                      return(approx(x = density_df[, "x"],
                      y = density_df[, "y"], xout = z)$y)
                    }
                  })
                  return(dens)
                }
                norm_const <- stats::integrate(
                  f = function(z) {
                    denstemp(z)
                  }, lower = min_x, upper = max_x,
                  stop.on.error = FALSE
                )$value
                return(norm_const)
              }

              if (rspde$est_nu) {
                norm_const_nu <- norm_const(result$marginals.nu$nu)
                result$marginals.nu$nu[, "y"] <-
                result$marginals.nu$nu[, "y"] / norm_const_nu

                result$summary.nu <- create_summary_from_density(result$marginals.nu$nu,
                name = "nu")
              }
            }


    }

  result$n_par <- length(rspde$start.theta)

  class(result) <- "rspde_result"
  result$stationary <- stationary
  result$parameterization <- parameterization
  if(stationary){
    result$params <- c(name_theta1,name_theta2)
   if(rspde$est_nu){
    result$params <- c(result$params, "nu")
  }
  } else {
    result$params <- row_names
  }

  return(result)
}


#' @name gg_df
#' @title Data frame for result objects from R-INLA fitted models to be used in ggplot2
#' @param result a result object for which the data frame is desired
#' @param ... further arguments passed to or from other methods.
#' @return A data frame containing the posterior densities.
#'
#' @rdname gg_df
#' @export
gg_df <- function(result, ...){
UseMethod("gg_df", result)
}




#' Data frame for rspde_result objects to be used in ggplot2
#'
#' Returns a ggplot-friendly data-frame with the marginal posterior densities.
#' 
#' @name gg_df.rspde_result
#' @param result An rspde_result object.
#' @param parameter Vector. Which parameters to get the posterior density in the data.frame? The options are `std.dev`, `range`, `tau`, `kappa` and `nu`.
#' @param transform Should the posterior density be given in the original scale?
#' @param restrict_x_axis Variables to restrict the range of x axis based on quantiles.
#' @param restrict_quantiles Named list of quantiles to restrict x axis. It should contain the name of the parameter
#' along with a vector with two elements specifying the lower and upper quantiles. The names should be
#' match the ones in result$params. For example, if we want to restrict nu to the 0.05 and 0.95 quantiles
#' we do `restrict_quantiles = c(0.05, 0.95)`.
#' @param ... currently not used.
#'
#' @return A data frame containing the posterior densities.
#' @export
gg_df.rspde_result <- function(result, 
                          parameter = result$params,
                          transform = TRUE,
                          restrict_x_axis = NULL,
                          restrict_quantiles = NULL,
                          ...) {
      rspde_result <- result
      parameter <- intersect(parameter, result$params)
      if(length(parameter) == 0){
        stop("You should choose at least one of the parameters. The available parameters are given in result$params!")
      }
    if ("nu" %in% parameter) {
    if (is.null(rspde_result$marginals.nu)) {
      parameter <- parameter[parameter != "nu"]
    }
  }

  stationary <- result$stationary

  if(stationary){
        param <- parameter[[1]]
        if(transform){
          param <- paste0("marginals.", param)
        } else{
          if(param != "nu"){
            param <- paste0("marginals.log.", param)
          } else{
            param <- paste0("marginals.logit.", param)
          }
        }
        ret_df <- data.frame(x = rspde_result[[param]][[parameter[1]]][,1], 
        y = rspde_result[[param]][[parameter[1]]][,2], 
        parameter = parameter[[1]])

        if(parameter[[1]] %in% restrict_x_axis){
          if(is.null( restrict_quantiles[[parameter[[1]]]])){
             restrict_quantiles[[parameter[[1]]]] <- c(0,1)
          }
          d_t <- c(0,diff(ret_df$x))
          emp_cdf <- cumsum(d_t*ret_df$y)
          lower_quant <- restrict_quantiles[[parameter[[1]]]][1]
          upper_quant <- restrict_quantiles[[parameter[[1]]]][2]
          filter_coord <- (emp_cdf >= lower_quant) * (emp_cdf <= upper_quant)
          filter_coord <- as.logical(filter_coord)
          ret_df <- ret_df[filter_coord, ]
        }

        if(length(parameter) > 1){
        for(i in 2:length(parameter)){
        param <- parameter[[i]]
        if(transform){
          param <- paste0("marginals.", param)
        } else{
          if(param != "nu"){
            param <- paste0("marginals.log.", param)
          } else{
            param <- paste0("marginals.logit.", param)
          }
        }
          tmp <- data.frame(x = rspde_result[[param]][[parameter[i]]][,1], 
            y = rspde_result[[param]][[parameter[i]]][,2], 
            parameter = parameter[[i]])

          if(parameter[[i]] %in% restrict_x_axis){
          if(is.null( restrict_quantiles[[parameter[[i]]]])){
             restrict_quantiles[[parameter[[i]]]] <- c(0,1)
          }
          d_t <- c(0,diff(tmp$x))
          emp_cdf <- cumsum(d_t*tmp$y)
          lower_quant <- restrict_quantiles[[parameter[[i]]]][1]
          upper_quant <- restrict_quantiles[[parameter[[i]]]][2]
          filter_coord <- (emp_cdf >= lower_quant) * (emp_cdf <= upper_quant)
          filter_coord <- as.logical(filter_coord)
          tmp <- tmp[filter_coord, ]
        }

          ret_df <- rbind(ret_df, tmp)
        }
        }
  }  else{

        param <- parameter[[1]]
        if(transform && param == "nu"){
          param <- paste0("marginals.", param)
        } else if (param == "nu"){
            param <- paste0("marginals.logit.", param)
        } else{
            param <- paste0("marginals.", param)
        }
        ret_df <- data.frame(x = rspde_result[[param]][[parameter[1]]][,1], 
        y = rspde_result[[param]][[parameter[1]]][,2], 
        parameter = parameter[[1]])

        if(parameter[[1]] %in% restrict_x_axis){
          if(is.null( restrict_quantiles[[parameter[[1]]]])){
             restrict_quantiles[[parameter[[1]]]] <- c(0,1)
          }
          d_t <- c(0,diff(ret_df$x))
          emp_cdf <- cumsum(d_t*ret_df$y)
          lower_quant <- restrict_quantiles[[parameter[[1]]]][1]
          upper_quant <- restrict_quantiles[[parameter[[1]]]][2]
          filter_coord <- (emp_cdf >= lower_quant) * (emp_cdf <= upper_quant)
          filter_coord <- as.logical(filter_coord)
          ret_df <- ret_df[filter_coord, ]
        }

        if(length(parameter) > 1){
        for(i in 2:length(parameter)){
        param <- parameter[[i]]
        if(transform && param == "nu"){
          param <- paste0("marginals.", param)
        } else if (param == "nu"){
            param <- paste0("marginals.logit.", param)
        } else{
            param <- paste0("marginals.", param)
        }
          tmp <- data.frame(x = rspde_result[[param]][[parameter[i]]][,1], 
            y = rspde_result[[param]][[parameter[i]]][,2], 
            parameter = parameter[[i]])

          if(parameter[[i]] %in% restrict_x_axis){
          if(is.null( restrict_quantiles[[parameter[[i]]]])){
             restrict_quantiles[[parameter[[i]]]] <- c(0,1)
          }
          d_t <- c(0,diff(tmp$x))
          emp_cdf <- cumsum(d_t*tmp$y)
          lower_quant <- restrict_quantiles[[parameter[[i]]]][1]
          upper_quant <- restrict_quantiles[[parameter[[i]]]][2]
          filter_coord <- (emp_cdf >= lower_quant) * (emp_cdf <= upper_quant)
          filter_coord <- as.logical(filter_coord)
          tmp <- tmp[filter_coord, ]
        }

          ret_df <- rbind(ret_df, tmp)
        }
        }

  }

  return(ret_df)
}


#' @name summary.rspde_result
#' @title Summary for posteriors of field parameters for an `inla_rspde`
#' model from a `rspde_result` object
#' @description Summary for posteriors of rSPDE field parameters in
#' their original scales.
#' @param object A `rspde_result` object.
#' @param digits integer, used for number formatting with signif()
#' @param ... Currently not used.
#' @return Returns a `data.frame`
#' containing the summary.
#' @export
#' @method summary rspde_result
#' @examples
#' \donttest{ #devel version
#' if (requireNamespace("INLA", quietly = TRUE)){
#' library(INLA)
#' 
#' set.seed(123)
#'
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
#' Abar <- rspde.make.A(mesh = mesh_2d, loc = loc_2d_mesh)
#' mesh.index <- rspde.make.index(name = "field", mesh = mesh_2d)
#' st.dat <- inla.stack(
#'   data = list(y = as.vector(y)),
#'   A = Abar,
#'   effects = mesh.index
#' )
#' rspde_model <- rspde.matern(
#'   mesh = mesh_2d,
#'   nu.upper.bound = 2
#' )
#' f <- y ~ -1 + f(field, model = rspde_model)
#' rspde_fit <- inla(f,
#'   data = inla.stack.data(st.dat),
#'   family = "gaussian",
#'   control.predictor =
#'     list(A = inla.stack.A(st.dat))
#' )
#' result <- rspde.result(rspde_fit, "field", rspde_model)
#' summary(result)
#' }
#' #devel.tag
#' }
#'
summary.rspde_result <- function(object,
                                 digits = 6,
                                 ...) {

  if (is.null(object[[paste0("summary.",object$params[1])]])) {
    warning("The summary was not computed, rerun rspde_result with
    compute.summary set to TRUE.")
  } else {
    n_par <- object$n_par
    out <- object[[paste0("summary.",object$params[1])]]
    if(n_par > 1){
      for(i in 2:n_par){
        out <- rbind(out, object[[paste0("summary.",object$params[i])]])
      }
    }

    if (!is.null(object$summary.nu)) {
      out <- rbind(out, object$summary.nu)
    }
    return(signif(out, digits = digits))
  }
}




#' @name rspde.mesh.project
#' @title Calculate a lattice projection to/from an `inla.mesh` for
#' rSPDE objects
#' @aliases rspde.mesh.project rspde.mesh.projector rspde.mesh.project.inla.mesh
#' rspde.mesh.project.rspde.mesh.projector rspde.mesh.project.inla.mesh.1d
#' @description Calculate a lattice projection to/from an `inla.mesh` for
#' rSPDE objects
#' @param mesh An `inla.mesh` or `inla.mesh.1d` object.
#' @param nu The smoothness parameter. If `NULL`, it will be assumed that
#' nu was estimated.
#' @param rspde.order The order of the rational approximation.
#' @param loc	Projection locations. Can be a matrix or a SpatialPoints or a
#' SpatialPointsDataFrame object.
#' @param field Basis function weights, one per mesh basis function, describing
#' the function to be evaluated at the projection locations.
#' @param projector A `rspde.mesh.projector` object.
#' @param lattice An `inla.mesh.lattice` object.
#' @param xlim X-axis limits for a lattice. For R2 meshes, defaults to covering
#' the domain.
#' @param ylim Y-axis limits for a lattice. For R2 meshes, defaults to covering
#' the domain.
#' @param dims Lattice dimensions.
#' @param projection One of c("default", "longlat", "longsinlat", "mollweide").
#' @param ... Additional parameters.
#' @return A list with projection information for rspde.mesh.project. For
#' rspde.mesh.projector(mesh, ...),
#' a rspde.mesh.projector object. For rspde.mesh.project(projector, field, ...),
#' a field projected from the mesh onto the locations
#' given by the projector object.
#' @details This function is built upon the inla.mesh.project and
#' inla.mesh.projector functions from INLA.
#' @rdname rspde.mesh.project
#' @export
#'
rspde.mesh.project <- function(...) {
  UseMethod("rspde.mesh.project")
}

#' @rdname rspde.mesh.project
#' @export

rspde.mesh.projector <- function(mesh,
                                 nu = NULL,
                                 rspde.order = 2,
                                 loc = NULL,
                                 lattice = NULL,
                                 xlim = NULL,
                                 ylim = NULL,
                                 dims = c(100, 100),
                                 projection = NULL,
                                 ...) {
  args_list <- list()
  args_list[["mesh"]] <- mesh
  if (!is.null(loc)) {
    args_list[["loc"]] <- loc
  }
  if (!is.null(lattice)) {
    args_list[["lattice"]] <- lattice
  }
  if (!is.null(xlim)) {
    args_list[["xlim"]] <- xlim
  }
  if (!is.null(ylim)) {
    args_list[["ylim"]] <- ylim
  }
  if (!is.null(projection)) {
    args_list[["projection"]] <- projection
  }
  args_list[["dims"]] <- dims
  # out <- do.call(INLA::inla.mesh.projector, args_list)
  out <- do.call(fmesher::fm_evaluator, args_list)
  dim <- get_inla_mesh_dimension(mesh)

  out$proj$A <- rspde.make.A(
    A = out$proj$A, rspde.order = rspde.order, dim = dim,
    nu = nu
  )

  class(out) <- c("rspde.mesh.projector",class(out))
  return(out)
}


#' @rdname rspde.mesh.project
#' @method rspde.mesh.project inla.mesh
#' @export

rspde.mesh.project.inla.mesh <- function(mesh, loc = NULL,
                                         field = NULL, rspde.order = 2,
                                         nu = NULL, ...) {
  cond1 <- inherits(mesh, "inla.mesh.1d")
  cond2 <- inherits(mesh, "inla.mesh")
  stopifnot(cond1 || cond2)

  if (!missing(field) && !is.null(field)) {
    proj <- rspde.mesh.projector(mesh,
      loc = loc, rspde.order = rspde.order, nu = nu,
      ...
    )
    # return(INLA::inla.mesh.project(proj, field = field))
    return(fmesher::fm_evaluate(proj, field = field))
  }
  jj <- which(rowSums(matrix(is.na(as.vector(loc)),
    nrow = nrow(loc),
    ncol = ncol(loc)
  )) == 0)
  # smorg <- (INLA::inla.fmesher.smorg(mesh$loc,
  # mesh$graph$tv, points2mesh = loc[jj, ,
  #   drop = FALSE
  # ]))
  smorg <- fmesher::fm_bary(mesh,loc=mesh$loc)
  ti <- matrix(0L, nrow(loc), 1)
  b <- matrix(0, nrow(loc), 3)
  ti[jj, 1L] <- smorg$p2m.t
  b[jj, ] <- smorg$p2m.b
  ok <- (ti[, 1L] > 0L)
  ti[ti[, 1L] == 0L, 1L] <- NA
  ii <- which(ok)
  A <- (sparseMatrix(dims = c(nrow(loc), mesh$n), i = rep(
    ii,
    3
  ), j = as.vector(mesh$graph$tv[ti[ii, 1L], ]), x = as.vector(b[ii, ])))

  if (!is.null(nu)) {
    if (!is.numeric(nu)) {
      stop("nu must be numeric!")
    }
  }

  fixed_nu <- !is.null(nu)
  if (fixed_nu) {
    alpha <- nu + 1
    integer_alpha <- (alpha %% 1 == 0)
    if (integer_alpha) {
      Abar <- A
    } else {
      Abar <- kronecker(matrix(1, 1, rspde.order + 1), A)
    }
  } else {
    Abar <- kronecker(matrix(1, 1, rspde.order + 1), A)
  }

  list(t = ti, bary = b, A = Abar, ok = ok)
}


#' @rdname rspde.mesh.project
#' @method rspde.mesh.project rspde.mesh.projector
#' @export
#'

rspde.mesh.project.rspde.mesh.projector <- function(projector, field, ...) {
  # return(INLA::inla.mesh.project(projector = projector, field = field, ...))
  return(fmesher::fm_evaluate(projector = projector, field = field, ...))
}



#' @rdname rspde.mesh.project
#' @method rspde.mesh.project inla.mesh.1d
#' @export
#'

rspde.mesh.project.inla.mesh.1d <- function(mesh, loc, field = NULL,
                                            rspde.order = 2, nu = NULL, ...) {
  stopifnot(inherits(mesh, "inla.mesh.1d"))
  if (!missing(field) && !is.null(field)) {
    proj <- rspde.mesh.projector(mesh, loc,
    rspde.order = rspde.order, nu = nu, ...)
    # return(INLA::inla.mesh.project(proj, field))
    return(fmesher::fm_evaluate(proj, field))
  }
  # A <- INLA::inla.mesh.1d.A(mesh, loc = loc)
  A <- fmesher::fm_basis(mesh, loc = loc)
  if (!is.null(nu)) {
    if (!is.numeric(nu)) {
      stop("nu must be numeric!")
    }
  }

  fixed_nu <- !is.null(nu)
  if (fixed_nu) {
    alpha <- nu + 1 / 2
    integer_alpha <- (alpha %% 1 == 0)
    if (integer_alpha) {
      Abar <- A
    } else {
      Abar <- kronecker(matrix(1, 1, rspde.order + 1), A)
    }
  } else {
    Abar <- kronecker(matrix(1, 1, rspde.order + 1), A)
  }
  return(list(A = Abar, ok = (loc >= mesh$interval[1]) & (loc <=
    mesh$interval[2])))
}



#' @name rspde.matern.precision.opt
#' @title Optimized precision matrix of the covariance-based rational
#' approximation
#' @description `rspde.matern.precision` is used for computing the
#' optimized version of the precision matrix of the
#' covariance-based rational SPDE approximation of a stationary Gaussian random
#' fields on \eqn{R^d} with a Matern covariance function
#' \deqn{C(h) = \frac{\sigma^2}{2^{\nu-1}\Gamma(\nu)}(\kappa h)^\nu
#' K_\nu(\kappa h).}
#' @param kappa Range parameter of the covariance function.
#' @param tau Scale parameter of the covariance function.
#' @param nu Shape parameter of the covariance function.
#' @param rspde.order The order of the rational approximation
#' @param dim The dimension of the domain
#' @param fem_matrices A list containing the FEM-related matrices.
#' The list should contain elements C, G, G_2, G_3, etc.
#' @param graph The sparsity graph of the matrices. If NULL, only a vector
#' of the elements will be returned, if non-NULL, a sparse matrix will
#' be returned.
#' @param sharp The sparsity graph should have the correct sparsity (costs
#' more to perform a sparsity analysis) or an upper bound for the sparsity?
#' @param type_rational_approx Which type of rational approximation
#' should be used? The current types are "chebfun", "brasil" or "chebfunLB".
#' @return The precision matrix
#' @export

rspde.matern.precision.opt <- function(kappa, nu, tau, rspde.order,
dim, fem_matrices, graph = NULL, sharp, type_rational_approx) {
  n_m <- rspde.order

  mt <- get_rational_coefficients(n_m, type_rational_approx)

  beta <- nu / 2 + dim / 4

  m_alpha <- floor(2 * beta)

  # r <- sapply(1:(n_m), function(i) {
  #   approx(mt$alpha, mt[[paste0("r", i)]], cut_decimals(2 * beta))$y
  # })

  # p <- sapply(1:(n_m), function(i) {
  #   approx(mt$alpha, mt[[paste0("p", i)]], cut_decimals(2 * beta))$y
  # })

  # k <- approx(mt$alpha, mt$k, cut_decimals(2 * beta))$y

  row_nu <- round(1000*cut_decimals(2*beta))
  r <- unlist(mt[row_nu, 2:(1+rspde.order)])
  p <- unlist(mt[row_nu, (2+rspde.order):(1+2*rspde.order)])
  k <- unlist(mt[row_nu, 2+2*rspde.order])


  if (m_alpha == 0) {
    L <- (fem_matrices[["C"]] + fem_matrices[["G"]] / (kappa^2))
    Q <- (L - p[1] * fem_matrices[["C"]]) / r[1]
    if (length(r) > 1) {
      for (i in 2:length(r)) {
        Q <- c(Q, (L - p[i] * fem_matrices[["C"]]) / r[i])
      }
    }
  } else {
    if (m_alpha == 1) {
      Malpha <- (fem_matrices[["C"]] + fem_matrices[["G"]] / (kappa^2))
    } else if (m_alpha > 1) {
      Malpha <- fem_matrices[["C"]] + m_alpha * fem_matrices[["G"]] / (kappa^2)
      for (j in 2:m_alpha) {
        Malpha <- Malpha + choose(m_alpha, j) *
        fem_matrices[[paste0("G_", j)]] / (kappa^(2 * j))
      }
    } else {
      stop("Something is wrong with the value of nu!")
    }


    if (m_alpha == 1) {
      Malpha2 <- (fem_matrices[["G"]] + fem_matrices[["G_2"]] / (kappa^2))
    } else if (m_alpha > 1) {
      Malpha2 <- fem_matrices[["G"]] + m_alpha *
      fem_matrices[["G_2"]] / (kappa^2)
      for (j in 2:m_alpha) {
        Malpha2 <- Malpha2 + choose(m_alpha, j) *
        fem_matrices[[paste0("G_", j + 1)]] / (kappa^(2 * j))
      }
    } else {
      stop("Something is wrong with the value of nu!")
    }

    Q <- 1 / r[1] * (Malpha + Malpha2 / kappa^2 - p[1] * Malpha)

    if (length(r) > 1) {
      for (i in 2:length(r)) {
        Q <- c(Q, 1 / r[i] * (Malpha + Malpha2 / kappa^2 - p[i] * Malpha))
      }
    }
  }

  # add k_part into Q

  if (sharp) {
    if (m_alpha == 0) {
      Kpart <- fem_matrices[["C_less"]]
      idx_nonzero <- (Kpart != 0)
      Kpart[idx_nonzero] <- 1/Kpart[idx_nonzero]
      Kpart <- Kpart / k
    } else {
      if (m_alpha == 1) {
        Kpart <- 1 / k * (fem_matrices[["C_less"]] +
        fem_matrices[["G_less"]] / (kappa^2))
      } else if (m_alpha > 1) {
        Kpart <- 1 / k * fem_matrices[["C_less"]] +
        1 / k * m_alpha * fem_matrices[["G_less"]] / (kappa^2)
        for (j in 2:m_alpha) {
          Kpart <- Kpart + 1 / k * choose(m_alpha, j) *
          fem_matrices[[paste0("G_", j, "_less")]] / (kappa^(2 * j))
        }
      } else {
        stop("Something is wrong with the value of nu!")
      }
    }
  } else {
    if (m_alpha == 0) {
      Kpart <- fem_matrices[["C"]]
      idx_nonzero <- (Kpart != 0)
      Kpart[idx_nonzero] <- 1/Kpart[idx_nonzero]      
      Kpart <- Kpart / k
    } else {
      if (m_alpha == 1) {
        Kpart <- 1 / k * (fem_matrices[["C"]] + fem_matrices[["G"]] / (kappa^2))
      } else if (m_alpha > 1) {
        Kpart <- 1 / k * fem_matrices[["C"]] +
        1 / k * m_alpha * fem_matrices[["G"]] / (kappa^2)
        for (j in 2:m_alpha) {
          Kpart <- Kpart + 1 / k * choose(m_alpha, j) *
          fem_matrices[[paste0("G_", j)]] / (kappa^(2 * j))
        }
      } else {
        stop("Something is wrong with the value of nu!")
      }
    }
  }




  Q <- c(Q, Kpart)

  Q <- Q * kappa^(4 * beta)

  Q <- tau^2 * Q

  if (!is.null(graph)) {
    # graph <- as(graph, "dgTMatrix")
    graph <- as(graph,"TsparseMatrix")
    idx <- which(graph@i <= graph@j)
    Q <- Matrix::sparseMatrix(
      i = graph@i[idx], j = graph@j[idx], x = Q,
      symmetric = TRUE, index1 = FALSE
    )
  }

  return(Q)
}

#' @name rspde.matern.precision
#' @title Precision matrix of the covariance-based rational approximation of
#' stationary Gaussian Matern random fields
#' @description `rspde.matern.precision` is used for computing the
#' precision matrix of the
#' covariance-based rational SPDE approximation of a stationary Gaussian random
#' fields on \eqn{R^d} with a Matern covariance function
#' \deqn{C(h) = \frac{\sigma^2}{2^(\nu-1)\Gamma(\nu)}(\kappa h)^\nu
#' K_\nu(\kappa h)}{C(h) = (\sigma^2/(2^(\nu-1)\Gamma(\nu))(\kappa h)^\nu
#' K_\nu(\kappa h)}
#' @param kappa Range parameter of the covariance function.
#' @param tau Scale parameter of the covariance function. If sigma is
#' not provided, tau should be provided.
#' @param sigma Standard deviation of the covariance function. If tau is
#' not provided, sigma should be provided.
#' @param nu Shape parameter of the covariance function.
#' @param rspde.order The order of the rational approximation
#' @param dim The dimension of the domain
#' @param fem_mesh_matrices A list containing the FEM-related matrices. The
#' list should contain elements c0, g1, g2, g3, etc.
#' @param only_fractional Logical. Should only the fractional-order part of
#' the precision matrix be returned?
#' @param return_block_list Logical. For `type = "covariance"`, should the
#' block parts of the precision matrix be returned separately as a list?
#' @param type_rational_approx Which type of rational approximation should be
#' used? The current types are "chebfun", "brasil" or "chebfunLB".
#'
#' @return The precision matrix
#' @export
#' @examples
#' set.seed(123)
#' nobs <- 101
#' x <- seq(from = 0, to = 1, length.out = nobs)
#' fem <- rSPDE.fem1d(x)
#' kappa <- 40
#' sigma <- 1
#' d <- 1
#' nu <- 2.6
#' tau <- sqrt(gamma(nu) / (kappa^(2 * nu) * (4 * pi)^(d / 2) *
#' gamma(nu + d / 2)))
#' range <- sqrt(8*nu)/kappa
#' op_cov <- matern.operators(
#'   loc_mesh = x, nu = nu, range = range, sigma = sigma,
#'   d = 1, m = 2, compute_higher_order = TRUE,
#'   parameterization = "matern"
#' )
#' v <- t(rSPDE.A1d(x, 0.5))
#' c.true <- matern.covariance(abs(x - 0.5), kappa, nu, sigma)
#' Q <- rspde.matern.precision(
#'   kappa = kappa, nu = nu, tau = tau, rspde.order = 2, d = 1,
#'   fem_mesh_matrices = op_cov$fem_mesh_matrices
#' )
#' A <- Diagonal(nobs)
#' Abar <- cbind(A, A, A)
#' w <- rbind(v, v, v)
#' c.approx_cov <- (Abar) %*% solve(Q, w)
#'
#' # plot the result and compare with the true Matern covariance
#' plot(x, matern.covariance(abs(x - 0.5), kappa, nu, sigma),
#'   type = "l", ylab = "C(h)",
#'   xlab = "h", main = "Matern covariance and rational approximations"
#' )
#' lines(x, c.approx_cov, col = 2)
rspde.matern.precision <- function(kappa, nu, tau = NULL, sigma = NULL,
rspde.order, dim, fem_mesh_matrices,
only_fractional = FALSE, return_block_list = FALSE,
type_rational_approx = "chebfun") {
  if (is.null(tau) && is.null(sigma)) {
    stop("You should provide either tau or sigma!")
  }

  if (is.null(tau)) {
    tau <- sqrt(gamma(nu) / (sigma^2 * kappa^(2 * nu) *
    (4 * pi)^(dim / 2) * gamma(nu + dim / 2)))
  }

  n_m <- rspde.order

  mt <- get_rational_coefficients(n_m, type_rational_approx)

  beta <- nu / 2 + dim / 4

  m_alpha <- floor(2 * beta)

  row_nu <- round(1000*cut_decimals(2*beta))
  r <- unlist(mt[row_nu, 2:(1+rspde.order)])
  p <- unlist(mt[row_nu, (2+rspde.order):(1+2*rspde.order)])
  k <- unlist(mt[row_nu, 2+2*rspde.order])

  if (!only_fractional) {
    if (m_alpha == 0) {
      L <- ((kappa^2) * fem_mesh_matrices[["c0"]] +
      fem_mesh_matrices[["g1"]]) / kappa^2
      Q <- (L - p[1] * fem_mesh_matrices[["c0"]]) / r[1]
      if (length(r) > 1) {
        for (i in 2:length(r)) {
          Q <- bdiag(Q, (L - p[i] * fem_mesh_matrices[["c0"]]) / r[i])
        }
      }
    } else {
      if (m_alpha == 1) {
        Malpha <- (fem_mesh_matrices[["c0"]] +
        fem_mesh_matrices[["g1"]] / (kappa^2))
      } else if (m_alpha > 1) {
        Malpha <- fem_mesh_matrices[["c0"]] + m_alpha *
        fem_mesh_matrices[["g1"]] / (kappa^2)
        for (j in 2:m_alpha) {
          Malpha <- Malpha + choose(m_alpha, j) *
          fem_mesh_matrices[[paste0("g", j)]] / (kappa^(2 * j))
        }
      } else {
        stop("Something is wrong with the value of nu!")
      }


      if (m_alpha == 1) {
        Malpha2 <- (fem_mesh_matrices[["g1"]] +
        fem_mesh_matrices[["g2"]] / (kappa^2))
      } else if (m_alpha > 1) {
        Malpha2 <- fem_mesh_matrices[["g1"]] + m_alpha *
        fem_mesh_matrices[["g2"]] / (kappa^2)
        for (j in 2:m_alpha) {
          Malpha2 <- Malpha2 + choose(m_alpha, j) *
          fem_mesh_matrices[[paste0("g", j + 1)]] / (kappa^(2 * j))
        }
      } else {
        stop("Something is wrong with the value of nu!")
      }

      Q <- 1 / r[1] * (Malpha + Malpha2 / kappa^2 - p[1] * Malpha)

      if (length(r) > 1) {
        for (i in 2:length(r)) {
          Q <- bdiag(Q, 1 / r[i] * (Malpha + Malpha2 / kappa^2 - p[i] * Malpha))
        }
      }
    }


    # add k_part into Q

    if (m_alpha == 0) {
      C <- fem_mesh_matrices[["c0"]]
      Kpart <- Matrix::Diagonal(dim(C)[1], 1 / rowSums(C))
      
      Kpart <- Kpart / k
    } else {
      if (m_alpha == 1) {
        Kpart <- 1 / k * (fem_mesh_matrices[["c0"]] +
        fem_mesh_matrices[["g1"]] / (kappa^2))
      } else if (m_alpha > 1) {
        Kpart <- 1 / k * fem_mesh_matrices[["c0"]] +
        1 / k * m_alpha * fem_mesh_matrices[["g1"]] / (kappa^2)
        for (j in 2:m_alpha) {
          Kpart <- Kpart + 1 / k * choose(m_alpha, j) *
          fem_mesh_matrices[[paste0("g", j)]] / (kappa^(2 * j))
        }
      } else {
        stop("Something is wrong with the value of nu!")
      }
    }

    Q <- bdiag(Q, Kpart)

    Q <- Q * kappa^(4 * beta)

    Q <- tau^2 * Q



    return(Q)
  } else {
    L <- ((kappa^2) * fem_mesh_matrices[["c0"]] +
    fem_mesh_matrices[["g1"]]) / kappa^2

    if (return_block_list) {
      Q <- list()

      Q[[length(Q) + 1]] <- kappa^(4 * beta) * tau^2 *
      (L - p[1] * fem_mesh_matrices[["c0"]]) / r[1]


      if (n_m > 1) {
        for (i in 2:(n_m)) {
          Q[[length(Q) + 1]] <- kappa^(4 * beta) * tau^2 *
          (L - p[i] * fem_mesh_matrices[["c0"]]) / r[i]
        }
      }
      if(m_alpha==0) {
        C <- fem_mesh_matrices[["c0"]]
        Kpart <- Matrix::Diagonal(dim(C)[1], 1 / rowSums(C))
      } else {
        Kpart <- fem_mesh_matrices[["c0"]]
      }
      
      Q[[length(Q) + 1]] <- kappa^(4 * beta) * tau^2 *
      Kpart / k

      return(Q)
    } else {
      Q <- (L - p[1] * fem_mesh_matrices[["c0"]]) / r[1]

      if (n_m > 1) {
        for (i in 2:(n_m)) {
          temp <- (L - p[i] * fem_mesh_matrices[["c0"]]) / r[i]
          Q <- bdiag(Q, temp)
        }
      }
      
      if(m_alpha==0) {
        C <- fem_mesh_matrices[["c0"]]
        Kpart <- Matrix::Diagonal(dim(C)[1], 1 / rowSums(C))
      } else {
        Kpart <- fem_mesh_matrices[["c0"]]
      }
      Q <- bdiag(Q, Kpart / k)


      Q <- Q * kappa^(4 * beta)

      Q <- tau^2 * Q
      return(Q)
    }
  }
}


#' @name rspde.matern.precision.integer.opt
#' @title Optimized precision matrix of stationary Gaussian Matern
#' random fields with integer covariance exponent
#' @description `rspde.matern.precision.integer.opt` is used
#' for computing the optimized version of the precision matrix of
#' stationary Gaussian random fields on \eqn{R^d} with a Matern
#' covariance function
#' \deqn{C(h) = \frac{\sigma^2}{2^{\nu-1}\Gamma(\nu)}(\kappa h)^\nu
#' K_\nu(\kappa h),}
#' where \eqn{\alpha = \nu + d/2} is a natural number.
#' @param kappa Range parameter of the covariance function.
#' @param tau Scale parameter of the covariance function.
#' @param nu Shape parameter of the covariance function.
#' @param d The dimension of the domain
#' @param fem_matrices A list containing the FEM-related matrices.
#' The list should contain elements C, G, G_2, G_3, etc.
#' @param graph The sparsity graph of the matrices. If NULL, only a vector
#' of the elements will be returned, if non-NULL, a sparse matrix will
#' be returned.
#' @return The precision matrix
#' @export

rspde.matern.precision.integer.opt <- function(kappa,
                                               nu,
                                               tau,
                                               d,
                                               fem_matrices,
                                               graph = NULL) {
  beta <- nu / 2 + d / 4

  n_beta <- floor(2 * beta)

  if (n_beta == 1) {
    Q <- (kappa^2 * fem_matrices[["C_less"]] + fem_matrices[["G_less"]])
  } else if (n_beta > 1) {
    Q <- kappa^(2 * n_beta) * fem_matrices[["C_less"]] + n_beta *
    kappa^(2 * (n_beta - 1)) * fem_matrices[["G_less"]]
    for (j in 2:n_beta) {
      Q <- Q + kappa^(2 * (n_beta - j)) * choose(n_beta, j) *
      fem_matrices[[paste0("G_", j, "_less")]]
    }
  } else {
    stop("Something is wrong with the value of nu!")
  }

  Q <- tau^2 * Q

  if (!is.null(graph)) {
    # graph <- as(graph, "dgTMatrix")
    graph <- as(graph,"TsparseMatrix")
    idx <- which(graph@i <= graph@j)
    Q <- Matrix::sparseMatrix(
      i = graph@i[idx], j = graph@j[idx], x = Q,
      symmetric = TRUE, index1 = FALSE
    )
  }

  return(Q)
}

#' @name rspde.matern.precision.integer
#' @title Precision matrix of stationary Gaussian Matern
#' random fields with integer covariance exponent
#' @description `rspde.matern.precision.integer.opt` is
#' used for computing the precision matrix of stationary
#' Gaussian random fields on \eqn{R^d} with a Matern
#' covariance function
#' \deqn{C(h) = \frac{\sigma^2}{2^(\nu-1)\Gamma(\nu)}
#' (\kappa h)^\nu K_\nu(\kappa h)}{C(h) =
#' (\sigma^2/(2^(\nu-1)\Gamma(\nu))(\kappa h)^\nu K_\nu(\kappa h)},
#' where \eqn{\alpha = \nu + d/2} is a natural number.
#' @param kappa Range parameter of the covariance function.
#' @param tau Scale parameter of the covariance function.
#' @param sigma Standard deviation of the covariance function.
#' If tau is not provided, sigma should be provided.
#' @param nu Shape parameter of the covariance function.
#' @param dim The dimension of the domain
#' @param fem_mesh_matrices A list containing the FEM-related
#' matrices. The list should contain elements c0, g1, g2, g3, etc.
#' @return The precision matrix
#' @export
#' @examples
#' set.seed(123)
#' nobs <- 101
#' x <- seq(from = 0, to = 1, length.out = nobs)
#' fem <- rSPDE.fem1d(x)
#' kappa <- 40
#' sigma <- 1
#' d <- 1
#' nu <- 0.5
#' tau <- sqrt(gamma(nu) / (kappa^(2 * nu) *
#' (4 * pi)^(d / 2) * gamma(nu + d / 2)))
#' range <- sqrt(8*nu)/kappa
#' op_cov <- matern.operators(
#'   loc_mesh = x, nu = nu, range = range, sigma = sigma,
#'   d = 1, m = 2, parameterization = "matern"
#' )
#' v <- t(rSPDE.A1d(x, 0.5))
#' c.true <- matern.covariance(abs(x - 0.5), kappa, nu, sigma)
#' Q <- rspde.matern.precision.integer(
#'   kappa = kappa, nu = nu, tau = tau, d = 1,
#'   fem_mesh_matrices = op_cov$fem_mesh_matrices
#' )
#' A <- Diagonal(nobs)
#' c.approx_cov <- A %*% solve(Q, v)
#'
#' # plot the result and compare with the true Matern covariance
#' plot(x, matern.covariance(abs(x - 0.5), kappa, nu, sigma),
#'   type = "l", ylab = "C(h)",
#'   xlab = "h", main = "Matern covariance and rational approximations"
#' )
#' lines(x, c.approx_cov, col = 2)
rspde.matern.precision.integer <- function(kappa, nu, tau = NULL,
sigma = NULL, dim, fem_mesh_matrices) {
  if (is.null(tau) && is.null(sigma)) {
    stop("You should provide either tau or sigma!")
  }

  if (is.null(tau)) {
    tau <- sqrt(gamma(nu) / (sigma^2 * kappa^(2 * nu) *
    (4 * pi)^(dim / 2) * gamma(nu + dim / 2)))
  }

  beta <- nu / 2 + dim / 4

  n_beta <- floor(2 * beta)

  if (n_beta == 1) {
    Q <- (kappa^2 * fem_mesh_matrices[["c0"]] + fem_mesh_matrices[["g1"]])
  } else if (n_beta > 1) {
    Q <- kappa^(2 * n_beta) * fem_mesh_matrices[["c0"]] + n_beta *
    kappa^(2 * (n_beta - 1)) * fem_mesh_matrices[["g1"]]
    for (j in 2:n_beta) {
      Q <- Q + kappa^(2 * (n_beta - j)) * choose(n_beta, j) *
      fem_mesh_matrices[[paste0("g", j)]]
    }
  } else {
    stop("Something is wrong with the value of nu!")
  }

  Q <- tau^2 * Q

  return(Q)
}


#' @name rspde.metric_graph
#' @title Matern rSPDE model object for metric graphs in INLA
#' @description Creates an INLA object for a stationary Matern model on a metric graph with
#' general smoothness parameter.
#' @param graph_obj The graph object to build the model. Needs to be of class `metric_graph`. It should have a built mesh.
#' If the mesh is not built, one will be built using h=0.01 as default.
#' @param h The width of the mesh in case the mesh was not built.
#' @param nu.upper.bound Upper bound for the smoothness parameter.
#' @param rspde.order The order of the covariance-based rational SPDE approach.
#' @param nu If nu is set to a parameter, nu will be kept fixed and will not
#' be estimated. If nu is `NULL`, it will be estimated.
#' @param B.sigma Matrix with specification of log-linear model for \eqn{\sigma}. Will be used if `parameterization = 'matern'`.
#' @param B.range Matrix with specification of log-linear model for \eqn{\rho}, which is a range-like parameter (it is exactly the range parameter in the stationary case). Will be used if `parameterization = 'matern'`.
#' @param parameterization Which parameterization to use? `matern` uses range, std. deviation and nu (smoothness). `spde` uses kappa, tau and nu (smoothness). The default is `matern`.
#' @param B.tau Matrix with specification of log-linear model for \eqn{\tau}. Will be used if `parameterization = 'spde'`.
#' @param B.kappa Matrix with specification of log-linear model for \eqn{\kappa}. Will be used if `parameterization = 'spde'`.
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
#' distribution. Check details below.
#' @param type.rational.approx Which type of rational approximation
#' should be used? The current types are "chebfun", "brasil" or "chebfunLB".
#' @param debug INLA debug argument
#' @param shared_lib Which shared lib to use for the cgeneric implementation? 
#' If "INLA", it will use the shared lib from INLA's installation. If 'rSPDE', then
#' it will use the local installation (does not work if your installation is from CRAN).
#' Otherwise, you can directly supply the path of the .so (or .dll) file.
#'
#' @return An INLA model.
#' @export

rspde.metric_graph <- function(graph_obj,
                         h = NULL,
                         nu.upper.bound = 2, rspde.order = 2,
                         nu = NULL, 
                         debug = FALSE,
                         B.sigma = matrix(c(0, 1, 0), 1, 3), 
                         B.range = matrix(c(0, 0, 1), 1, 3), 
                         parameterization = c("matern", "spde"),
                         B.tau = matrix(c(0, 1, 0), 1, 3), 
                         B.kappa = matrix(c(0, 0, 1), 1, 3), 
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
                         prior.nu.dist = c("lognormal", "beta"),
                         nu.prec.inc = 1,
                         type.rational.approx = c("chebfun",
                         "brasil", "chebfunLB"),
                         shared_lib = "INLA") {
    if(!inherits(graph_obj, "metric_graph")){
      stop("The graph object should be of class metric_graph!")
    }

  prior.theta.param <- prior.theta.param[[1]]

  if(!(prior.theta.param %in% c("theta", "spde"))){
    stop("theta.theta.param should be either 'theta' or 'spde'!")
  }

  type.rational.approx <- type.rational.approx[[1]]

  parameterization <- parameterization[[1]]

  prior.nu.dist <- prior.nu.dist[[1]]
  if (!prior.nu.dist %in% c("beta", "lognormal")) {
    stop("prior.nu.dist should be either beta or lognormal!")
  }

  if (!parameterization %in% c("matern", "spde")) {
    stop("parameterization should be either matern or spde!")
  }

  if (!type.rational.approx %in% c("chebfun", "brasil", "chebfunLB")) {
    stop("type.rational.approx should be either chebfun, brasil or chebfunLB!")
  }
    if(is.null(graph_obj$mesh)){
      if(is.null(h)){
        graph_obj$build_mesh(h = 0.01)
      } else{
        graph_obj$build_mesh(h = h)
      }
    }
    
    if(is.null(graph_obj$mesh$C)){
      graph_obj$compute_fem()
    }

  fixed_nu <- !is.null(nu)
    if (fixed_nu) {
    nu_order <- nu
    start.nu <- nu
  } else {
    nu_order <- nu.upper.bound
  }

  

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

    param <- get_parameters_rSPDE_graph(graph_obj, 2 * beta, 
            B.tau, 
            B.kappa, 
            B.sigma,
            B.range, 
            start.nu,
            start.nu + 1/2,
            parameterization,
            prior.std.dev.nominal, 
            prior.range.nominal, 
            prior.tau.mean, 
            prior.kappa.mean, 
            theta.prior.mean, 
            theta.prior.prec) 

    if(is.null(start.theta)){
      start.theta <- param$theta.prior.mean
    }

    theta.prior.mean <- param$theta.prior.mean
    theta.prior.prec <- param$theta.prior.prec
    
    mesh <- list(d = 1, C = graph_obj$mesh$C, 
                                G = graph_obj$mesh$G)
    class(mesh) <- "metric_graph"
    if(parameterization == "matern"){
        rspde_model <- rspde.matern(mesh = mesh,
                                nu.upper.bound = nu.upper.bound,
                                rspde.order = rspde.order,
                                nu = nu,
                                debug = debug,
                                B.sigma = B.sigma,
                                B.range = B.range,
                                start.theta = start.theta,
                                theta.prior.mean = theta.prior.mean,
                                theta.prior.prec = theta.prior.prec,
                                parameterization = parameterization,
                                prior.nu.dist = prior.nu.dist,
                                nu.prec.inc = nu.prec.inc,
                                type.rational.approx = type.rational.approx,
                                vec_param = param,
                                prior.theta.param = prior.theta.param
                                )
    } else{
        rspde_model <- rspde.matern(mesh = mesh,
                                nu.upper.bound = nu.upper.bound,
                                rspde.order = rspde.order,
                                nu = nu,
                                debug = debug,
                                B.tau = B.tau,
                                B.kappa = B.kappa,
                                start.theta = start.theta,
                                theta.prior.mean = theta.prior.mean,
                                theta.prior.prec = theta.prior.prec,
                                parameterization = parameterization,
                                prior.nu.dist = prior.nu.dist,
                                nu.prec.inc = nu.prec.inc,
                                type.rational.approx = type.rational.approx,
                                vec_param = param,
                                prior.theta.param = prior.theta.param
                                )
    }

        
        rspde_model$mesh <- graph_obj$clone()
        # rspde_model$n.spde <- nrow(graph_obj$mesh$E)
        rspde_model$n.spde <- nrow(graph_obj$mesh$VtE)

  class(rspde_model) <- c("rspde_metric_graph", class(rspde_model))
  return(rspde_model)
                         }

#' @name precision.inla_rspde
#' @title Get the precision matrix of `inla_rspde` objects
#' @description Function to get the precision matrix of an `inla_rspde` object created with the `rspde.matern()` function.
#' @param object The `inla_rspde` object obtained with the `rspde.matern()` function.
#' @param theta If null, the starting values for theta will be used. Otherwise, it must be suplied as a vector. 
#' For stationary models, we have `theta = c(log(tau), log(kappa), nu)`. For nonstationary models, we have
#' `theta = c(theta_1, theta_2, ..., theta_n, nu)`.
#' @param ... Currently not used.
#' @return The precision matrix.
#' @method precision inla_rspde
#' @seealso [precision.CBrSPDEobj()], [matern.operators()]
#' @export
#'
precision.inla_rspde <- function(object,
                                 theta = NULL,
                                 ...) {
  mesh_model <- object$mesh

  rspde_order <- object$rspde.order

  if(is.null(theta)){
    theta <- object$start.theta
    nu <- object$start.nu
  } else{
    n_tmp <- length(theta)
    nu <- theta[n_tmp]
    theta <- theta[-n_tmp]
  }

  if(!object$integer.nu){
    nu <- nu + 1e-10
  }

  alpha <- nu + object$dim/2

  if(object$stationary){
    if(object$parameterization == "spde"){
      tau <- exp(theta[1])
      kappa <- exp(theta[2])
      op <- matern.operators(mesh = mesh_model, 
                         alpha = alpha, kappa = kappa, 
                         tau = tau, 
                         m = rspde_order,
                         parameterization = "spde",
                         type = "covariance",
                         type_rational_approximation = object$type.rational.approx)
  } else{
      sigma <- exp(theta[1])
      range <- exp(theta[2])
      op <- matern.operators(mesh = mesh_model, 
                         nu = nu, range = range, 
                         sigma = sigma, 
                         m = rspde_order,
                         type = "covariance",
                         parameterization = "matern",
                         type_rational_approximation = object$type.rational.approx)
  }
  } else {
    B_tau_vec <- object$f$cgeneric$data$matrices$B_tau
    B_kappa_vec <- object$f$cgeneric$data$matrices$B_kappa
    n_total <- length(B_tau_vec)
    dim_B_matrices <- B_tau_vec[1:2]
    B_tau <- matrix(B_tau_vec[3:n_total], dim_B_matrices[1], dim_B_matrices[2], byrow = TRUE)
    B_kappa <- matrix(B_kappa_vec[3:n_total], dim_B_matrices[1], dim_B_matrices[2], byrow = TRUE)

    op <- spde.matern.operators(B.tau = B_tau, B.kappa = B_kappa, theta = theta, alpha = alpha, parameterization = "spde",
                                  mesh = mesh_model, m = rspde_order, type = "covariance",
                                   type_rational_approximation = object$type.rational.approx)
  }

  Q <- op$Q
  return(Q)
}