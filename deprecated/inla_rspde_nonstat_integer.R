
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
#' @param debug INLA debug argument
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
#' @param start.ltau Starting value for log of tau.
#' @param start.lrange Starting value for log of range. Will not be used if start.lkappa is non-null.
#' @param start.lstd.dev Starting value for log of std. deviation. Will not be used if start.ltau is non-null.
#' @param parameterization Which parameterization to use? `matern` uses range, std. deviation and nu (smoothness). `spde` uses kappa, tau and nu (smoothness). The default is `matern`.
#' @param prior.nu.dist The distribution of the smoothness parameter.
#' The current options are "beta" or "lognormal". The default is "beta".
#' @param nu.prec.inc Amount to increase the precision in the beta prior
#' distribution. Check details below.
#' @param type.rational.approx Which type of rational approximation
#' should be used? The current types are "chebfun", "brasil" or "chebfunLB".
#' @param shared_lib Which shared lib to use for the cgeneric implementation? 
#' If "INLA", it will use the shared lib from INLA's installation. If 'rSPDE', then
#' it will use the local installation (does not work if your installation is from CRAN).
#' Otherwise, you can directly supply the path of the .so (or .dll) file.
#' 
#' @return An INLA model.
#' @export

rspde.matern4 <- function(mesh,
                         nu.upper.bound = 4, rspde.order = 2,
                         nu = NULL, 
                         debug = FALSE,
                         prior.nu = NULL,
                         start.nu = NULL,
                         start.theta = NULL,
                         B_tau = c(0,1,0),
                         B_kappa = c(0,0,1),
                         prior.variance.nominal = 1,
                         prior.range.nominal = NULL,
                         prior.tau = NULL,
                         prior.kappa = NULL,
                         theta.prior.mean = NULL,
                         theta.prior.prec = 0.1, 
                         parameterization = c("matern", "spde"),
                         prior.nu.dist = c("lognormal", "beta"),
                         nu.prec.inc = 1,
                         type.rational.approx = c("chebfun",
                         "brasil", "chebfunLB"),
                         shared_lib = "INLA") {
  type.rational.approx <- type.rational.approx[[1]]

  parameterization <- parameterization[[1]]

  ### Handle parameterization = "matern" in R
  # Create a new B.tau and B.kappa

        #   nu = alpha - d/2
        # kappa0 = log(8 * nu)/2
        # tau0 = 0.5 * (lgamma(nu) - lgamma(nu + d/2) - d/2 * log(4 * 
        #     pi)) - nu * kappa0
        #  B.tau = cbind(tau0, 
        #     nu, -1)
        #     B.kappa = cbind(kappa0, -1, 0)

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


    if(inherits(mesh, c("inla.mesh", "inla.mesh.1d"))){
      if (integer_alpha) {
        integer.nu <- TRUE
        if (d == 1) {
          fem_mesh <- fem_mesh_order_1d(mesh, m_order = m_alpha + 1)
        } else {
          fem_mesh <- INLA::inla.mesh.fem(mesh, order = m_alpha)
        }
      } else {
        if (d == 1) {
          fem_mesh <- fem_mesh_order_1d(mesh, m_order = m_alpha + 2)
        } else {
          fem_mesh <- INLA::inla.mesh.fem(mesh, order = m_alpha + 1)
        }
      }
    } else{
      if(is.null(mesh$C) || is.null(mesh$G)){
        stop("If mesh is not an inla.mesh object, you should manually supply a list with elements c0, g1, g2...")
      }
      fem_mesh <- generic_fem_mesh_order(mesh, m_order = m_alpha + 2)
    }

    C <- fem_mesh[["c0"]]
    G <- fem_mesh[["g1"]]

    n_cgeneric <- ncol(fem_mesh[["c0"]])

    B_tau <- cbind(0,rep(1,n_cgeneric),0)
    B_kappa <- cbind(0, 0, rep(1,n_cgeneric))

    if(!is.null(start.theta)){
      if(length(start.theta) != ncol(B_kappa) - 1){
        stop("The length of starting values for theta is incorrect!")
      }
    }

    param <- INLA::param2.matern.orig(mesh, 2*beta, B_tau, B_kappa, 
            prior.variance.nominal, prior.range.nominal, prior.tau, 
            prior.kappa, theta.prior.mean, theta.prior.prec)

    if(is.null(start.theta)){
      start.theta <- param$theta.prior.mean
    }

  
  if(shared_lib == "INLA"){
    rspde_lib <- dirname(INLA:::inla.call.builtin())
		rspde_lib <- paste0(rspde_lib, "/external/rSPDE/librSPDE.so")
  } else if(shared_lib == "rSPDE"){
    rspde_lib <- system.file('src', package='rSPDE')
    if(Sys.info()['sysname']=='Windows') {
		rspde_lib <- paste0(rspde_lib, "/rspde_cgeneric_models.dll")
            } else {
		rspde_lib <- paste0(rspde_lib, "/rSPDE.so")
            }
  }

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
            B_tau = B_tau,
            B_kappa = B_kappa,
            start.theta = start.theta,
            theta.prior.mean = param$theta.prior.mean,
            theta.prior.prec = param$theta.prior.prec
            ))
    
    model$cgeneric_type <- "general"
 
  model$nu <- nu
  model$prior.nu <- prior.nu
  model$start.nu <- start.nu
  model$rspde.order <- rspde.order
  class(model) <- c("inla_rspde", class(model))
  model$dim <- d
  model$n.spde <- mesh$n
  model$nu.upper.bound <- nu.upper.bound
  model$debug <- debug
  model$type.rational.approx <- type.rational.approx
  model$mesh <- mesh
  model$fem_mesh <- fem_mesh
  model$param <- param
  return(model)
}
