utils::globalVariables(c(
  "C", "C_inv", "C_inv_G", "G", "d", "loc", "n",
  "n_m", "nu", "nu_upper_bound",
  "do_optimize", "idx_symmetric", "n_Q", "rspde_order",
  "graph_opt", "fem_matrices", "sharp",
  "prior.kappa", "prior.nu", "prior.tau",
  "start.lkappa", "start.ltau", "start.nu",
  "prior.nu.dist", "type.rational.approx"
))

#' @importFrom stats dnorm pnorm dbeta
#' @importFrom methods as
#' @name 'inla.rgeneric.cov_rspde_general'
#' @title Generic INLA method for the covariance-based rSPDE approach
#' @description Generic INLA method for the covariance-based rSPDE approach
#' @param cmd INLA parameter
#' @param theta Parameters of the model
#' @param ... Additional arguments to the call, or arguments with changed
#' values. Use name = NULL to remove the argument name.
#' @param args Arguments.
#' @return An INLA model
#' @noRd
"inla.rgeneric.cov_rspde_general" <- function(cmd = c(
                                                "graph",
                                                "Q",
                                                "mu",
                                                "initial",
                                                "log.norm.const",
                                                "log.prior",
                                                "quit"
                                              ),
                                              theta = NULL,
                                              args = NULL,
                                              ...) {
  initial <- function(n, theta) {
    return(c(start.ltau, start.lkappa, log(start.nu /
    (nu_upper_bound - start.nu))))
  }

  ######## parameter
  interpret.theta <- function(n, theta) {
    return(list(
      ltau = theta[1L],
      lkappa = theta[2L],
      lnu = theta[3L]
    ))
  }

  ######## precision matrix
  Q <- function(n, theta) {
    param <- interpret.theta(n, theta)

    nu <- (exp(param$lnu) / (1 + exp(param$lnu))) * nu_upper_bound
    tau <- exp(param$ltau)
    kappa <- exp(param$lkappa)

    if (do_optimize) {
      return(rSPDE::rspde.matern.precision.opt(
        kappa = kappa, nu = nu, tau = tau,
        rspde_order = rspde_order,
        d = d, fem_matrices = fem_matrices,
        sharp = sharp, graph = NULL,
        type_rational_approx = type.rational.approx
      ))
    } else {
      return(rSPDE::rspde.matern.precision(
        kappa = kappa, nu = nu, tau = tau,
        rspde_order = rspde_order,
        d = d, fem_mesh_matrices = fem_matrices,
        type_rational_approx = type.rational.approx
      ))
    }
  }


  ############################# mean
  mu <- function(n, theta) {
    return(numeric(0))
  }
  ###################### log normal constant
  log.norm.const <- function(n, theta) {
    return(numeric(0))
  }

  ############################# graph skeleton
  graph <- function(n, theta) {
    return(graph_opt)
  }



  ######################## log prior
  log.prior <- function(n, theta) {
    param <- interpret.theta(n, theta)

    if (prior.nu.dist == "lognormal") {
      nu <- (exp(param$lnu) / (1 + exp(param$lnu))) * nu_upper_bound

      tdnorm_nu <- dnorm(log(nu), 0, 1, log = TRUE) - log(nu) -
        pnorm(log(nu_upper_bound), prior.nu$loglocation,
          prior.nu$logscale,
          log.p = TRUE
        )

      res <- tdnorm_nu + dnorm(param$lkappa, prior.kappa$meanlog,
        prior.kappa$sdlog,
        log = TRUE
      ) +
        dnorm(param$ltau, prior.tau$meanlog, prior.tau$sdlog, log = TRUE) -
        param$ltau - param$lkappa
    } else if (prior.nu.dist == "beta") {
      s_1 <- (prior.nu[["mean"]] / nu_upper_bound) * prior.nu$prec
      s_2 <- (1 - prior.nu[["mean"]] / nu_upper_bound) * prior.nu$prec

      nu <- (exp(param$lnu) / (1 + exp(param$lnu))) * nu_upper_bound

      res <- dbeta(nu / nu_upper_bound, shape1 = s_1, shape2 = s_2,
      log = TRUE) - log(nu_upper_bound) +
        dnorm(param$lkappa, prior.kappa$meanlog, prior.kappa$sdlog,
        log = TRUE) +
        dnorm(param$ltau, prior.tau$meanlog, prior.tau$sdlog,
        log = TRUE) -
        param$lkappa - param$ltau
    } else {
      stop("You must choose prior.nu.dist between beta and lognormal!")
    }

    return(res)
  }


  quit <- function(n, theta) {
    return(invisible())
  }

  if (!length(theta)) {
    theta <- initial(n, theta)
  }
  res <-
    do.call(match.arg(cmd), args = list(n = as.integer(args$n), theta = theta))
  return(res)
}

#' @name 'inla.rgeneric.cov_rspde_frac_alpha'
#' @title Non-integer fixed smoothmess generic INLA method
#' @description Non-integer fixed smoothmess generic INLA method
#' for the covariance-based rSPDE approach
#' @param cmd INLA parameter
#' @param theta Parameters of the model
#' @param ... Additional arguments to the call, or arguments with
#' changed values. Use name = NULL to remove the argument name.
#' @param args Arguments.
#' @return An INLA model
#' @noRd
"inla.rgeneric.cov_rspde_frac_alpha" <- function(cmd = c(
                                                   "graph",
                                                   "Q",
                                                   "mu",
                                                   "initial",
                                                   "log.norm.const",
                                                   "log.prior",
                                                   "quit"
                                                 ),
                                                 theta = NULL,
                                                 args = NULL,
                                                 ...) {
  initial <- function(n, theta) {
    return(c(start.ltau, start.lkappa))
  }

  ######## parameter
  interpret.theta <- function(n, theta) {
    return(list(
      ltau = theta[1L],
      lkappa = theta[2L]
    ))
  }

  ######## precision matrix
  Q <- function(n, theta) {
    param <- interpret.theta(n, theta)
    tau <- exp(param$ltau)
    kappa <- exp(param$lkappa)

    if (do_optimize) {
      return(rSPDE::rspde.matern.precision.opt(
        kappa = kappa, nu = nu, tau = tau,
        rspde_order = rspde_order,
        d = d, fem_matrices = fem_matrices,
        sharp = sharp, graph = NULL,
        type_rational_approx = type.rational.approx
      ))
    } else {
      return(rSPDE::rspde.matern.precision(
        kappa = kappa, nu = nu, tau = tau,
        rspde_order = rspde_order,
        d = d, fem_mesh_matrices = fem_matrices,
        type_rational_approx = type.rational.approx
      ))
    }
  }


  ############################# mean
  mu <- function(n, theta) {
    return(numeric(0))
  }
  ###################### log normal constant
  log.norm.const <- function(n, theta) {
    return(numeric(0))
  }

  ############################# graph skeleton
  graph <- function(n, theta) {
    return(graph_opt)
  }

  ######################## log prior
  log.prior <- function(n, theta) {
    param <- interpret.theta(n, theta)

    res <- dnorm(param$lkappa, prior.kappa$meanlog,
      prior.kappa$sdlog,
      log = TRUE
    ) +
      dnorm(param$ltau, prior.tau$meanlog,
        prior.tau$sdlog,
        log = TRUE
      ) -
      param$lkappa - param$ltau
    return(res)
  }


  quit <- function(n, theta) {
    return(invisible())
  }

  if (!length(theta)) {
    theta <- initial(n, theta)
  }
  res <-
    do.call(match.arg(cmd), args = list(n = as.integer(args$n), theta = theta))
  return(res)
}


#' @name 'inla.rgeneric.cov_rspde_int_alpha'
#' @title Integer fixed smoothmess generic INLA method
#' @description Integer fixed smoothmess generic INLA method
#' for the covariance-based rSPDE approach
#' @param cmd INLA parameter
#' @param theta Parameters of the model
#' @param ... Additional arguments to the call, or arguments with
#' changed values. Use name = NULL to remove the argument name.
#' @param args Arguments.
#' @return An INLA model
#' @noRd
"inla.rgeneric.cov_rspde_int_alpha" <- function(cmd = c(
                                                  "graph",
                                                  "Q",
                                                  "mu",
                                                  "initial",
                                                  "log.norm.const",
                                                  "log.prior",
                                                  "quit"
                                                ),
                                                theta = NULL,
                                                args = NULL,
                                                ...) {
  initial <- function(n, theta) {
    return(c(start.ltau, start.lkappa))
  }

  ######## parameter
  interpret.theta <- function(n, theta) {
    return(list(
      ltau = theta[1L],
      lkappa = theta[2L]
    ))
  }

  ######## precision matrix
  Q <- function(n, theta) {
    param <- interpret.theta(n, theta)
    tau <- exp(param$ltau)
    kappa <- exp(param$lkappa)



    if (do_optimize) {
      return(rSPDE::rspde.matern.precision.integer.opt(kappa, nu, tau,
        d, fem_matrices,
        graph = NULL
      ))
    } else {
      return(rSPDE::rspde.matern.precision.integer(
        kappa = kappa, nu = nu, tau = tau,
        d = d, fem_mesh_matrices = fem_matrices
      ))
    }
  }


  ############################# mean
  mu <- function(n, theta) {
    return(numeric(0))
  }
  ###################### log normal constant
  log.norm.const <- function(n, theta) {
    return(numeric(0))
  }

  ############################# graph skeleton
  graph <- function(n, theta) {
    return(graph_opt)
  }



  ######################## log prior
  log.prior <- function(n, theta) {
    param <- interpret.theta(n, theta)

    res <- dnorm(param$lkappa, prior.kappa$meanlog,
      prior.kappa$sdlog,
      log = TRUE
    ) +
      dnorm(param$ltau, prior.tau$meanlog,
        prior.tau$sdlog,
        log = TRUE
      ) -
      param$lkappa - param$ltau
    return(res)
  }


  quit <- function(n, theta) {
    return(invisible())
  }

  if (!length(theta)) {
    theta <- initial(n, theta)
  }
  res <-
    do.call(match.arg(cmd), args = list(n = as.integer(args$n), theta = theta))
  return(res)
}



#' @name rspde.matern
#' @title Matern rSPDE model object for INLA
#' @description Creates an INLA object for a stationary Matern model with
#' general smoothness parameter.
#' @param mesh The mesh to build the model. Should be an \code{inla.mesh} or
#' an \code{inla.mesh.1d} object.
#' @param nu_upper_bound Upper bound for the smoothness parameter.
#' @param rspde_order The order of the covariance-based rational SPDE approach.
#' @param nu If nu is set to a parameter, nu will be kept fixed and will not
#' be estimated. If nu is \code{NULL}, it will be estimated.
#' @param sharp The sparsity graph should have the correct sparsity (costs
#' more to perform a sparsity analysis) or an upper bound for the sparsity? If
#' \code{TRUE}, the graph will have the correct sparsity.
#' @param debug INLA debug argument.
#' @param optimize Should the model be optimized? In this case the sparsities of
#' the matrices will be analyzed.
#' @param prior.kappa a \code{list} containing the elements \code{meanlog} and
#' \code{sdlog}, that is, the mean and standard deviation on the log scale.
#' @param prior.nu a list containing the elements \code{mean} and \code{prec}
#' for beta distribution, or \code{loglocation} and \code{logscale} for a
#' truncated lognormal distribution. \code{loglocation} stands for
#' the location parameter of the truncated lognormal distribution in the log
#' scale. \code{prec} stands for the precision of a beta distribution.
#' \code{logscale} stands for the scale of the truncated lognormal
#' distribution on the log scale. Check details below.
#' @param prior.tau a list containing the elements \code{meanlog} and
#' \code{sdlog}, that is, the mean and standard deviation on the log scale.
#' @param start.lkappa Starting value for log of kappa.
#' @param start.nu Starting value for nu.
#' @param start.ltau Starting value for log of tau.
#' @param prior.nu.dist The distribution of the smoothness parameter.
#' The current options are "beta" or "lognormal". The default is "beta".
#' @param nu.prec.inc Amount to increase the precision in the beta prior
#' distribution. Check details below.
#' @param type.rational.approx Which type of rational approximation
#' should be used? The current types are "chebfun", "brasil" or "chebfunLB".
#'
#' @return An INLA model.
#' @export
#' @details This function constructs a stationary Matern rSPDE model to
#' be used with the INLA interface. The parameters are the range
#' parameter \eqn{\kappa}, the smoothness parameter
#' \eqn{\nu} and the variance rescaling parameter \eqn{\tau}.
#'
#' For this model, an upper bound for the smoothness parameter
#' \eqn{\nu} should be given. It is given by the
#' \code{nu_upper_bound} argument.
#'
#' It is very important to notice that the larger the value
#' of \code{nu_upper_bound} the higher the computational cost
#' to fit the model. So, it is generally best to initially fit
#' a model with a small value of \code{nu_upper_bound} and
#' increase it only if it is really needed (for instance, if
#' the estimated smoothness parameter was very close to
#' \code{nu_upper_bound}).
#'
#' The following parameterization is used:
#' \deqn{\log(\tau) = \theta_1,}
#' \deqn{\log(\kappa) = \theta_2}
#' and for \eqn{\theta_3} we can have a beta prior
#' or a truncated lognormal prior distribution. In each case,
#' the prior distribution has support on the interval
#' \eqn{(0,\nu_{UB})}, where \eqn{\nu_{UB}} is
#' \code{nu_upper_bound}. Then, the following parameterization
#' is considered:
#' \deqn{\log\Big(\frac{\nu}{\nu_{UB}-\nu}\Big) = \theta_3.}
#'
#' By default, an optimized version of this model is considered. The optimized
#' version is generally much faster for larger datasets, however it takes more
#' time to build the model as the sparsity of the graph should be analyzed.
#' However, for small datasets, it is possible that the time taken to
#' analyze sparsity plus fitting the model is larger than the time taken
#' to fit an unoptimized model. So, for a small dataset it might be
#' convenient to set \code{optimize=FALSE}.
#'
#' A way to use the optimized version but reduce the cost of sparsity analysis
#' is to set \code{sharp} to \code{FALSE}. However, it should increase
#' the cost of fitting the model. Therefore, one
#' usually would not benefit from setting the \code{sharp} argument to
#' \code{FALSE} when fitting the model to large datasets.
#'
#' Finally, when considering a beta prior, the beta distribution will be
#' parameterized in terms of its mean, say \eqn{\mu} and a precision
#' parameter \eqn{\phi}, which is such that the variance of the beta
#' distribution is given by \eqn{\mu(\nu_{UB}-\mu)/(1+\phi)}.
#' The mean of the beta prior is determined by the \code{prior.nu$mean}, whereas
#' the precision parameter is determined by the \code{prior.nu$prec}. If
#' \code{prior.nu$prec} is \code{NULL} (which is the default case), the
#' precision parameter is taken
#' as
#' \deqn{\phi = \max\Big\{\frac{\nu_{UB}}{\mu},
#' \frac{\nu_{UB}}{\nu_{UB}-\mu}\Big\} + \textrm{nu.prec.inc},}
#' where \eqn{\mu} is the prior mean of the smoothness parameter.
#'
#' This choice of precision parameter is to ensure that the prior beta density
#' has boundary values equal to zero (where the boundary values are defined
#' either by continuity or by limits).
#'
#' Hence, the higher the value of \code{nu.prec.inc} the more informative
#' the prior is.
#'
#' @examples
#' \donttest{
#' # devel version
#' library(INLA)
#'
#' # Organizing the data
#' data(PRprec)
#' data(PRborder)
#'
#' Y <- rowMeans(PRprec[, 3 + 1:31])
#' ind <- !is.na(Y)
#' Y <- Y[ind]
#' coords <- as.matrix(PRprec[ind, 1:2])
#' alt <- PRprec$Altitude[ind]
#'
#' seaDist <- apply(
#'   spDists(coords, PRborder[1034:1078, ], longlat = TRUE),
#'   1, min
#' )
#'
#' # Creating INLA mesh
#' prdomain <- inla.nonconvex.hull(coords, -0.03, -0.05,
#' resolution = c(80, 80))
#' prmesh <- inla.mesh.2d(boundary = prdomain,
#' max.edge = c(0.6, 1.2), cutoff = 0.3)
#'
#' # Building the A matrix
#' Abar <- rspde.make.A(mesh = prmesh, loc = coords)
#'
#' # Building the index
#' mesh.index <- rspde.make.index(name = "field", mesh = prmesh)
#'
#' # Creating the model
#' rspde_model <- rspde.matern(mesh = prmesh)
#'
#' # INLA stack
#' stk.dat <- inla.stack(
#'   data = list(y = Y), A = list(Abar, 1), tag = "est",
#'   effects = list(
#'     c(
#'       mesh.index,
#'       list(Intercept = 1)
#'     ),
#'     list(
#'       long = inla.group(coords[, 1]),
#'       lat = inla.group(coords[, 2]),
#'       seaDist = inla.group(seaDist)
#'     )
#'   )
#' )
#'
#' # INLA formula
#' f.s <- y ~ -1 + Intercept + f(seaDist, model = "rw1") +
#'   f(field, model = rspde_model)
#'
#' # Fitting the model
#' rspde_fit <- inla(f.s,
#'   family = "Gamma", data = inla.stack.data(stk.dat),
#'   control.inla = list(int.strategy = "eb"),
#'   control.predictor = list(A = inla.stack.A(stk.dat)),
#'            inla.mode = "experimental"
#' )
#'
#' # The result
#' summary(rspde_fit)
#' # devel.tag
#' }
#'
rspde.matern <- function(mesh,
                         nu_upper_bound = 4, rspde_order = 2,
                         nu = NULL, sharp = TRUE,
                         debug = FALSE,
                         optimize = TRUE,
                         prior.kappa = NULL,
                         prior.nu = NULL,
                         prior.tau = NULL,
                         start.lkappa = NULL,
                         start.nu = NULL,
                         start.ltau = NULL,
                         prior.nu.dist = c("beta", "lognormal"),
                         nu.prec.inc = 1,
                         type.rational.approx = c("chebfun",
                         "brasil", "chebfunLB")) {
  type.rational.approx <- type.rational.approx[[1]]

  prior.nu.dist <- prior.nu.dist[[1]]
  if (!prior.nu.dist %in% c("beta", "lognormal")) {
    stop("prior.nu.dist should be either beta or lognormal!")
  }

  integer.nu <- FALSE

  if (mesh$manifold == "R1") {
    d <- 1
  } else if (mesh$manifold == "R2") {
    d <- 2
  } else {
    stop("The domain must be flat manifolds of dimension 1 or 2, that is,
         the domain must be a line or a plane.")
  }

  if (nu_upper_bound - floor(nu_upper_bound) == 0) {
    nu_upper_bound <- nu_upper_bound - 1e-5
  }

  if (!is.null(nu)) {
    if (!is.numeric(nu)) {
      stop("nu must be numeric!")
    }
  }

  if (d == 1) {
    if (nu_upper_bound > 2) {
      warning("In dimension 1 you can have unstable results
      for nu_upper_bound > 2. Consider changing
      nu_upper_bound to 2 or 1.")
    }
  }

  fixed_nu <- !is.null(nu)

  if (fixed_nu) {
    alpha <- nu + d / 2
    integer_alpha <- (alpha %% 1 == 0)
  } else {
    integer_alpha <- FALSE
  }

  if (fixed_nu) {
    nu_order <- nu
  } else {
    nu_order <- nu_upper_bound
  }

  if (optimize) {
    beta <- nu_order / 2 + d / 4

    m_alpha <- floor(2 * beta)

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
  } else {
    fem_mesh <- NULL
    beta <- nu_order / 2 + d / 4

    m_alpha <- floor(2 * beta)

    if (integer_alpha) {
      integer.nu <- TRUE
      if (d == 1) {
        fem_matrices <- fem_mesh_order_1d(mesh, m_order = m_alpha + 1)
      } else {
        fem_matrices <- INLA::inla.mesh.fem(mesh, order = m_alpha)
      }
    } else {
      if (d == 1) {
        fem_matrices <- fem_mesh_order_1d(mesh, m_order = m_alpha + 2)
      } else {
        fem_matrices <- INLA::inla.mesh.fem(mesh, order = m_alpha + 1)
      }
    }
  }

  if (optimize) {
    if (integer_alpha) {
      result_sparsity <- analyze_sparsity_rspde(
        nu_upper_bound = nu_order, dim = d,
        rspde_order = rspde_order,
        fem_mesh_matrices = fem_mesh,
        include_higher_order = FALSE
      )
    } else {
      if (sharp) {
        result_sparsity <- analyze_sparsity_rspde(
          nu_upper_bound = nu_order, dim = d,
          rspde_order = rspde_order,
          fem_mesh_matrices = fem_mesh
        )
      } else {
        result_sparsity <- analyze_sparsity_rspde(
          nu_upper_bound = nu_order, dim = d,
          rspde_order = rspde_order,
          fem_mesh_matrices = fem_mesh,
          include_lower_order = FALSE
        )
      }
      positions_matrices <- result_sparsity$positions_matrices
    }

    idx_matrices <- result_sparsity$idx_matrices
    positions_matrices_less <- result_sparsity$positions_matrices_less
  } else {
    positions_matrices <- NULL
    idx_matrices <- NULL
    positions_matrices_less <- NULL
  }

  if (optimize) {
    fem_matrices <- list()

    if (sharp || integer_alpha) {
      if (m_alpha > 0) {
        fem_matrices[[paste0("G_", m_alpha, "_less")]] <-
        fem_mesh[[paste0("g", m_alpha)]]@x[idx_matrices[[m_alpha + 1]]]

        fem_matrices[["C_less"]] <-
        rep(0, length(fem_matrices[[paste0("G_", m_alpha, "_less")]]))
        fem_matrices[["C_less"]][positions_matrices_less[[1]]] <-
        fem_mesh$c0@x[idx_matrices[[1]]]

        fem_matrices[["G_less"]] <-
        rep(0, length(fem_matrices[[paste0("G_", m_alpha, "_less")]]))
        fem_matrices[["G_less"]][positions_matrices_less[[2]]] <-
        fem_mesh$g1@x[idx_matrices[[2]]]
      } else {
        fem_matrices[["C_less"]] <- fem_mesh[["c0"]]@x
        fem_matrices[["G_less"]] <- fem_mesh[["g1"]]@x
      }


      # The case m_alpha=2 already uses G_2_less defined above
      if (m_alpha > 2) {
        for (j in 2:(m_alpha - 1)) {
          fem_matrices[[paste0("G_", j, "_less")]] <-
          rep(0, length(fem_matrices[[paste0("G_", m_alpha, "_less")]]))
          fem_matrices[[paste0("G_",
          j, "_less")]][positions_matrices_less[[j + 1]]] <-
          fem_mesh[[paste0("g", j)]]@x[idx_matrices[[j + 1]]]
        }
      }
    }


    if (!integer_alpha) {
      if (m_alpha == 0) {
        fem_matrices[["G"]] <- fem_mesh$g1@x
        fem_matrices[["C"]] <- fem_mesh$c0@x
      } else if (m_alpha > 0) {
        fem_matrices[[paste0("G_", m_alpha + 1)]] <-
        fem_mesh[[paste0("g", m_alpha + 1)]]@x[idx_matrices[[m_alpha + 2]]]

        fem_matrices[["G"]] <-
        rep(0, length(fem_matrices[[paste0("G_", m_alpha + 1)]]))
        fem_matrices[["G"]][positions_matrices[[2]]] <-
        fem_mesh$g1@x[idx_matrices[[2]]]

        fem_matrices[["C"]] <-
        rep(0, length(fem_matrices[[paste0("G_", m_alpha + 1)]]))
        fem_matrices[["C"]][positions_matrices[[1]]] <-
        fem_mesh$c0@x[idx_matrices[[1]]]
      }
      if (m_alpha > 1) {
        for (j in 2:(m_alpha)) {
          fem_matrices[[paste0("G_", j)]] <-
          rep(0, length(fem_matrices[[paste0("G_", m_alpha + 1)]]))
          fem_matrices[[paste0("G_",
          j)]][positions_matrices[[j + 1]]] <-
          fem_mesh[[paste0("g", j)]]@x[idx_matrices[[j + 1]]]
        }
      }
    }
  }






  if (is.null(prior.nu$loglocation)) {
    prior.nu$loglocation <- log(min(1, nu_upper_bound / 2))
  }

  if (is.null(prior.nu[["mean"]])) {
    prior.nu[["mean"]] <- min(1, nu_upper_bound / 2)
  }

  if (is.null(prior.kappa$meanlog)) {
    mesh.range <- ifelse(d == 2, (max(c(diff(range(mesh$loc[
      ,
      1
    ])), diff(range(mesh$loc[, 2])), diff(range(mesh$loc[
      ,
      3
    ]))))), diff(mesh$interval))
    prior.range.nominal <- mesh.range * 0.2
    if (prior.nu.dist == "lognormal") {
      prior.kappa$meanlog <- log(sqrt(8 *
      exp(prior.nu[["loglocation"]])) / prior.range.nominal)
    } else if (prior.nu.dist == "beta") {
      prior.kappa$meanlog <- log(sqrt(8 *
      prior.nu[["mean"]]) / prior.range.nominal)
    }
  }

  if (is.null(prior.tau$meanlog)) {
    if (prior.nu.dist == "lognormal") {
      prior.tau$meanlog <- log(sqrt(gamma(exp(prior.nu[["loglocation"]])) /
      gamma(exp(prior.nu[["loglocation"]]) + d / 2) / (4 *
        pi * exp(prior.kappa$meanlog)^(2 * exp(prior.nu[["loglocation"]])))))
    } else if (prior.nu.dist == "beta") {
      prior.tau$meanlog <- log(sqrt(gamma(prior.nu[["mean"]]) /
      gamma(prior.nu[["mean"]] + d / 2) / (4 *
        pi * exp(prior.kappa$meanlog)^(2 * prior.nu[["mean"]]))))
    }
  }
  if (is.null(prior.kappa$sdlog)) {
    prior.kappa$sdlog <- sqrt(10)
  }
  if (is.null(prior.nu$prec)) {
    mu_temp <- prior.nu[["mean"]] / nu_upper_bound
    prior.nu$prec <- max(1 / mu_temp, 1 / (1 - mu_temp)) + nu.prec.inc
  }

  if (is.null(prior.nu[["logscale"]])) {
    prior.nu[["logscale"]] <- 1
  }

  if (is.null(prior.tau$sdlog)) {
    prior.tau$sdlog <- sqrt(10)
  }

  if (is.null(start.lkappa)) {
    start.lkappa <- prior.kappa$meanlog
  }
  if (is.null(start.ltau)) {
    start.ltau <- prior.tau$meanlog
  }
  if (is.null(start.nu)) {
    if (prior.nu.dist == "beta") {
      start.nu <- prior.nu[["mean"]]
    } else if (prior.nu.dist == "lognormal") {
      start.nu <- exp(prior.nu[["loglocation"]])
    } else {
      stop("prior.nu.dist should be either beta or lognormal!")
    }
  } else if (start.nu > nu_upper_bound || start.nu < 0) {
    stop("start.nu should be a number between 0 and nu_upper_bound!")
  }


  if (!fixed_nu) {
    if (optimize) {
      graph_opt <- rSPDE::get.sparsity.graph.rspde(
        fem_mesh_matrices = fem_mesh, dim = d,
        nu = nu_upper_bound,
        rspde_order = rspde_order,
        sharp = sharp,
        force_non_integer = TRUE
      )
    } else {
      graph_opt <- rSPDE::get.sparsity.graph.rspde(
        fem_mesh_matrices = fem_matrices, dim = d,
        nu = nu_upper_bound,
        rspde_order = rspde_order,
        sharp = TRUE,
        force_non_integer = TRUE
      )
    }
    model <- INLA::inla.rgeneric.define(inla.rgeneric.cov_rspde_general,
      nu_upper_bound = nu_upper_bound,
      fem_matrices = fem_matrices,
      graph_opt = graph_opt,
      sharp = sharp,
      prior.kappa = prior.kappa,
      prior.nu = prior.nu,
      prior.tau = prior.tau,
      start.lkappa = start.lkappa,
      start.nu = start.nu,
      start.ltau = start.ltau,
      type.rational.approx = type.rational.approx,
      d = d, rspde_order = rspde_order,
      prior.nu.dist = prior.nu.dist,
      n = ncol(C) * (rspde_order + 1),
      debug = debug,
      do_optimize = optimize, optimize = optimize
    )
    model$rgeneric_type <- "general"
    model$optimize <- optimize
  } else if (!integer_alpha) {
    if (optimize) {
      graph_opt <- rSPDE::get.sparsity.graph.rspde(
        fem_mesh_matrices = fem_mesh, dim = d,
        nu = nu,
        rspde_order = rspde_order,
        sharp = sharp
      )
    } else {
      graph_opt <- rSPDE::get.sparsity.graph.rspde(
        fem_mesh_matrices = fem_matrices, dim = d,
        nu = nu,
        rspde_order = rspde_order,
        sharp = TRUE, force_non_integer = TRUE
      )
    }

    model <- INLA::inla.rgeneric.define(inla.rgeneric.cov_rspde_frac_alpha,
      nu = nu,
      fem_matrices = fem_matrices,
      graph_opt = graph_opt,
      sharp = sharp,
      prior.kappa = prior.kappa,
      prior.nu = prior.nu,
      prior.tau = prior.tau,
      start.lkappa = start.lkappa,
      start.ltau = start.ltau,
      type.rational.approx = type.rational.approx,
      d = d, rspde_order = rspde_order,
      n = ncol(C) * (rspde_order + 1),
      debug = debug,
      do_optimize = optimize, optimize = optimize
    )
    model$rgeneric_type <- "frac_alpha"
    model$optimize <- optimize
  } else {
    if (optimize) {
      graph_opt <- rSPDE::get.sparsity.graph.rspde(
        fem_mesh_matrices = fem_mesh, dim = d,
        nu = nu,
        rspde_order = rspde_order,
        sharp = sharp
      )
    } else {
      graph_opt <- rSPDE::get.sparsity.graph.rspde(
        fem_mesh_matrices = fem_matrices, dim = d,
        nu = nu,
        rspde_order = rspde_order,
        force_non_integer = FALSE
      )
    }
    model <- INLA::inla.rgeneric.define(inla.rgeneric.cov_rspde_int_alpha,
      nu = nu,
      fem_matrices = fem_matrices,
      graph_opt = graph_opt,
      prior.kappa = prior.kappa,
      prior.nu = prior.nu,
      prior.tau = prior.tau,
      type.rational.approx = type.rational.approx,
      start.lkappa = start.lkappa,
      start.ltau = start.ltau,
      d = d,
      n = ncol(C),
      debug = debug,
      do_optimize = optimize, optimize = optimize
    )
    model$rgeneric_type <- "int_alpha"
    model$optimize <- optimize
  }

  model$nu <- nu
  model$prior.kappa <- prior.kappa
  model$prior.nu <- prior.nu
  model$prior.tau <- prior.tau
  model$start.lkappa <- start.lkappa
  model$start.ltau <- start.ltau
  model$start.nu <- start.nu
  model$integer.nu <- integer.nu
  if (integer.nu) {
    rspde_order <- 0
  }
  model$rspde_order <- rspde_order
  class(model) <- c(class(model), "inla.rspde")
  model$dim <- d
  model$est_nu <- !fixed_nu
  model$n.spde <- mesh$n
  model$nu_upper_bound <- nu_upper_bound
  model$prior.nu.dist <- prior.nu.dist
  model$sharp <- sharp
  model$debug <- debug
  model$type.rational.approx <- type.rational.approx
  model$mesh <- mesh
  model$fem_mesh <- fem_mesh
  model$fem_matrices <- fem_matrices

  return(model)
}




#' @name rspde.matern.precision.opt
#' @title Optimized precision matrix of the covariance-based rational
#' approximation
#' @description \code{rspde.matern.precision} is used for computing the
#' optimized version of the precision matrix of the
#' covariance-based rational SPDE approximation of a stationary Gaussian random
#' fields on \eqn{R^d} with a Matern covariance function
#' \deqn{C(h) = \frac{\sigma^2}{2^{\nu-1}\Gamma(\nu)}(\kappa h)^\nu
#' K_\nu(\kappa h).}
#' @param kappa Range parameter of the covariance function.
#' @param tau Scale parameter of the covariance function.
#' @param nu Shape parameter of the covariance function.
#' @param rspde_order The order of the rational approximation
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

rspde.matern.precision.opt <- function(kappa, nu, tau, rspde_order,
dim, fem_matrices, graph = NULL, sharp, type_rational_approx) {
  n_m <- rspde_order

  mt <- get_rational_coefficients(n_m, type_rational_approx)

  beta <- nu / 2 + dim / 4

  m_alpha <- floor(2 * beta)

  r <- sapply(1:(n_m), function(i) {
    approx(mt$alpha, mt[[paste0("r", i)]], cut_decimals(2 * beta))$y
  })

  p <- sapply(1:(n_m), function(i) {
    approx(mt$alpha, mt[[paste0("p", i)]], cut_decimals(2 * beta))$y
  })

  k <- approx(mt$alpha, mt$k, cut_decimals(2 * beta))$y


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
    graph <- as(graph, "dgTMatrix")
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
#' @description \code{rspde.matern.precision} is used for computing the
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
#' @param rspde_order The order of the rational approximation
#' @param dim The dimension of the domain
#' @param fem_mesh_matrices A list containing the FEM-related matrices. The
#' list should contain elements c0, g1, g2, g3, etc.
#' @param only_fractional Logical. Should only the fractional-order part of
#' the precision matrix be returned?
#' @param return_block_list Logical. For \code{type = "covariance"}, should the
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
#' op_cov <- matern.operators(
#'   C = fem$C, G = fem$G, nu = nu, kappa = kappa, sigma = sigma,
#'   d = 1, m = 2, compute_higher_order = TRUE
#' )
#' v <- t(rSPDE.A1d(x, 0.5))
#' c.true <- matern.covariance(abs(x - 0.5), kappa, nu, sigma)
#' Q <- rspde.matern.precision(
#'   kappa = kappa, nu = nu, tau = tau, rspde_order = 2, d = 1,
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
rspde_order, dim, fem_mesh_matrices,
only_fractional = FALSE, return_block_list = FALSE,
type_rational_approx = "chebfun") {
  if (is.null(tau) && is.null(sigma)) {
    stop("You should provide either tau or sigma!")
  }

  if (is.null(tau)) {
    tau <- sqrt(gamma(nu) / (sigma^2 * kappa^(2 * nu) *
    (4 * pi)^(dim / 2) * gamma(nu + dim / 2)))
  }

  n_m <- rspde_order

  mt <- get_rational_coefficients(n_m, type_rational_approx)

  beta <- nu / 2 + dim / 4

  m_alpha <- floor(2 * beta)

  r <- sapply(1:(n_m), function(i) {
    approx(mt$alpha, mt[[paste0("r", i)]], cut_decimals(2 * beta))$y
  })

  p <- sapply(1:(n_m), function(i) {
    approx(mt$alpha, mt[[paste0("p", i)]], cut_decimals(2 * beta))$y
  })

  k <- approx(mt$alpha, mt$k, cut_decimals(2 * beta))$y

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
      Kpart <- fem_mesh_matrices[["c0"]]
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

      Q[[length(Q) + 1]] <- kappa^(4 * beta) * tau^2 *
      fem_mesh_matrices[["c0"]] / k

      return(Q)
    } else {
      Q <- (L - p[1] * fem_mesh_matrices[["c0"]]) / r[1]

      if (n_m > 1) {
        for (i in 2:(n_m)) {
          temp <- (L - p[i] * fem_mesh_matrices[["c0"]]) / r[i]
          Q <- bdiag(Q, temp)
        }
      }

      Q <- bdiag(Q, fem_mesh_matrices[["c0"]] / k)


      Q <- Q * kappa^(4 * beta)

      Q <- tau^2 * Q
      return(Q)
    }
  }
}


#' @name rspde.matern.precision.integer.opt
#' @title Optimized precision matrix of stationary Gaussian Matern
#' random fields with integer covariance exponent
#' @description \code{rspde.matern.precision.integer.opt} is used
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
    graph <- as(graph, "dgTMatrix")
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
#' @description \code{rspde.matern.precision.integer.opt} is
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
#' op_cov <- matern.operators(
#'   C = fem$C, G = fem$G, nu = nu, kappa = kappa, sigma = sigma,
#'   d = 1, m = 2
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

#' @name rspde.make.A
#' @title Observation/prediction matrices for rSPDE models.
#' @description Constructs observation/prediction weight matrices
#' for rSPDE models based on \code{inla.mesh} or
#' \code{inla.mesh.1d} objects.
#' @param mesh An \code{inla.mesh} or
#' an \code{inla.mesh.1d} object.
#' @param loc Locations, needed if an INLA mesh is provided
#' @param A The A matrix from the standard SPDE approach, such as the matrix
#' returned by \code{inla.spde.make.A}. Should only be provided if
#' \code{mesh} is not provided.
#' @param dim the dimension. Should only be provided if an
#' \code{mesh} is not provided.
#' @param rspde_order The order of the covariance-based rational SPDE approach.
#' @param nu If \code{NULL}, then the model will assume that nu will
#' be estimated. If nu is fixed, you should provide the value of nu.
#' @param index For each observation/prediction value, an index into loc.
#' Default is \code{seq_len(nrow(A.loc))}.
#' @param group For each observation/prediction value, an index into
#' the group model.
#' @param repl For each observation/prediction value, the replicate index.
#' @param n.group The size of the group model.
#' @param n.repl The total number of replicates.
#' @return The \eqn{A} matrix for rSPDE models.
#' @export
#' @examples
#' \donttest{
#' # devel version
#' library(INLA)
#'
#' set.seed(123)
#' loc <- matrix(runif(100 * 2) * 100, 100, 2)
#' mesh <- inla.mesh.2d(
#'   loc = loc,
#'   cutoff = 50,
#'   max.edge = c(50, 500)
#' )
#' A <- rspde.make.A(mesh, loc = loc, rspde_order = 3)
#' # devel.tag
#' }
rspde.make.A <- function(mesh = NULL,
                         loc = NULL,
                         A = NULL,
                         dim = NULL,
                         rspde_order = 2, nu = NULL,
                         index = NULL,
                         group = NULL,
                         repl = 1L,
                         n.group = NULL,
                         n.repl = NULL) {
  if (!is.null(mesh)) {
    cond1 <- inherits(mesh, "inla.mesh.1d")
    cond2 <- inherits(mesh, "inla.mesh")
    stopifnot(cond1 || cond2)
    if (mesh$manifold == "R1") {
      dim <- 1
    } else if (mesh$manifold == "R2") {
      dim <- 2
    } else {
      stop("The domain must be flat manifolds of dimension 1 or 2, that is,
         the domain must be a line or a plane.")
    }
  } else if (is.null(dim)) {
    stop("If mesh is not provided, then you should provide the dimension d!")
  }

  if (!is.null(mesh)) {
    if (is.null(loc)) {
      stop("If you provided mesh, you should also provide the locations, loc.")
    }
  }

  if (!is.null(mesh)) {
    A <- INLA::inla.spde.make.A(
      mesh = mesh, loc = loc,
      index = index, group = group,
      repl = repl, n.group = n.group,
      n.repl = n.repl
    )
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
      Abar <- kronecker(matrix(1, 1, rspde_order + 1), A)
      integer_nu <- FALSE
    }
  } else {
    Abar <- kronecker(matrix(1, 1, rspde_order + 1), A)
    integer_nu <- FALSE
  }

  if (integer_nu) {
    rspde_order <- 0
  }


  attr(Abar, "inla.rspde.Amatrix") <- TRUE
  attr(Abar, "rspde_order") <- rspde_order
  attr(Abar, "integer_nu") <- integer_nu
  return(Abar)
}


#' @name rspde.make.index
#' @title rSPDE model index vector generation
#' @description Generates a list of named index vectors for an rSPDE model.
#' @param name A character string with the base name of the effect.
#' @param mesh An \code{inla.mesh} or
#' an \code{inla.mesh.1d} object.
#' @param rspde_order The order of the rational approximation
#' @param nu If \code{NULL}, then the model will assume that nu will
#' be estimated. If nu is fixed, you should provide the value of nu.
#' @param n.spde The number of basis functions in the mesh model.
#' @param n.group The size of the group model.
#' @param n.repl The total number of replicates.
#' @param dim the dimension of the domain. Should only be provided if
#' \code{mesh} is not provided.
#' @return A list of named index vectors.
#' \item{name}{Indices into the vector of latent variables}
#' \item{name.group}{'group' indices}
#' \item{name.repl}{Indices for replicates}
#' @export
#' @examples
#' \donttest{
#' # devel version
#' library(INLA)
#' set.seed(123)
#'
#' m <- 100
#' loc_2d_mesh <- matrix(runif(m * 2), m, 2)
#' mesh_2d <- inla.mesh.2d(
#'   loc = loc_2d_mesh,
#'   cutoff = 0.05,
#'   max.edge = c(0.1, 0.5)
#' )
#' sigma <- 0.01
#' range <- 0.2
#' nu <- 0.8
#' kappa <- sqrt(8 * nu) / range
#' op <- matern.operators(
#'   mesh = mesh_2d, nu = nu,
#'   kappa = kappa, sigma = sigma, m = 2
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
#'   nu_upper_bound = 1
#' )
#' f <- y ~ -1 + f(field, model = rspde_model)
#' rspde_fit <- inla(f,
#'   data = inla.stack.data(st.dat),
#'   family = "gaussian",
#'   control.predictor =
#'     list(A = inla.stack.A(st.dat)),
#'            inla.mode = "experimental"
#' )
#' result <- rspde.result(rspde_fit, "field", rspde_model)
#' plot(result)
#' # devel.tag
#' }
rspde.make.index <- function(name, n.spde = NULL, n.group = 1,
                             n.repl = 1, mesh = NULL,
                             rspde_order = 2, nu = NULL, dim = NULL) {
  if (is.null(n.spde) && is.null(mesh)) {
    stop("You should provide either n.spde or mesh!")
  }

  if (!is.null(mesh)) {
    n_mesh <- mesh$n

    if (mesh$manifold == "R1") {
      dim <- 1
    } else if (mesh$manifold == "R2") {
      dim <- 2
    } else {
      stop("The domain must be flat manifolds of dimension 1 or 2, that is,
         the domain must be a line or a plane.")
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
      factor_rspde <- rspde_order + 1
      integer_nu <- FALSE
    }
  } else {
    factor_rspde <- rspde_order + 1
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
  class(out) <- c(class(out), "inla.rspde.index")
  if (integer_nu) {
    rspde_order <- 0
  }
  attr(out, "rspde_order") <- rspde_order
  attr(out, "integer_nu") <- integer_nu
  attr(out, "n.mesh") <- n_mesh
  attr(out, "name") <- name
  attr(out, "n.group") <- n.group
  attr(out, "n.repl") <- n.repl
  return(out)
}

#' @name rspde.precision
#' @title Precision matrices for \code{inla.rspde} objects
#' @description Precision matrices for rSPDE models
#'
#' Calculates the precision matrix
#' for given parameter values based on an \code{inla.rspde} model object.
#' @param rspde An \code{inla.rspde} object.
#' @param theta The parameter vector. See the details in
#' \code{\link{rspde.matern}} to see the parameterizations.
#' @param optimized Logical indicating if only the elements
#' (the \code{x} slot) of the precision
#' matrix should be returned.
#' @return A sparse precision matrix.
#' @export
#' @examples
#' \donttest{
#' # devel version
#' library(INLA)
#'
#' set.seed(1)
#' n <- 10
#'
#' coords <- cbind(long = sample(1:n), lat = sample(1:n))
#'
#' mesh <- inla.mesh.2d(coords, max.edge = c(20, 40))
#' rspde_model_int <- rspde.matern(mesh = mesh, nu = 1)
#'
#' prec_int <- rspde.precision(rspde_model_int, theta = log(c(1, 3)))
#'
#' rspde_model <- rspde.matern(mesh)
#' prec <- rspde.precision(rspde_model, theta = log(c(1, 3, 1.2)))
#' # devel.tag
#' }
rspde.precision <- function(rspde,
                            theta,
                            optimized = FALSE) {
  check_class_inla_rspde(rspde)
  stopifnot(is.logical(optimized))
  if (length(theta) != length(rspde$f$rgeneric$definition(cmd = "initial"))) {
    stop("Length of theta is incorrect!")
  }
  if (optimized) {
    if (rspde$rspde$f$rgeneric$optimize) {
      return(rspde$f$rgeneric$definition(cmd = "Q", theta = theta))
    } else {
      return(rspde$f$rgeneric$definition(cmd = "Q", theta = theta)@x)
    }
  } else {
    if (rspde$f$rgeneric$optimize) {
      graph <- rspde$f$rgeneric$definition(cmd = "graph")
      entries <- rspde$f$rgeneric$definition(cmd = "Q", theta = theta)
      return(build_sparse_matrix_rspde(entries = entries, graph = graph))
    } else {
      return(rspde$f$rgeneric$definition(cmd = "Q", theta = theta))
    }
  }
}

#' @name rspde.result
#' @title rSPDE result extraction from INLA estimation results
#' @description Extract field and parameter values and distributions
#' for an rspde effect from an inla result object.
#' @param inla An \code{inla} object obtained from a call to
#' \code{inla()}.
#' @param name A character string with the name of the rSPDE effect
#' in the inla formula.
#' @param rspde The \code{inla.rspde} object used for the effect in
#' the inla formula.
#' @param compute.summary Should the summary be computed?
#' @return Returns a list containing:
#' \item{marginals.kappa}{Marginal densities for kappa}
#' \item{marginals.log.kappa}{Marginal densities for log(kappa)}
#' \item{marginals.log.tau}{Marginal densities for log(tau)}
#' \item{marginals.tau}{Marginal densities for tau}
#' \item{marginals.values}{Marginal densities for the field values}
#' \item{summary.log.kappa}{Summary statistics for log(kappa)}
#' \item{summary.log.tau}{Summary statistics for log(tau)}
#' \item{summary.values}{Summary statistics for the field values}
#' If nu was estimated, then the list will also contain
#' \item{marginals.nu}{Marginal densities for nu}
#' If nu was estimated and a beta prior was used, then the list will
#' also contain
#' \item{marginals.logit.nu}{Marginal densities for logit(nu)}
#' \item{summary.logit.nu}{Marginal densities for logit(nu)}
#' If nu was estimated and a truncated lognormal prior was used,
#' then the list will also contain
#' \item{marginals.log.nu}{Marginal densities for log(nu)}
#' \item{summary.log.nu}{Marginal densities for log(nu)}
#' If \code{compute.summary} is \code{TRUE}, then the list will also contain
#' \item{summary.kappa}{Summary statistics for kappa}
#' \item{summary.tau}{Summary statistics for tau}
#' If nu was estimated and \code{compute.summary} is \code{TRUE},
#' then the list will also contain
#' \item{summary.nu}{Summary statistics for nu}
#' @export
#' @examples
#' \donttest{
#' # devel version
#' library(INLA)
#' set.seed(123)
#'
#' m <- 100
#' loc_2d_mesh <- matrix(runif(m * 2), m, 2)
#' mesh_2d <- inla.mesh.2d(
#'   loc = loc_2d_mesh,
#'   cutoff = 0.05,
#'   max.edge = c(0.1, 0.5)
#' )
#' sigma <- 0.01
#' range <- 0.2
#' nu <- 0.8
#' kappa <- sqrt(8 * nu) / range
#' op <- matern.operators(
#'   mesh = mesh_2d, nu = nu,
#'   kappa = kappa, sigma = sigma, m = 2
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
#'   nu_upper_bound = 1
#' )
#' f <- y ~ -1 + f(field, model = rspde_model)
#' rspde_fit <- inla(f,
#'   data = inla.stack.data(st.dat),
#'   family = "gaussian",
#'   control.predictor =
#'     list(A = inla.stack.A(st.dat)),
#'            inla.mode = "experimental"
#' )
#' result <- rspde.result(rspde_fit, "field", rspde_model)
#' summary(result)
#' plot(result)
#' # devel.tag
#' }
rspde.result <- function(inla, name, rspde, compute.summary = TRUE) {
  check_class_inla_rspde(rspde)

  nu_upper_bound <- rspde$nu_upper_bound
  result <- list()

  if (!rspde$est_nu) {
    row_names <- c("tau", "kappa")
  } else {
    row_names <- c("tau", "kappa", "nu")
  }


  result$summary.values <- inla$summary.random[[name]]

  if (!is.null(inla$marginals.random[[name]])) {
    result$marginals.values <- inla$marginals.random[[name]]
  }


  result$summary.log.tau <- INLA::inla.extract.el(
    inla$summary.hyperpar,
    paste("Theta1 for ", name, "$", sep = "")
  )
  rownames(result$summary.log.tau) <- "log(tau)"
  result$summary.log.kappa <- INLA::inla.extract.el(
    inla$summary.hyperpar,
    paste("Theta2 for ", name, "$", sep = "")
  )
  rownames(result$summary.log.kappa) <- "log(kappa)"
  if (rspde$est_nu) {
    result$summary.logit.nu <- INLA::inla.extract.el(
      inla$summary.hyperpar,
      paste("Theta3 for ", name, "$", sep = "")
    )
    rownames(result$summary.logit.nu) <- "logit(nu)"
  }

  if (!is.null(inla$marginals.hyperpar[[paste0("Theta1 for ", name)]])) {
    result$marginals.log.tau <- INLA::inla.extract.el(
      inla$marginals.hyperpar,
      paste("Theta1 for ", name, "$", sep = "")
    )
    names(result$marginals.log.tau) <- "tau"
    result$marginals.log.kappa <- INLA::inla.extract.el(
      inla$marginals.hyperpar,
      paste("Theta2 for ", name, "$", sep = "")
    )
    names(result$marginals.log.kappa) <- "kappa"

    if (rspde$est_nu) {
      result$marginals.logit.nu <- INLA::inla.extract.el(
        inla$marginals.hyperpar,
        paste("Theta3 for ", name, "$", sep = "")
      )
      names(result$marginals.logit.nu) <- "nu"
    }

    result$marginals.tau <- lapply(
      result$marginals.log.tau,
      function(x) {
        INLA::inla.tmarginal(
          function(y) exp(y),
          x
        )
      }
    )
    result$marginals.kappa <- lapply(
      result$marginals.log.kappa,
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
              nu_upper_bound * exp(y) / (1 + exp(y))
            },
            x
          )
        }
      )
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
        subdivisions = nrow(density_df)
      )$value
      return(norm_const)
    }

    norm_const_tau <- norm_const(result$marginals.tau$tau)
    result$marginals.tau$tau[, "y"] <-
    result$marginals.tau$tau[, "y"] / norm_const_tau

    norm_const_kappa <- norm_const(result$marginals.kappa$kappa)
    result$marginals.kappa$kappa[, "y"] <-
    result$marginals.kappa$kappa[, "y"] / norm_const_kappa




    result$summary.tau <- create_summary_from_density(result$marginals.tau$tau,
    name = "tau")
    result$summary.kappa <-
    create_summary_from_density(result$marginals.kappa$kappa, name = "kappa")
    if (rspde$est_nu) {
      norm_const_nu <- norm_const(result$marginals.nu$nu)
      result$marginals.nu$nu[, "y"] <-
      result$marginals.nu$nu[, "y"] / norm_const_nu

      result$summary.nu <- create_summary_from_density(result$marginals.nu$nu,
      name = "nu")
    }
  }

  class(result) <- "rspde.result"
  return(result)
}

#' @name plot.rspde.result
#' @title Posterior plots for field parameters for an \code{inla.rspde} model
#' from a \code{rspde.result} object
#' @description Posterior plots for rSPDE field parameters in their
#' original scales.
#' @param x A \code{rspde.result} object.
#' @param which For which parameters the posterior should be plotted?
#' @param caption captions to appear above the plots; character
#' vector or list of
#' valid graphics annotations. Can be set to "" or NA to suppress all captions.
#' @param sub.caption	common title-above the figures if there are more than one.
#' @param type_plot what type of plot should be drawn. The default is 'l'.
#' @param ask logical; if \code{TRUE}, the user is asked before each plot.
#' @param main character; title to be placed at each plot additionally
#' (and above) all captions.
#' @param cex.caption	controls the size of caption.
#' @param cex.oma.main controls the size of the sub.caption only if
#' that is above the figures when there is more than one.
#' @param ylab Label for y axis.
#' @param xlab Label for x axis.
#' @param ... Additional arguments.
#' @return Called for its side effects.
#' @export
#' @method plot rspde.result
#' @examples
#' \donttest{
#' # devel version
#' library(INLA)
#' set.seed(123)
#'
#' m <- 100
#' loc_2d_mesh <- matrix(runif(m * 2), m, 2)
#' mesh_2d <- inla.mesh.2d(
#'   loc = loc_2d_mesh,
#'   cutoff = 0.05,
#'   max.edge = c(0.1, 0.5)
#' )
#' sigma <- 0.01
#' range <- 0.2
#' nu <- 0.8
#' kappa <- sqrt(8 * nu) / range
#' op <- matern.operators(
#'   mesh = mesh_2d, nu = nu,
#'   kappa = kappa, sigma = sigma, m = 2
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
#'   nu_upper_bound = 1
#' )
#' f <- y ~ -1 + f(field, model = rspde_model)
#' rspde_fit <- inla(f,
#'   data = inla.stack.data(st.dat),
#'   family = "gaussian",
#'   control.predictor =
#'     list(A = inla.stack.A(st.dat)),
#'            inla.mode = "experimental"
#' )
#' result <- rspde.result(rspde_fit, "field", rspde_model)
#' plot(result)
#' # devel.tag
#' }
plot.rspde.result <- function(x, which = c("tau", "kappa", "nu"),
                              caption = list(
                                "Posterior density for tau",
                                "Posterior density for kappa",
                                "Posterior density for nu"
                              ),
                              sub.caption = NULL,
                              type_plot = "l",
                              ask = prod(graphics::par("mfcol")) <
                                length(which) && grDevices::dev.interactive(),
                              main = "",
                              cex.oma.main = 1.25,
                              cex.caption = 1,
                              ylab = "Density",
                              xlab = "x",
                              ...) {
  result <- x
  which <- which[which %in% c("tau", "kappa", "nu")]
  stopifnot(!is.null(which))

  one.fig <- prod(graphics::par("mfcol")) == 1

  if (ask) {
    oask <- grDevices::devAskNewPage(TRUE)
    on.exit(grDevices::devAskNewPage(oask))
  }

  getCaption <- function(k) {
    if (length(caption) < k) {
      NA_character_
    } else {
      grDevices::as.graphicsAnnot(caption[[k]])
    }
  }

  if ("nu" %in% which) {
    if (is.null(result$marginals.nu)) {
      which <- which[which != "nu"]
    }
  }

  param <- c("tau", "kappa", "nu")
  for (i in 1:3) {
    if (param[i] %in% which) {
      graph_temp <- result[[paste0("marginals.", param[i])]][[param[i]]]
      graphics::plot(graph_temp,
        type = type_plot, main = main, ylab = ylab, xlab = xlab,
        ...
      )
      graphics::mtext(getCaption(i), side = 3, cex = cex.caption)
    }


    if (one.fig) {
      graphics::title(sub = sub.caption, ...)
    }
  }

  if (!one.fig && graphics::par("oma")[3L] >= 1) {
    graphics::mtext(sub.caption, outer = TRUE, cex = 1.25)
  }

  grDevices::dev.flush()

  invisible()
}


#' @name summary.rspde.result
#' @title Summary for posteriors of field parameters for an \code{inla.rspde}
#' model from a \code{rspde.result} object
#' @description Summary for posteriors of rSPDE field parameters in
#' their original scales.
#' @param object A \code{rspde.result} object.
#' @param digits integer, used for number formatting with signif()
#' @param ... Currently not used.
#' @return Returns a \code{data.frame}
#' containing the summary.
#' @export
#' @method summary rspde.result
#' @examples
#' \donttest{
#' # devel version
#' library(INLA)
#' set.seed(123)
#'
#' m <- 100
#' loc_2d_mesh <- matrix(runif(m * 2), m, 2)
#' mesh_2d <- inla.mesh.2d(
#'   loc = loc_2d_mesh,
#'   cutoff = 0.05,
#'   max.edge = c(0.1, 0.5)
#' )
#' sigma <- 0.01
#' range <- 0.2
#' nu <- 0.8
#' kappa <- sqrt(8 * nu) / range
#' op <- matern.operators(
#'   mesh = mesh_2d, nu = nu,
#'   kappa = kappa, sigma = sigma, m = 2
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
#'   nu_upper_bound = 1
#' )
#' f <- y ~ -1 + f(field, model = rspde_model)
#' rspde_fit <- inla(f,
#'   data = inla.stack.data(st.dat),
#'   family = "gaussian",
#'   control.predictor =
#'     list(A = inla.stack.A(st.dat)),
#'            inla.mode = "experimental"
#' )
#' result <- rspde.result(rspde_fit, "field", rspde_model)
#' summary(result)
#' # devel.tag
#' }
#'
summary.rspde.result <- function(object,
                                 digits = 6,
                                 ...) {
  if (is.null(object$summary.tau)) {
    warning("The summary was not computed, rerun rspde.result with
    compute.summary set to TRUE.")
  } else {
    out <- object$summary.tau
    out <- rbind(out, object$summary.kappa)
    if (!is.null(object$summary.nu)) {
      out <- rbind(out, object$summary.nu)
    }
    return(signif(out, digits = digits))
  }
}




#' @name rspde.mesh.project
#' @title Calculate a lattice projection to/from an \code{inla.mesh} for
#' rSPDE objects
#' @aliases rspde.mesh.project rspde.mesh.projector rspde.mesh.project.inla.mesh
#' rspde.mesh.project.rspde.mesh.projector rspde.mesh.project.inla.mesh.1d
#' @description Calculate a lattice projection to/from an \code{inla.mesh} for
#' rSPDE objects
#' @param mesh An \code{inla.mesh} or \code{inla.mesh.1d} object.
#' @param nu The smoothness parameter. If \code{NULL}, it will be assumed that
#' nu was estimated.
#' @param rspde_order The order of the rational approximation.
#' @param loc	Projection locations. Can be a matrix or a SpatialPoints or a
#' SpatialPointsDataFrame object.
#' @param field Basis function weights, one per mesh basis function, describing
#' the function to be evaluated at the projection locations.
#' @param projector A \code{rspde.mesh.projector} object.
#' @param lattice An \code{inla.mesh.lattice} object.
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
                                 rspde_order = 2,
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
  out <- do.call(INLA::inla.mesh.projector, args_list)
  if (mesh$manifold == "R1") {
    dim <- 1
  } else if (mesh$manifold == "R2") {
    dim <- 2
  } else {
    stop("The domain must be flat manifolds of dimension 1 or 2, that is,
         the domain must be a line or a plane.")
  }

  out$proj$A <- rspde.make.A(
    A = out$proj$A, rspde_order = rspde_order, dim = dim,
    nu = nu
  )

  class(out) <- c(class(out), "rspde.mesh.projector")
  return(out)
}


#' @rdname rspde.mesh.project
#' @method rspde.mesh.project inla.mesh
#' @export

rspde.mesh.project.inla.mesh <- function(mesh, loc = NULL,
                                         field = NULL, rspde_order = 2,
                                         nu = NULL, ...) {
  cond1 <- inherits(mesh, "inla.mesh.1d")
  cond2 <- inherits(mesh, "inla.mesh")
  stopifnot(cond1 || cond2)

  if (!missing(field) && !is.null(field)) {
    proj <- rspde.mesh.projector(mesh,
      loc = loc, rspde_order = rspde_order, nu = nu,
      ...
    )
    return(INLA::inla.mesh.project(proj, field = field))
  }
  jj <- which(rowSums(matrix(is.na(as.vector(loc)),
    nrow = nrow(loc),
    ncol = ncol(loc)
  )) == 0)
  smorg <- (INLA::inla.fmesher.smorg(mesh$loc,
  mesh$graph$tv, points2mesh = loc[jj, ,
    drop = FALSE
  ]))
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
      Abar <- kronecker(matrix(1, 1, rspde_order + 1), A)
    }
  } else {
    Abar <- kronecker(matrix(1, 1, rspde_order + 1), A)
  }

  list(t = ti, bary = b, A = Abar, ok = ok)
}


#' @rdname rspde.mesh.project
#' @method rspde.mesh.project rspde.mesh.projector
#' @export
#'

rspde.mesh.project.rspde.mesh.projector <- function(projector, field, ...) {
  return(INLA::inla.mesh.project(projector = projector, field = field, ...))
}



#' @rdname rspde.mesh.project
#' @method rspde.mesh.project inla.mesh.1d
#' @export
#'

rspde.mesh.project.inla.mesh.1d <- function(mesh, loc, field = NULL,
                                            rspde_order = 2, nu = NULL, ...) {
  stopifnot(inherits(mesh, "inla.mesh.1d"))
  if (!missing(field) && !is.null(field)) {
    proj <- rspde.mesh.projector(mesh, loc,
    rspde_order = rspde_order, nu = nu, ...)
    return(INLA::inla.mesh.project(proj, field))
  }
  A <- INLA::inla.mesh.1d.A(mesh, loc = loc)
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
      Abar <- kronecker(matrix(1, 1, rspde_order + 1), A)
    }
  } else {
    Abar <- kronecker(matrix(1, 1, rspde_order + 1), A)
  }
  return(list(A = Abar, ok = (loc >= mesh$interval[1]) & (loc <=
    mesh$interval[2])))
}



