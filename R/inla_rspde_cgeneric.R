
#' @name rspde.matern_cgeneric
#' @title Matern rSPDE model object for INLA
#' @description Creates an INLA object for a stationary Matern model with
#' general smoothness parameter.
#' @param mesh The mesh to build the model. Should be an `inla.mesh` or
#' an `inla.mesh.1d` object.
#' @param nu_upper_bound Upper bound for the smoothness parameter.
#' @param rspde_order The order of the covariance-based rational SPDE approach.
#' @param nu If nu is set to a parameter, nu will be kept fixed and will not
#' be estimated. If nu is `NULL`, it will be estimated.
#' @param sharp The sparsity graph should have the correct sparsity (costs
#' more to perform a sparsity analysis) or an upper bound for the sparsity? If
#' `TRUE`, the graph will have the correct sparsity.
#' @param debug INLA debug argument.
#' @param optimize Should the model be optimized? In this case the sparsities of
#' the matrices will be analyzed.
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

rspde.matern_cgeneric <- function(mesh,
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

  d <- get_inla_mesh_dimension(mesh)

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

    n_cgeneric <- ncol(fem_mesh[["c0"]])

    fem_mesh_orig <- fem_mesh

    fem_mesh <- fem_mesh[setdiff(names(fem_mesh),c("ta","va"))]

    fem_mesh <- lapply(fem_mesh, transpose_cgeneric)

    if (integer_alpha) {
      result_sparsity <- analyze_sparsity_rspde(
        nu_upper_bound = nu_order, dim = d,
        rspde_order = rspde_order,
        fem_mesh_matrices = fem_mesh,
        include_higher_order = FALSE
      )
    } else {
        result_sparsity <- analyze_sparsity_rspde(
          nu_upper_bound = nu_order, dim = d,
          rspde_order = rspde_order,
          fem_mesh_matrices = fem_mesh
        )
      positions_matrices <- result_sparsity$positions_matrices
    }

    idx_matrices <- result_sparsity$idx_matrices
    positions_matrices_less <- result_sparsity$positions_matrices_less

    if (integer_alpha) {
    n_tmp <- length(
        fem_mesh[[paste0("g", m_alpha)]]@x[idx_matrices[[m_alpha + 1]]]
        )
    # tmp <- fem_mesh[[paste0("g", m_alpha)]]
    # ii <- tmp@i
    # tmp@i <- tmp@j
    # tmp@j <- ii
    # idx <- which(tmp@i <= tmp@j)
    # rm(ii)

    tmp <-
    rep(0, n_tmp)
    tmp[positions_matrices_less[[1]]] <-
    fem_mesh$c0@x[idx_matrices[[1]]]

    # matrices_less <- tmp[idx]

    matrices_less <- tmp

    tmp <-
    rep(0, n_tmp)
    tmp[positions_matrices_less[[2]]] <-
    fem_mesh$g1@x[idx_matrices[[2]]]

    # matrices_less <- c(matrices_less, tmp[idx])

    matrices_less <- c(matrices_less, tmp)

    if (m_alpha > 2) {
    for (j in 2:(m_alpha - 1)) {
        tmp <-
        rep(0, n_tmp)
        tmp[positions_matrices_less[[j + 1]]] <-
        fem_mesh[[paste0("g", j)]]@x[idx_matrices[[j + 1]]]
        # matrices_less <- c(matrices_less, tmp[idx])
        matrices_less <- c(matrices_less, tmp)
    }
    }

    tmp <- fem_mesh[[paste0("g", m_alpha)]]@x[idx_matrices[[m_alpha + 1]]]

    # matrices_less <- c(matrices_less, tmp[idx])
    matrices_less <- c(matrices_less, tmp)
    }


    if (!integer_alpha) {
      if (m_alpha == 0) {
        # tmp <- fem_mesh$g1
        # ii <- tmp@i
        # tmp@i <- tmp@j
        # tmp@j <- ii
        # idx <- which(tmp@i <= tmp@j)
        # rm(ii)
        n_tmp <- length(fem_mesh$g1@x)
        tmp <-
        rep(0, n_tmp)
        tmp[positions_matrices[[1]]] <-
        fem_mesh$c0@x[idx_matrices[[1]]]
        # matrices_full <- tmp[idx]
        # matrices_full <- c(matrices_full, fem_mesh$g1@x[idx])
        matrices_full <- tmp
        matrices_full <- c(matrices_full, fem_mesh$g1@x)
      } else if (m_alpha > 0) {
        # tmp <- fem_mesh[[paste0("g", m_alpha + 1)]]
        # ii <- tmp@i
        # tmp@i <- tmp@j
        # tmp@j <- ii
        # idx <- which(tmp@i <= tmp@j)
        # rm(ii)

        n_tmp <- length(
            fem_mesh[[paste0("g", m_alpha + 1)]]@x[idx_matrices[[m_alpha + 2]]]
        )

        tmp <- 
        rep(0, n_tmp)
        tmp[positions_matrices[[1]]] <-
        fem_mesh$c0@x[idx_matrices[[1]]]

        # matrices_full <- tmp[idx]
        
        matrices_full <- tmp
        
        tmp <-
        rep(0, n_tmp)
        tmp[positions_matrices[[2]]] <-
        fem_mesh$g1@x[idx_matrices[[2]]]

        # matrices_full <- c(matrices_full, tmp[idx])
        matrices_full <- c(matrices_full, tmp)

      if (m_alpha > 1) {
        for (j in 2:(m_alpha)) {
          tmp <-
          rep(0, length(fem_matrices[[paste0("G_", m_alpha + 1)]]))
          tmp[positions_matrices[[j + 1]]] <-
          fem_mesh[[paste0("g", j)]]@x[idx_matrices[[j + 1]]]
          # matrices_full <- c(matrices_full, tmp[idx])
          matrices_full <- c(matrices_full, tmp)
        }
      }

        tmp <-
        fem_mesh[[paste0("g", m_alpha + 1)]]@x[idx_matrices[[m_alpha + 2]]]
        # matrices_full <- c(matrices_full, tmp[idx])
        matrices_full <- c(matrices_full, tmp)
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
      graph_opt <- rSPDE::get.sparsity.graph.rspde(
        fem_mesh_matrices = fem_mesh, dim = d,
        nu = nu_upper_bound,
        rspde_order = rspde_order,
        sharp = sharp,
        force_non_integer = TRUE
      )

    # graph_opt <- transpose_cgeneric(graph_opt) # Already transposed everything above
    
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
      n = n_cgeneric * (rspde_order + 1),
      debug = debug,
      do_optimize = optimize, optimize = optimize
    )
    model$rgeneric_type <- "general"
    model$optimize <- optimize
  } else if (!integer_alpha) {
      graph_opt <- rSPDE::get.sparsity.graph.rspde(
        fem_mesh_matrices = fem_mesh, dim = d,
        nu = nu,
        rspde_order = rspde_order,
        sharp = sharp
      )
    graph_opt <- transpose_cgeneric(graph_opt)

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
      n = n_cgeneric * (rspde_order + 1),
      debug = debug,
      do_optimize = optimize, optimize = optimize
    )
    model$rgeneric_type <- "frac_alpha"
    model$optimize <- optimize
  } else {
      graph_opt <- rSPDE::get.sparsity.graph.rspde(
        fem_mesh_matrices = fem_mesh_orig, dim = d,
        nu = nu,
        rspde_order = rspde_order,
        sharp = sharp
      )
      graph_opt <- transpose_cgeneric(graph_opt) 

    rspde_lib <- system.file('libs', package='rSPDE')

    old_matrices_less <- matrices_less

    matrices_less <- restructure_matrices_less(matrices_less, m_alpha)

    model <- do.call(
        'inla.cgeneric.define',
        list(model="inla_cgeneric_rspde_stat_int_model",
            shlib=paste0(rspde_lib, '/rspde_cgeneric_models.so'),
            n=as.integer(n_cgeneric), debug=debug,
            matrices_less = as.double(matrices_less),
            m_alpha = as.integer(m_alpha),
            graph_opt_i = graph_opt@i,
            graph_opt_j = graph_opt@j,
            prior.kappa.meanlog = prior.kappa$meanlog,
            prior.kappa.sdlog = prior.kappa$sdlog,
            prior.tau.meanlog = prior.tau$meanlog,
            prior.tau.sdlog = prior.tau$sdlog,
            start.lkappa = start.lkappa,
            start.ltau = start.ltau))
    model$rgeneric_type <- "int_alpha"
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
  class(model) <- c("inla_rspde", class(model))
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
  model$matrices_less <- matrices_less
  model$old_matrices_less <- old_matrices_less

  return(model)
}

#' @noRd 

transpose_cgeneric <- function(Cmatrix){
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