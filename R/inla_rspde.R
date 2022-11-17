
#' @name rspde.matern
#' @title Matern rSPDE model object for INLA
#' @description Creates an INLA object for a stationary Matern model with
#' general smoothness parameter.
#' @param mesh The mesh to build the model. It can be an `inla.mesh` or
#' an `inla.mesh.1d` object. Otherwise, should be a list containing elements d, the dimension, C, the mass matrix,
#' and G, the stiffness matrix.
#' @param nu_upper_bound Upper bound for the smoothness parameter.
#' @param rspde_order The order of the covariance-based rational SPDE approach.
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
#'
#' @return An INLA model.
#' @export

rspde.matern <- function(mesh,
                         nu_upper_bound = 4, rspde_order = 2,
                         nu = NULL, 
                         debug = FALSE,
                         prior.kappa = NULL,
                         prior.nu = NULL,
                         prior.tau = NULL,
                         prior.range = NULL,
                         prior.std.dev = NULL,
                         start.lkappa = NULL,
                         start.nu = NULL,
                         start.ltau = NULL,
                         start.lrange = NULL,
                         start.lstd.dev = NULL,
                         parameterization = c("matern", "spde"),
                         prior.nu.dist = c("lognormal", "beta"),
                         nu.prec.inc = 1,
                         type.rational.approx = c("chebfun",
                         "brasil", "chebfunLB")) {
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

  integer.nu <- FALSE

  if(inherits(mesh, c("inla.mesh", "inla.mesh.1d"))){
    d <- get_inla_mesh_dimension(mesh)
  } else if (!is.null(mesh$d)){
    d <- mesh$d
  } else{
    stop("The mesh object should either be an INLA mesh object or contain d, the dimension!")
  }



  if (nu_upper_bound - floor(nu_upper_bound) == 0) {
    nu_upper_bound <- nu_upper_bound - 1e-5
  }
  fixed_nu <- !is.null(nu)
    if (fixed_nu) {
    nu_order <- nu
  } else {
    nu_order <- nu_upper_bound
  }

    beta <- nu_order / 2 + d / 4

    m_alpha <- floor(2 * beta)

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


  if (fixed_nu) {
    alpha <- nu + d / 2
    integer_alpha <- (alpha %% 1 == 0)
    if(!integer_alpha){
      if(rspde_order > 0){
      n_m <- rspde_order
            mt <- get_rational_coefficients(rspde_order, type.rational.approx)
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
    if(rspde_order > 0){
      rational_table <- get_rational_coefficients(rspde_order, type.rational.approx)
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
    } else if(rspde_order > 0) {
        result_sparsity <- analyze_sparsity_rspde(
          nu_upper_bound = nu_order, dim = d,
          rspde_order = rspde_order,
          fem_mesh_matrices = fem_mesh
        )
      positions_matrices <- result_sparsity$positions_matrices
    } else{
            result_sparsity <- analyze_sparsity_rspde(
          nu_upper_bound = nu_order, dim = d,
          rspde_order = rspde_order,
          fem_mesh_matrices = fem_mesh, 
          include_lower_order = FALSE
        )
      positions_matrices <- result_sparsity$positions_matrices
    }

    idx_matrices <- result_sparsity$idx_matrices

    if(rspde_order > 0 || integer_alpha){
      positions_matrices_less <- result_sparsity$positions_matrices_less
    }

    # if (integer_alpha) {
    if(rspde_order > 0 || integer_alpha){  
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

  # Prior nu

  if (is.null(prior.nu$loglocation)) {
    prior.nu$loglocation <- log(min(1, nu_upper_bound / 2))
  }

  if (is.null(prior.nu[["mean"]])) {
    prior.nu[["mean"]] <- min(1, nu_upper_bound / 2)
  }

  if (is.null(prior.nu$prec)) {
    mu_temp <- prior.nu[["mean"]] / nu_upper_bound
    prior.nu$prec <- max(1 / mu_temp, 1 / (1 - mu_temp)) + nu.prec.inc
  }

  if (is.null(prior.nu[["logscale"]])) {
    prior.nu[["logscale"]] <- 1
  }

  # Prior kappa and prior range

  if (is.null(prior.kappa$meanlog) && is.null(prior.range$meanlog)) {
    mesh.range <- ifelse(d == 2, (max(c(diff(range(mesh$loc[
      ,
      1
    ])), diff(range(mesh$loc[, 2])), diff(range(mesh$loc[
      ,
      3
    ]))))), diff(mesh$interval))
    prior.range.nominal <- mesh.range * 0.2

    prior.range$meanlog <- log(prior.range.nominal)

    if (prior.nu.dist == "lognormal") {
      prior.kappa$meanlog <- log(sqrt(8 *
      exp(prior.nu[["loglocation"]])) / prior.range.nominal)
    } else if (prior.nu.dist == "beta") {
      prior.kappa$meanlog <- log(sqrt(8 *
      prior.nu[["mean"]]) / prior.range.nominal)
    }
  } else if(is.null(prior.kappa$meanlog)){
    prior.range.nominal <- exp(prior.range$meanlog)
    if (prior.nu.dist == "lognormal") {
      prior.kappa$meanlog <- log(sqrt(8 *
      exp(prior.nu[["loglocation"]])) / prior.range.nominal)
    } else if (prior.nu.dist == "beta") {
      prior.kappa$meanlog <- log(sqrt(8 *
      prior.nu[["mean"]]) / prior.range.nominal)
    }
  } else{
      if (prior.nu.dist == "lognormal") {
      prior.range$meanlog <- log(sqrt(8 *
      exp(prior.nu[["loglocation"]]))) - prior.kappa$meanlog 
    } else if (prior.nu.dist == "beta") {
      prior.range$meanlog <- log(sqrt(8 *
      prior.nu[["mean"]])) - prior.kappa$meanlog
    }
  }

  if(is.null(prior.range$sdlog)){
    prior.range$sdlog <- sqrt(10)
  }

  if (is.null(prior.kappa$sdlog)) {
    prior.kappa$sdlog <- sqrt(10)
  }

  # Prior tau and prior std. dev

  if (is.null(prior.tau$meanlog) && is.null(prior.std.dev$meanlog)) {
    if (prior.nu.dist == "lognormal") {
      prior.tau$meanlog <- log(sqrt(gamma(exp(prior.nu[["loglocation"]])) /
      gamma(exp(prior.nu[["loglocation"]]) + d / 2) / (4 *
        pi * exp(prior.kappa$meanlog)^(2 * exp(prior.nu[["loglocation"]])))))
    } else if (prior.nu.dist == "beta") {
      prior.tau$meanlog <- log(sqrt(gamma(prior.nu[["mean"]]) /
      gamma(prior.nu[["mean"]] + d / 2) / (4 *
        pi * exp(prior.kappa$meanlog)^(2 * prior.nu[["mean"]]))))
    }

    prior.std.dev$meanlog <- 0
  } else if (is.null(prior.tau$meanlog)){
        if (prior.nu.dist == "lognormal") {
      prior.tau$meanlog <-  - prior.std.dev$meanlog + log(sqrt(gamma(exp(prior.nu[["loglocation"]])) /
      gamma(exp(prior.nu[["loglocation"]]) + d / 2) / (4 *
        pi * exp(prior.kappa$meanlog)^(2 * exp(prior.nu[["loglocation"]])))))
    } else if (prior.nu.dist == "beta") {
      prior.tau$meanlog <-  - prior.std.dev$meanlog + log(sqrt(gamma(prior.nu[["mean"]]) /
      gamma(prior.nu[["mean"]] + d / 2) / (4 *
        pi * exp(prior.kappa$meanlog)^(2 * prior.nu[["mean"]]))))
    }
  } else{
        prior.std.dev$meanlog <- 0
  }

  if (is.null(prior.tau$sdlog)) {
    prior.tau$sdlog <- sqrt(10)
  }


  if(is.null(prior.std.dev$sdlog)){
    prior.std.dev$sdlog <- sqrt(10)
  }

  # Starting values

  if (is.null(start.lkappa)) {
    start.lkappa <- prior.kappa$meanlog
  }
  if (is.null(start.ltau)) {
    start.ltau <- prior.tau$meanlog
  }

  if(is.null(start.lrange)){
    start.lrange <- prior.range$meanlog
  }
  if(is.null(start.lstd.dev)){
    start.lstd.dev <- prior.std.dev$meanlog
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
  
  rspde_lib <- system.file('shared', package='rSPDE')

  if(parameterization == "spde"){
    prior.theta1 <- prior.tau
    prior.theta2 <- prior.kappa
    start.theta1 <- start.ltau
    start.theta2 <- start.lkappa
  } else{
    prior.theta1 <- prior.std.dev
    prior.theta2 <- prior.range
    start.theta1 <- start.lstd.dev
    start.theta2 <- start.lrange
  }

  if (!fixed_nu) {

    if(rspde_order == 0){

      # fem_mesh has already been transposed
        graph_opt <-  fem_mesh[[paste0("g", m_alpha + 1)]]

  model <- do.call(
        'inla.cgeneric.define',
        list(model="inla_cgeneric_rspde_stat_parsim_gen_model",
            shlib=paste0(rspde_lib, '/rspde_cgeneric_models.so'),
            n=as.integer(n_cgeneric), debug=debug,
            d = as.double(d),
            nu_upper_bound = nu_upper_bound,
            matrices_full = as.double(matrices_full),
            graph_opt_i = graph_opt@i,
            graph_opt_j = graph_opt@j,
            prior.theta1.meanlog = prior.theta1$meanlog,
            prior.theta1.sdlog = prior.theta1$sdlog,
            prior.theta2.meanlog = prior.theta2$meanlog,
            prior.theta2.sdlog = prior.theta2$sdlog,
            prior.nu.loglocation = prior.nu$loglocation,
            prior.nu.mean = prior.nu$mean,
            prior.nu.prec = prior.nu$prec,
            prior.nu.logscale = prior.nu$logscale,
            start.theta1 = start.theta1,
            start.theta2 = start.theta2,
            start.nu = start.nu,
            prior.nu.dist = prior.nu.dist,
            parameterization = parameterization))


    } else{
          graph_opt <- get.sparsity.graph.rspde(
        fem_mesh_matrices = fem_mesh_orig, dim = d,
        nu = nu_upper_bound,
        rspde_order = rspde_order,
        force_non_integer = TRUE
      )


    graph_opt <- transpose_cgeneric(graph_opt) 


    # matrices_less <- restructure_matrices_less(matrices_less, m_alpha)
    # matrices_full <- restructure_matrices_full(matrices_full, m_alpha)

    model <- do.call(
        'inla.cgeneric.define',
        list(model="inla_cgeneric_rspde_stat_general_model",
            shlib=paste0(rspde_lib, '/rspde_cgeneric_models.so'),
            n=as.integer(n_cgeneric)*(rspde_order+1), debug=debug,
            d = as.double(d),
            nu_upper_bound = nu_upper_bound,
            matrices_less = as.double(matrices_less),
            matrices_full = as.double(matrices_full),
            rational_table = as.matrix(rational_table),
            graph_opt_i = graph_opt@i,
            graph_opt_j = graph_opt@j,
            prior.theta1.meanlog = prior.theta1$meanlog,
            prior.theta1.sdlog = prior.theta1$sdlog,
            prior.theta2.meanlog = prior.theta2$meanlog,
            prior.theta2.sdlog = prior.theta2$sdlog,
            prior.nu.loglocation = prior.nu$loglocation,
            prior.nu.mean = prior.nu$mean,
            prior.nu.prec = prior.nu$prec,
            prior.nu.logscale = prior.nu$logscale,
            start.theta1 = start.theta1,
            start.theta2 = start.theta2,
            start.nu = start.nu,
            rspde_order = as.integer(rspde_order),
            prior.nu.dist = prior.nu.dist,
            parameterization = parameterization))
    }
    
    model$cgeneric_type <- "general"
  } else if (!integer_alpha) {

      if(rspde_order == 0){
        graph_opt <- fem_mesh[[paste0("g", m_alpha + 1)]]

        model <- do.call(
        'inla.cgeneric.define',
        list(model="inla_cgeneric_rspde_stat_parsim_fixed_model",
            shlib=paste0(rspde_lib, '/rspde_cgeneric_models.so'),
            n=as.integer(n_cgeneric), debug=debug,
            d = as.double(d),
            nu = nu,
            matrices_full = as.double(matrices_full),
            graph_opt_i = graph_opt@i,
            graph_opt_j = graph_opt@j,
            prior.theta1.meanlog = prior.theta1$meanlog,
            prior.theta1.sdlog = prior.theta1$sdlog,
            prior.theta2.meanlog = prior.theta2$meanlog,
            prior.theta2.sdlog = prior.theta2$sdlog,
            start.theta1 = start.theta1,
            start.theta2 = start.theta2,
            parameterization = parameterization))

      } else{

      graph_opt <- get.sparsity.graph.rspde(
        fem_mesh_matrices = fem_mesh_orig, dim = d,
        nu = nu,
        rspde_order = rspde_order,
        force_non_integer = TRUE
      )


    graph_opt <- transpose_cgeneric(graph_opt) 

    # matrices_less <- restructure_matrices_less(matrices_less, m_alpha)
    # matrices_full <- restructure_matrices_full(matrices_full, m_alpha)

    model <- do.call(
        'inla.cgeneric.define',
        list(model="inla_cgeneric_rspde_stat_frac_model",
            shlib=paste0(rspde_lib, '/rspde_cgeneric_models.so'),
            n=as.integer(n_cgeneric)*(rspde_order+1), debug=debug,
            nu = nu,
            matrices_less = as.double(matrices_less),
            matrices_full = as.double(matrices_full),
            r_ratapprox = as.vector(r),
            p_ratapprox = as.vector(p),
            k_ratapprox = k,
            graph_opt_i = graph_opt@i,
            graph_opt_j = graph_opt@j,
            prior.theta1.meanlog = prior.theta1$meanlog,
            prior.theta1.sdlog = prior.theta1$sdlog,
            prior.theta2.meanlog = prior.theta2$meanlog,
            prior.theta2.sdlog = prior.theta2$sdlog,
            start.theta1 = start.theta1,
            start.theta2 = start.theta2,
            rspde_order = as.integer(rspde_order),
            parameterization = parameterization,
            d = as.integer(d)))
      }
    
    model$cgeneric_type <- "frac_alpha"
  } else {
      graph_opt <- get.sparsity.graph.rspde(
        fem_mesh_matrices = fem_mesh_orig, dim = d,
        nu = nu,
        rspde_order = rspde_order
      )
      graph_opt <- transpose_cgeneric(graph_opt) 

    # matrices_less <- restructure_matrices_less(matrices_less, m_alpha)

    model <- do.call(
        'inla.cgeneric.define',
        list(model="inla_cgeneric_rspde_stat_int_model",
            shlib=paste0(rspde_lib, '/rspde_cgeneric_models.so'),
            n=as.integer(n_cgeneric), debug=debug,
            matrices_less = as.double(matrices_less),
            m_alpha = as.integer(m_alpha),
            graph_opt_i = graph_opt@i,
            graph_opt_j = graph_opt@j,
            prior.theta1.meanlog = prior.theta1$meanlog,
            prior.theta1.sdlog = prior.theta1$sdlog,
            prior.theta2.meanlog = prior.theta2$meanlog,
            prior.theta2.sdlog = prior.theta2$sdlog,
            start.theta1 = start.theta1,
            start.theta2 = start.theta2,
            nu = nu,
            parameterization = parameterization
            ))
    model$cgeneric_type <- "int_alpha"
  }

  model$nu <- nu
  model$prior.theta1 <- prior.theta1
  model$prior.nu <- prior.nu
  model$prior.theta2 <- prior.theta2
  model$start.theta1 <- start.theta1
  model$start.theta2 <- start.theta2
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
  model$debug <- debug
  model$type.rational.approx <- type.rational.approx
  model$mesh <- mesh
  model$fem_mesh <- fem_mesh
  model$parameterization <- parameterization
  return(model)
}

#' @noRd 

transpose_cgeneric <- function(Cmatrix){
    Cmatrix <- INLA::inla.as.sparse(Cmatrix)
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

#' @name rspde.make.A
#' @title Observation/prediction matrices for rSPDE models.
#' @description Constructs observation/prediction weight matrices
#' for rSPDE models based on `inla.mesh` or
#' `inla.mesh.1d` objects.
#' @param mesh An `inla.mesh` or
#' an `inla.mesh.1d` object.
#' @param loc Locations, needed if an INLA mesh is provided
#' @param A The A matrix from the standard SPDE approach, such as the matrix
#' returned by `inla.spde.make.A`. Should only be provided if
#' `mesh` is not provided.
#' @param dim the dimension. Should only be provided if an
#' `mesh` is not provided.
#' @param rspde_order The order of the covariance-based rational SPDE approach.
#' @param nu If `NULL`, then the model will assume that nu will
#' be estimated. If nu is fixed, you should provide the value of nu.
#' @param index For each observation/prediction value, an index into loc.
#' Default is `seq_len(nrow(A.loc))`.
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
#' A <- rspde.make.A(mesh, loc = loc, rspde_order = 3)
#' }
#' #devel.tag
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
    dim <- get_inla_mesh_dimension(mesh)  
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
    if(rspde_order > 0){
      Abar <- kronecker(matrix(1, 1, rspde_order + 1), A)
    } else{
      Abar <- A
    }
      integer_nu <- FALSE
    }
  } else {
    if(rspde_order > 0){
      Abar <- kronecker(matrix(1, 1, rspde_order + 1), A)
    } else{
      Abar <- A
    }
    integer_nu <- FALSE
  }

  if (integer_nu) {
    rspde_order <- 0
  }


  attr(Abar, "inla_rspde_Amatrix") <- TRUE
  attr(Abar, "rspde_order") <- rspde_order
  attr(Abar, "integer_nu") <- integer_nu
  return(Abar)
}


#' @name rspde.make.index
#' @title rSPDE model index vector generation
#' @description Generates a list of named index vectors for an rSPDE model.
#' @param name A character string with the base name of the effect.
#' @param mesh An `inla.mesh` or
#' an `inla.mesh.1d` object.
#' @param rspde_order The order of the rational approximation
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
#'   nu_upper_bound = 2
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
#' }
#' #devel.tag
#' }
rspde.make.index <- function(name, n.spde = NULL, n.group = 1,
                             n.repl = 1, mesh = NULL,
                             rspde_order = 2, nu = NULL, dim = NULL) {
  if (is.null(n.spde) && is.null(mesh)) {
    stop("You should provide either n.spde or mesh!")
  }

  if (!is.null(mesh)) {
    n_mesh <- mesh$n
    dim <- get_inla_mesh_dimension(mesh)
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
      if(rspde_order > 0){
        factor_rspde <- rspde_order + 1
      } else{
        factor_rspde <- 1
      }
      integer_nu <- FALSE
    }
  } else {
      if(rspde_order > 0){
        factor_rspde <- rspde_order + 1
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
#'   nu_upper_bound = 2
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
#' }
#' #devel.tag
#' }
rspde.result <- function(inla, name, rspde, compute.summary = TRUE) {
  check_class_inla_rspde(rspde)

  nu_upper_bound <- rspde$nu_upper_bound
  result <- list()

  parameterization <- rspde$parameterization

  if (!rspde$est_nu) {
    if(parameterization == "spde"){
      row_names <- c("tau", "kappa")
    } else{
      row_names <- c("std.dev", "range")
    }
  } else {
    if(parameterization == "spde"){
      row_names <- c("tau", "kappa", "nu")
    } else{
      row_names <- c("std.dev", "range", "nu")
    }
  }


  result$summary.values <- inla$summary.random[[name]]

  if (!is.null(inla$marginals.random[[name]])) {
    result$marginals.values <- inla$marginals.random[[name]]
  }

  if(parameterization == "spde"){
    name_theta1 <- "tau"
    name_theta2 <- "kappa"
  } else{
    name_theta1 <- "std.dev"
    name_theta2 <- "range"
  }


  result[[paste0("summary.log.",name_theta1)]] <- INLA::inla.extract.el(
    inla$summary.hyperpar,
    paste("Theta1 for ", name, "$", sep = "")
  )
  rownames(  result[[paste0("summary.log.",name_theta1)]]) <- paste0("log(",name_theta1,")")
  
  result[[paste0("summary.log.",name_theta2)]] <- INLA::inla.extract.el(
    inla$summary.hyperpar,
    paste("Theta2 for ", name, "$", sep = "")
  )
  rownames(result[[paste0("summary.log.",name_theta2)]]) <- paste0("log(", name_theta2,")")
  if (rspde$est_nu) {
    result$summary.logit.nu <- INLA::inla.extract.el(
      inla$summary.hyperpar,
      paste("Theta3 for ", name, "$", sep = "")
    )
    rownames(result$summary.logit.nu) <- "logit(nu)"
  }

  if (!is.null(inla$marginals.hyperpar[[paste0("Theta1 for ", name)]])) {
    result[[paste0("marginals.log.",name_theta1)]] <- INLA::inla.extract.el(
      inla$marginals.hyperpar,
      paste("Theta1 for ", name, "$", sep = "")
    )
    names(result[[paste0("marginals.log.",name_theta1)]]) <- name_theta1
    result[[paste0("marginals.log.",name_theta2)]] <- INLA::inla.extract.el(
      inla$marginals.hyperpar,
      paste("Theta2 for ", name, "$", sep = "")
    )
    names(result[[paste0("marginals.log.",name_theta2)]]) <- name_theta2

    if (rspde$est_nu) {
      result$marginals.logit.nu <- INLA::inla.extract.el(
        inla$marginals.hyperpar,
        paste("Theta3 for ", name, "$", sep = "")
      )
      names(result$marginals.logit.nu) <- "nu"
    }

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

  class(result) <- "rspde.result"
  result$params <- c(name_theta1,name_theta2)
  if(rspde$est_nu){
    result$params <- c(result$params, "nu")
  }
  return(result)
}

#' @name plot.rspde.result
#' @title Posterior plots for field parameters for an `inla_rspde` model
#' from a `rspde.result` object
#' @description Posterior plots for rSPDE field parameters in their
#' original scales.
#' @param x A `rspde.result` object.
#' @param which For which parameters the posterior should be plotted?
#' @param caption captions to appear above the plots; character
#' vector or list of
#' valid graphics annotations. Can be set to "" or NA to suppress all captions.
#' @param sub.caption	common title-above the figures if there are more than one.
#' @param type_plot what type of plot should be drawn. The default is 'l'.
#' @param ask logical; if `TRUE`, the user is asked before each plot.
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
#'   nu_upper_bound = 2
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
#' }
#' #devel.tag
#' }
plot.rspde.result <- function(x, which = x$params,
                              caption = list(paste(
                                "Posterior density for", x$params[1]),
                                paste(
                                "Posterior density for", x$params[2]),
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
  which <- which[which %in% c("tau", "kappa", "nu", "range", "std.dev")]
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

  param <- result$params
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


#' Data frame for rspde.result objects to be used in ggplot2
#'
#' Returns a ggplot-friendly data-frame with the marginal posterior densities.
#'
#' @param rspde_result An rspde.result object.
#' @param parameter Vector. Which parameters to get the posterior density in the data.frame? The options are `std.dev`, `range`, `tau`, `kappa` and `nu`.
#' @param transform Should the posterior density be given in the original scale?
#' @param restrict_x_axis Variables to restrict the range of x axis based on quantiles.
#' @param restrict_quantiles List of quantiles to restrict x axis.
#'
#' @return A data frame containing the posterior densities.
#' @examples
#' # ADD_EXAMPLES_HERE
#' @export
gg_df <- function(rspde_result, 
                          parameter = rspde_result$params,
                          transform = TRUE,
                          restrict_x_axis = parameter,
                          restrict_quantiles = list(std.dev = c(0,1),
                          range = c(0,1),
                          nu = c(0,1),
                          kappa = c(0,1),
                          tau = c(0,1))) {
      if(!inherits(rspde_result, "rspde.result")){
        stop("The argument rspde_result should be of class rspde.result!")
      }
      parameter <- intersect(parameter, c("tau", "kappa", "nu", "range", "std.dev"))
      if(length(parameter) == 0){
        stop("You should choose at least one of the parameters 'tau', 'kappa', 'nu', 'range' or 'std.dev'!")
      }
    if ("nu" %in% parameter) {
    if (is.null(rspde_result$marginals.nu)) {
      parameter <- parameter[parameter != "nu"]
    }
  }
  
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
      warning("If you want to restrict x axis you should provide a quantile for the parameter!")
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
      warning(paste("No quantile for", parameter[[i]]))
      warning("If you want to restrict x axis you should provide a quantile for the parameter!")
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
  return(ret_df)
}


#' @name summary.rspde.result
#' @title Summary for posteriors of field parameters for an `inla_rspde`
#' model from a `rspde.result` object
#' @description Summary for posteriors of rSPDE field parameters in
#' their original scales.
#' @param object A `rspde.result` object.
#' @param digits integer, used for number formatting with signif()
#' @param ... Currently not used.
#' @return Returns a `data.frame`
#' containing the summary.
#' @export
#' @method summary rspde.result
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
#'   nu_upper_bound = 2
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
#' }
#' #devel.tag
#' }
#'
summary.rspde.result <- function(object,
                                 digits = 6,
                                 ...) {

  if (is.null(object[[paste0("summary.",object$params[1])]])) {
    warning("The summary was not computed, rerun rspde.result with
    compute.summary set to TRUE.")
  } else {
    out <- object[[paste0("summary.",object$params[1])]]
    out <- rbind(out, object[[paste0("summary.",object$params[2])]])
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
#' @param rspde_order The order of the rational approximation.
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
  dim <- get_inla_mesh_dimension(mesh)

  out$proj$A <- rspde.make.A(
    A = out$proj$A, rspde_order = rspde_order, dim = dim,
    nu = nu
  )

  class(out) <- c("rspde.mesh.projector",class(out))
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

  # r <- sapply(1:(n_m), function(i) {
  #   approx(mt$alpha, mt[[paste0("r", i)]], cut_decimals(2 * beta))$y
  # })

  # p <- sapply(1:(n_m), function(i) {
  #   approx(mt$alpha, mt[[paste0("p", i)]], cut_decimals(2 * beta))$y
  # })

  # k <- approx(mt$alpha, mt$k, cut_decimals(2 * beta))$y

  row_nu <- round(1000*cut_decimals(2*beta))
  r <- unlist(mt[row_nu, 2:(1+rspde_order)])
  p <- unlist(mt[row_nu, (2+rspde_order):(1+2*rspde_order)])
  k <- unlist(mt[row_nu, 2+2*rspde_order])


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
#' @param rspde_order The order of the rational approximation
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

#' @name rspde.precision
#' @title Precision matrices for `inla_rspde` objects
#' @description Precision matrices for rSPDE models
#'
#' Calculates the precision matrix
#' for given parameter values based on an `inla_rspde` model object.
#' @param rspde An `inla_rspde` object.
#' @param theta The parameter vector. See the details in
#' [rspde.matern()] to see the parameterizations.
#' @return A sparse precision matrix.
#' @export
#' @examples
#' \donttest{ #devel version
#' if (requireNamespace("INLA", quietly = TRUE)){
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
#' }
#' #devel.tag
#' }
rspde.precision <- function(rspde,
                            theta) {
  check_class_inla_rspde(rspde)
  if(!is.vector(theta)){
    stop("theta should be a vector!")
  }
  if(rspde$cgeneric_type == "general")  {
    if(length(theta)!=3){
      stop("The vector theta should have length 3!")
    }
    nu_upper_bound <- rspde$nu_upper_bound
    kappa = exp(theta[2])
    tau = exp(theta[1])
    nu <- (exp(theta[3]) / (1 + exp(theta[3]))) * nu_upper_bound
  return(rspde.matern.precision(kappa = kappa, nu = nu, tau=tau, rspde_order = rspde$rspde_order, fem_mesh_matrices = rspde$fem_mesh,
  dim = rspde$dim, type_rational_approx = rspde$type.rational.approx))
  } else{
    if(length(theta)!= 2){
      stop("The vector theta should have length 2!")
    }
    nu <- rspde$nu
    kappa = exp(theta[2])
    tau = exp(theta[1])
    if(rspde$cgeneric_type == "frac_alpha"){
    return(rspde.matern.precision(kappa = kappa, nu = nu, tau=tau, rspde_order = rspde$rspde_order, fem_mesh_matrices = rspde$fem_mesh,
  dim = rspde$dim, type_rational_approx = rspde$type.rational.approx))
    } else{
      return(rspde.matern.precision.integer(kappa = kappa, nu = nu, tau = tau, dim = rspde$dim, fem_mesh_matrices = rspde$fem_mesh))
    }
  }
  
}



#' @name rspde.metric_graph
#' @title Matern rSPDE model object for metric graphs in INLA
#' @description Creates an INLA object for a stationary Matern model on a metric graph with
#' general smoothness parameter.
#' @param graph_obj The graph object to build the model. Needs to be of class `GPGraph::graph`. It should have a built mesh.
#' If the mesh is not built, one will be built using h=0.01 as default.
#' @param h The width of the mesh in case the mesh was not built.
#' @param nu_upper_bound Upper bound for the smoothness parameter.
#' @param rspde_order The order of the covariance-based rational SPDE approach.
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
#'
#' @return An INLA model.
#' @export

rspde.metric_graph <- function(graph_obj,
                         h = NULL,
                         nu_upper_bound = 2, rspde_order = 2,
                         nu = NULL, 
                         debug = FALSE,
                         prior.kappa = NULL,
                         prior.nu = NULL,
                         prior.tau = NULL,
                         prior.range = NULL,
                         prior.std.dev = NULL,
                         start.lkappa = NULL,
                         start.nu = NULL,
                         start.ltau = NULL,
                         start.lrange = NULL,
                         start.lstd.dev = NULL,
                         parameterization = c("matern", "spde"),
                         prior.nu.dist = c("lognormal", "beta"),
                         nu.prec.inc = 1,
                         type.rational.approx = c("chebfun",
                         "brasil", "chebfunLB")) {
    if(!inherits(graph_obj, "GPGraph::graph")){
      stop("The graph object should be of class GPGraph::graph!")
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

    if(is.null(prior.range$meanlog)){
      if(is.null(graph_obj$geo.dist)){
        graph_obj$compute_geodist()
      }
      prior.range.nominal <- max(graph_obj$geo.dist) * 0.2
      prior.range$meanlog <- log(prior.range.nominal)
    }

    return(rspde.matern(mesh = list(d = 1, C = graph_obj$mesh$C, 
                                G = graph_obj$mesh$G),
                                nu_upper_bound = nu_upper_bound,
                                rspde_order = rspde_order,
                                nu = nu,
                                debug = debug,
                                prior.kappa = prior.kappa,
                                prior.nu = prior.nu,
                                prior.tau = prior.tau,
                                prior.range = prior.range,
                                prior.std.dev =  prior.std.dev,
                                start.lkappa = start.lkappa,
                                start.nu = start.nu,
                                start.ltau = start.ltau,
                                start.lrange = start.lrange,
                                start.lstd.dev = start.lstd.dev,
                                parameterization = parameterization,
                                prior.nu.dist = prior.nu.dist,
                                nu.prec.inc = nu.prec.inc,
                                type.rational.approx = type.rational.approx
                                ))

                         }