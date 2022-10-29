
#' @name rspde.matern
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

  d <- get_inla_mesh_dimension(mesh)


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
  } else {
    integer_alpha <- FALSE

    rational_table <- get_rational_coefficients(rspde_order, type.rational.approx)
  }



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

    # if (integer_alpha) {
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
      graph_opt <- get.sparsity.graph.rspde(
        fem_mesh_matrices = fem_mesh_orig, dim = d,
        nu = nu_upper_bound,
        rspde_order = rspde_order,
        sharp = sharp,
        force_non_integer = TRUE
      )


    graph_opt <- transpose_cgeneric(graph_opt) 

    rspde_lib <- system.file('libs', package='rSPDE')

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
            prior.kappa.meanlog = prior.kappa$meanlog,
            prior.kappa.sdlog = prior.kappa$sdlog,
            prior.tau.meanlog = prior.tau$meanlog,
            prior.tau.sdlog = prior.tau$sdlog,
            prior.nu.loglocation = prior.nu$loglocation,
            prior.nu.mean = prior.nu$mean,
            prior.nu.prec = prior.nu$prec,
            prior.nu.logscale = prior.nu$logscale,
            start.lkappa = start.lkappa,
            start.ltau = start.ltau,
            start.nu = start.nu,
            rspde_order = as.integer(rspde_order),
            prior.nu.dist = prior.nu.dist))
    
    model$cgeneric_type <- "general"
  } else if (!integer_alpha) {
      graph_opt <- get.sparsity.graph.rspde(
        fem_mesh_matrices = fem_mesh_orig, dim = d,
        nu = nu,
        rspde_order = rspde_order,
        sharp = sharp,
        force_non_integer = TRUE
      )


    graph_opt <- transpose_cgeneric(graph_opt) 

    rspde_lib <- system.file('libs', package='rSPDE')

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
            prior.kappa.meanlog = prior.kappa$meanlog,
            prior.kappa.sdlog = prior.kappa$sdlog,
            prior.tau.meanlog = prior.tau$meanlog,
            prior.tau.sdlog = prior.tau$sdlog,
            start.lkappa = start.lkappa,
            start.ltau = start.ltau,
            rspde_order = as.integer(rspde_order),
            d = as.integer(d)))
    
    model$cgeneric_type <- "frac_alpha"
  } else {
      graph_opt <- get.sparsity.graph.rspde(
        fem_mesh_matrices = fem_mesh_orig, dim = d,
        nu = nu,
        rspde_order = rspde_order,
        sharp = sharp
      )
      graph_opt <- transpose_cgeneric(graph_opt) 

    rspde_lib <- system.file('libs', package='rSPDE')

    old_matrices_less <- matrices_less

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
            prior.kappa.meanlog = prior.kappa$meanlog,
            prior.kappa.sdlog = prior.kappa$sdlog,
            prior.tau.meanlog = prior.tau$meanlog,
            prior.tau.sdlog = prior.tau$sdlog,
            start.lkappa = start.lkappa,
            start.ltau = start.ltau #,
            # positions_C = positions_matrices_less[[1]],
            # positions_G = positions_matrices_less[[2]]
            ))
    model$cgeneric_type <- "int_alpha"
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
#' If `compute.summary` is `TRUE`, then the list will also contain
#' \item{summary.kappa}{Summary statistics for kappa}
#' \item{summary.tau}{Summary statistics for tau}
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


#' Create a posterior data frame for rspde.result objects
#'
#' FUNCTION_DESCRIPTION
#'
#' @param rspde_result An rspde.result object.
#' @param parameter Vector. Which parameters to get the posterior density in the data.frame? The options are `tau`, `kappa` and `nu`.
#' @param transform Should the posterior density be given in the original scale?
#'
#' @return A data frame containing the posterior densities.
#' @examples
#' # ADD_EXAMPLES_HERE
#' @export
posterior_df <- function(rspde_result, 
                          parameter = c("tau", "kappa", "nu"),
                          transform = TRUE) {
      if(!inherits(rspde_result, "rspde.result")){
        stop("The argument rspde_result should be of class rspde.result!")
      }
      parameter <- intersect(parameter, c("tau", "kappa", "nu"))
      if(length(parameter) == 0){
        stop("You should choose at least one of the parameters 'tau', 'kappa' or 'nu'!")
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



