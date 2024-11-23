
#' Rational approximations of stationary anisotropic Gaussian Matern random fields
#'
#' `rspde.anistropic2d` computes a Finite Element Method (FEM) approximation of a
#' Gaussian random field defined as the solution to the stochastic partial
#' differential equation (SPDE):
#' \deqn{C(h) = \frac{\sigma^2}{2^{\nu-1}\Gamma(\nu)}(\sqrt{h^T H^{-1}h})^\nu K_\nu(\sqrt{h^T H^{-1}h})},
#' based on a SPDE representation of the form
#' \deqn{(I - \nabla\cdot(H\nabla))^{(\nu+1)/2}u = c\sigma W},
#' where $c>0$ is a constant. The matrix \eqn{H} is defined as
#' \deqn{\begin{bmatrix}
#' h_x^2 & h_xh_yh_{xy} \\
#' h_xh_yh_{xy} & h_y^2
#' \end{bmatrix}}
#'
#' @param mesh Spatial mesh for the FEM approximation.
#' @param nu If nu is set to a parameter, nu will be kept fixed and will not
#' be estimated. If nu is `NULL`, it will be estimated.
#' @param nu.upper.bound Upper bound for the smoothness parameter \eqn{\nu}. If `NULL`, it will be set to 2.
#' @param rspde.order The order of the covariance-based rational SPDE approach. The default order is 1.
#' @param prior.hx A list specifying the prior for the parameter \eqn{h_x} in the matrix \eqn{H}. This list may contain two elements: `mean` and/or `precision`, both of which must be numeric scalars. The precision refers to the prior on \eqn{\log(h_x)}. If `NULL`, default values will be used. The `mean` value is also used as starting value for hx.
#' @param prior.hy A list specifying the prior for the parameter \eqn{h_y} in the matrix \eqn{H}. This list may contain two elements: `mean` and/or `precision`, both of which must be numeric scalars. The precision refers to the prior on \eqn{\log(h_x)}. If `NULL`, default values will be used. The `mean` value is also used as starting value for hy.
#' @param prior.hxy A list specifying the prior for the parameter \eqn{h_x} in the matrix \eqn{H}. This list may contain two elements: `mean` and/or `precision`, both of which must be numeric scalars. The precision refers to the prior on \eqn{\log((h_{xy}+1)/(1-h_{xy}))}. If `NULL`, default values will be used. The `mean` value is also used as starting value for hxy.
#' @param prior.sigma A list specifying the prior for the variance parameter \eqn{\sigma}.
#' This list may contain two elements: `mean` and/or `precision`, both of which must
#' be numeric scalars. The precision refers to the prior on \eqn{\log(\sigma)}. If `NULL`,
#' default values will be used. The `mean` value is also used as starting value for sigma.
#' @param prior.precision A precision matrix for \eqn{\log(h_x), \log(h_y), \log((h_{xy}+1)/(1-h_{xy})), \log(\sigma)}. This matrix replaces the precision
#' element from `prior.kappa`, `prior.sigma`, `prior.gamma`, and `prior.rho` respectively. For dimension 1 `prior.precision` must be a 4x4 matrix. For dimension 2, \eqn{\rho} is a vector of length 2, so in this case `prior.precision` must be a 5x5 matrix. If `NULL`, a diagonal precision matrix with default values will be used.
#' @param prior.nu a list containing the elements `mean` and `prec`
#' for beta distribution, or `loglocation` and `logscale` for a
#' truncated lognormal distribution. `loglocation` stands for
#' the location parameter of the truncated lognormal distribution in the log
#' scale. `prec` stands for the precision of a beta distribution.
#' `logscale` stands for the scale of the truncated lognormal
#' distribution on the log scale. Check details below.
#' @param prior.nu.dist The distribution of the smoothness parameter.
#' The current options are "beta" or "lognormal". The default is "lognormal".
#' @param nu.prec.inc Amount to increase the precision in the beta prior
#' distribution. Check details below.
#' @param type.rational.approx Which type of rational approximation
#' should be used? The current types are "chebfun", "brasil" or "chebfunLB".
#' @param shared_lib String specifying which shared library to use for the Cgeneric
#' implementation. Options are "detect", "INLA", or "rSPDE". You may also specify the
#' direct path to a .so (or .dll) file.
#' @param debug Logical value indicating whether to enable INLA debug mode.
#' @param ... Additional arguments passed internally for configuration purposes.
#' @return An object of class `inla_rspde_spacetime` representing the FEM approximation of
#' the space-time Gaussian random field.
#' @export
#'
#' @examples
#' library(fmesher)
#' n_loc <- 2000
#' loc_2d_mesh <- matrix(runif(n_loc * 2), n_loc, 2)
#' mesh_2d <- fm_mesh_2d(loc = loc_2d_mesh, cutoff = 0.03, max.edge = c(0.1, 0.5))
#' aniso_model <- rspde.anistropic2d(mesh = mesh_2d)
#'

rspde.anistropic2d <- function(mesh,
                               nu = NULL,
                               nu.upper.bound = 2,
                               rspde.order = 1,
                               prior.hx = NULL,
                               prior.hy = NULL,
                               prior.hxy = NULL,
                               prior.sigma = NULL,
                               prior.precision = NULL,
                               prior.nu = NULL,
                               prior.nu.dist = "lognormal",
                               nu.prec.inc = 0.01,
                               type.rational.approx = "chebfun",
                               shared_lib = "detect",
                               debug = FALSE,
                               ...) {
  # Validate mesh input
  if (is.null(mesh) || !inherits(mesh, "fm_mesh_2d")) {
    stop("'mesh' must be a valid spatial mesh of class 'fm_mesh_2d'.")
  }

  if (nu.upper.bound - floor(nu.upper.bound) == 0) {
    nu.upper.bound <- nu.upper.bound - 1e-5
  }
  
  if(!is.null(nu)){
    nu.upper.bound <- nu
  }

  op <- matern2d.operators(
    hx = prior.hx$mean,
    hy = prior.hy$mean,
    hxy = prior.hxy$mean,
    sigma = prior.sigma$mean,
    mesh = mesh,
    nu = nu.upper.bound,
    m = rspde.order,
    type_rational_approximation = type.rational.approx,
    return_fem_matrices = TRUE
  )

  default_precision <- 0.1

  prior.hx <- set_prior(prior.hx, op$hx, default_precision, p = 1)
  prior.hy <- set_prior(prior.hy, op$hy, default_precision, p = 1)
  prior.hxy <- set_prior(prior.hxy, op$hxy, default_precision, p = 1)
  prior.sigma <- set_prior(prior.sigma, op$sigma, default_precision, p = 1)

  est_nu <- 0L

  if(is.null(nu)){
    est_nu <- 1L
    nu <- -1.0
  }

  result_nu <- handle_prior_nu(prior.nu, nu.upper.bound = nu.upper.bound, nu.prec.inc = nu.prec.inc, prior.nu.dist = prior.nu.dist)

  prior.nu <- result_nu$prior.nu
  start.nu <- result_nu$start.nu

  # Set default precision matrix if prior.precision is NULL
  if (is.null(prior.precision)) {
      prior.precision <- diag(c(
        prior.hx$precision,
        prior.hy$precision,
        prior.hxy$precision,
        prior.sigma$precision
      ))
  } else {

    if (!is.matrix(prior.precision) || !all(dim(prior.precision) == c(4,4))) {
      stop("prior.precision must be a 4x4 matrix.")
    }
  }

  rational_table <- as.matrix(get_rational_coefficients(rspde.order, type.rational.approx))

  rspde_lib <- get_shared_library(shared_lib)

  Cmatrix  <-  as(as(as(op$Q, "dMatrix"), "generalMatrix"), "TsparseMatrix")
  ii <- Cmatrix@i
  Cmatrix@i <- Cmatrix@j
  Cmatrix@j <- ii
  idx <- which(Cmatrix@i <= Cmatrix@j)
  Cmatrix@i <- Cmatrix@i[idx]
  Cmatrix@j <- Cmatrix@j[idx]
  Cmatrix@x <- Cmatrix@x[idx]

  list_args <- 
    list(
      model = "inla_cgeneric_rspde_anisotropic2d_model",
      shlib = rspde_lib,
      n = nrow(op$Q),
      debug = debug,
      Q = Cmatrix,
      C = op$C,
      Ci = op$Ci,
      Hxx = op$Hxx,
      Hyy = op$Hyy,
      Hxy = op$Hxy,
      prior.hx.mean = prior.hx$mean,
      prior.hy.mean = prior.hy$mean,
      prior.hxy.mean = prior.hxy$mean,
      prior.sigma.mean = prior.sigma$mean,
      prior.precision = prior.precision,
      nu = nu,
      est_nu = est_nu,
      start_nu = start.nu,
      prior.nu.loglocation = prior.nu$loglocation,
      prior.nu.mean = prior.nu$mean,
      prior.nu.prec = prior.nu$prec,
      prior.nu.logscale = prior.nu$logscale,
      prior.nu.dist = prior.nu.dist,
      nu_upper_bound = nu.upper.bound,
      rspde_order = as.integer(rspde.order),
      rational_table = rational_table
    )

  model <- do.call(INLA::inla.cgeneric.define, list_args)

  model$prior.sigma <- prior.sigma
  model$prior.hx <- prior.hx
  model$prior.hy <- prior.hy
  model$prior.hxy <- prior.hxy
  model$prior.precision <- prior.precision
  model$mesh <- mesh
  model$rspde.order <- rspde.order
  model$type_rational_approximation <- type.rational.approx
  model$est_nu <- est_nu
  model$nu <- nu
  model$nu_upper_bound <- nu.upper.bound
  model$rspde_version <- as.character(packageVersion("rSPDE"))

  ### The following objects are provided for backward compatibility
  if(!est_nu){
    model$integer.nu <- (nu %% 1) == 0
  } else{
    model$integer.nu <- FALSE
  }
  model$n.spde <- mesh$n
  ### 

  class(model) <- c("inla_rspde_anisotropic2d", class(model))

  return(model)
}



#' @title rSPDE anisotropic inlabru mapper
#' @name bru_get_mapper.inla_rspde_anisotropic2d
#' @param model An `inla_rspde_anisotropic2d` object for which to construct or extract a mapper
#' @param \dots Arguments passed on to other methods
#' @rdname bru_get_mapper.inla_rspde_anisotropic2d
#' @rawNamespace if (getRversion() >= "3.6.0") {
#'   S3method(inlabru::bru_get_mapper, inla_rspde_anisotropic2d)
#' }
#'

bru_get_mapper.inla_rspde_anisotropic2d <- function(model, ...) {
  stopifnot(requireNamespace("inlabru"))
  inlabru_version <- as.character(packageVersion("inlabru"))
  if(inlabru_version >= "2.11.1.9022"){
      n_rep <- model[["rspde.order"]] + 1
      if((model[["est_nu"]] == 0L) && (model[["nu"]] %% 1 == 0)){
          n_rep <- 1
      }
    inlabru::bru_mapper_repeat(inlabru::bru_mapper(model[["mesh"]]), n_rep = n_rep)
  } else{
    mapper <- list(model = model)
    inlabru::bru_mapper_define(mapper, new_class = "bru_mapper_inla_rspde")
  }
}



#####################################################
#####################################################
###### HELPER FUNCTIONS FOR CGENERIC DEBUGGING ######
#####################################################
#####################################################

# prepare_cgeneric_data <- function(data) {
#   # Preprocess integers
#   if (!is.null(data$ints)) {
#     data$ints <- lapply(data$ints, function(int_vec) {
#       if (length(int_vec) == 0 || is.null(int_vec)) {
#         as.integer(NA) # Replace empty or NULL values with NA
#       } else {
#         as.integer(int_vec)
#       }
#     })
#   } else {
#     data$ints <- list() # Ensure ints is at least an empty list
#   }

#   # Preprocess doubles
#   if (!is.null(data$doubles)) {
#     data$doubles <- lapply(data$doubles, function(double_vec) {
#       if (length(double_vec) == 0 || is.null(double_vec)) {
#         as.numeric(NA) # Replace empty or NULL values with NA
#       } else {
#         as.numeric(double_vec)
#       }
#     })
#   } else {
#     data$doubles <- list() # Ensure doubles is at least an empty list
#   }

#   # Preprocess sparse matrices
#   if (!is.null(data$smatrices)) {
#     data$smatrices <- lapply(data$smatrices, function(sm) {
#       nrow <- sm[1]
#       ncol <- sm[2]
#       n <- sm[3]
#       if (n == 0) {
#         return(NULL) # Handle completely empty sparse matrices
#       }
#       i <- sm[4:(4 + n - 1)]
#       j <- sm[(4 + n):(4 + 2 * n - 1)]
#       x <- sm[(4 + 2 * n):(4 + 3 * n - 1)]
#       list(
#         nrow = as.integer(nrow),
#         ncol = as.integer(ncol),
#         n = as.integer(n),
#         i = as.integer(i),
#         j = as.integer(j),
#         x = as.numeric(x)
#       )
#     })
#   } else {
#     data$smatrices <- list() # Ensure smatrices is at least an empty list
#   }

#   # Preprocess dense matrices
#   if (!is.null(data$matrices)) {
#     data$matrices <- lapply(data$matrices, function(mat) {
#       if (length(mat) < 2) {
#         return(NULL) # Handle completely empty matrices
#       }
#       nrow <- mat[1]
#       ncol <- mat[2]
#       x <- mat[-(1:2)] # Exclude the first two elements (nrow, ncol)
#       list(
#         nrow = as.integer(nrow),
#         ncol = as.integer(ncol),
#         x = as.numeric(x)
#       )
#     })
#   } else {
#     data$matrices <- list() # Ensure matrices is at least an empty list
#   }

#   return(data)
# }

# call_dynamic_cgeneric_model <- function(model, cmd, theta) {
#   # Validate model structure
#   if (!is.list(model) || is.null(model$f) || is.null(model$f$cgeneric)) {
#     stop("The model must be a list containing an $f$cgeneric structure.")
#   }

#   cgeneric <- model$f$cgeneric
#   if (!is.list(cgeneric) || is.null(cgeneric$shlib) || is.null(cgeneric$model) || is.null(cgeneric$data)) {
#     stop("The model$f$cgeneric structure must contain 'shlib', 'model', and 'data' elements.")
#   }

#   # Load the shared library explicitly
#   so_path <- cgeneric$shlib
#   if (!file.exists(so_path)) {
#     stop(sprintf("Shared library '%s' not found.", so_path))
#   }

#   tryCatch({
#     dyn.load(so_path)
#   }, error = function(e) {
#     stop(sprintf("Failed to load shared library '%s': %s", so_path, e$message))
#   })

#   # Preprocess data
#   cgeneric$data <- prepare_cgeneric_data(cgeneric$data)

#   # Convert cmd to integer
#   cmd <- switch(
#     cmd,
#     "void" = 0,
#     "Q" = 1,
#     "graph" = 2,
#     "mu" = 3,
#     "initial" = 4,
#     "log_norm_const" = 5,
#     "log_prior" = 6,
#     "quit" = 7,
#     stop(sprintf("Invalid command '%s'", cmd))
#   )

#   # Validate theta
#   if (!is.numeric(theta) || length(theta) == 0) {
#     stop("The 'theta' argument must be a non-empty numeric vector.")
#   }

#   # Call the C function
#   result <- tryCatch(
#     {
#       .Call(
#         "call_dynamic_inla_cgeneric", # Call the C wrapper function
#         as.integer(cmd),
#         as.numeric(theta),
#         cgeneric$data,
#         cgeneric$model,
#         cgeneric$shlib
#       )
#     },
#     error = function(e) {
#       stop(
#         sprintf(
#           "Error calling function '%s' in shared library '%s': %s",
#           cgeneric$model,
#           cgeneric$shlib,
#           e$message
#         )
#       )
#     }
#   )

#   # Optionally unload the shared library after use
#   tryCatch({
#     dyn.unload(so_path)
#   }, error = function(e) {
#     warning(sprintf("Failed to unload shared library '%s': %s", so_path, e$message))
#   })

#   return(result)
# }

# build_Q_cgeneric <- function(Q, graph) {
#   # Extract nrow (and ncol, since it's square) from the first element of graph
#   n <- graph[1]
  
#   # Extract the number of nonzero elements from the second element of graph
#   num_nonzero <- graph[2]
  
#   # Extract the row (i) and column (j) indices from the graph
#   i <- graph[3:(2 + num_nonzero)]
#   j <- graph[(3 + num_nonzero):(2 + 2 * num_nonzero)]
  
#   # Extract the nonzero values from Q (skipping the first two elements)
#   Q_values <- Q[-c(1, 2)]
  
#   # Ensure the lengths match the number of nonzero elements
#   if (length(Q_values) != num_nonzero) {
#     stop("Mismatch between the number of nonzero entries and the provided values in Q.")
#   }
#   if (length(i) != num_nonzero || length(j) != num_nonzero) {
#     stop("Mismatch between the graph indices and the number of nonzero entries.")
#   }
  
#   # Create a sparse matrix using the Matrix package
#   library(Matrix)
#   Q_sparse <- sparseMatrix(
#     i = i+1,
#     j = j+1,
#     x = Q_values,
#     dims = c(n, n),
#     symmetric = TRUE  # Assuming Q is symmetric
#   )
  
#   # Return the full dense matrix (or keep it sparse as needed)
#   Q_full <- as.matrix(Q_sparse)
#   return(Q_full)
# }

# library(fmesher)
# n_loc <- 20
# loc_2d_mesh <- matrix(runif(n_loc * 2), n_loc, 2)
# mesh_2d <- fm_mesh_2d(
#     loc = loc_2d_mesh,
#     cutoff = 0.01,
#     max.edge = c(0.1, 0.5)
# )
# library(inlabru)
# library(rSPDE)
# model_aniso <- rspde.anistropic2d(mesh = mesh_2d)


# Q <- call_dynamic_cgeneric_model(model_aniso, cmd = "Q", theta = c(-1,-1,0,0,0))
# graph <- call_dynamic_cgeneric_model(model_aniso, cmd = "graph", theta = c(-1,-1,0,0,0))

# Q_full <- build_Q_cgeneric(Q, graph)
