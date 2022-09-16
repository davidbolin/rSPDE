
#'
#' @title rSPDE inlabru mapper
#' @name bru_mapper.inla_rspde
#' @param model An `inla_rspde` for which to construct or extract a mapper
#' @param \dots Arguments passed on to other methods
#' @rdname bru_mapper.inla_rspde
#' @rawNamespace if (getRversion() >= "3.6.0") {
#'   S3method(inlabru::bru_mapper, inla_rspde)
#'   # S3method(inlabru::bru_get_mapper, inla_rspde)
#'   S3method(inlabru::ibm_n, bru_mapper_inla_rspde)
#'   S3method(inlabru::ibm_values, bru_mapper_inla_rspde)
#'   S3method(inlabru::ibm_amatrix, bru_mapper_inla_rspde)
#' }
#' 
#' @examples
#' \donttest{ #tryCatch version
#' tryCatch({
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
#' # For inlabru version 2.5.3.9002 or above:
#' cmp <- y ~ Intercept(1) + 
#'            field(coordinates, model = rspde_model)
#' 
#' #For inlabru version 2.5.3:
#' cmp <- y ~ Intercept(1) + 
#'            field(coordinates, model = rspde_model,
#'              mapper = bru_mapper(rspde_model))
#' 
#' rspde_fit <- bru(cmp, data = data_df)
#' }
#' #stable.tryCatch
#' }, error = function(e){print("Could not run the example")})
#' }
bru_mapper.inla_rspde <- function(model,...) {
  mapper <- list(model = model)
  # Note 1: From inlabru > 2.5.3, use bru_mapper_define instead.
  # Note 2: bru_mapper.default is not exported from inlabru, so
  # must call the generic bru_mapper()
  inlabru::bru_mapper(mapper, new_class = "bru_mapper_inla_rspde")
}

#' @param mapper A `bru_mapper.inla_rspde` object
#' @rdname bru_mapper.inla_rspde
ibm_n.bru_mapper_inla_rspde <- function(mapper, ...) {
  model <- mapper[["model"]]
  integer_nu <- model$integer.nu
  rspde_order <- model$rspde_order
  if(integer_nu){
            factor_rspde <- 1
  } else{
            factor_rspde <- rspde_order + 1
  }
  factor_rspde*model$n.spde
}
#' @rdname bru_mapper.inla_rspde
ibm_values.bru_mapper_inla_rspde <- function(mapper, ...) {
  seq_len(inlabru::ibm_n(mapper))
}
#' @param input The values for which to produce a mapping matrix
#' @rdname bru_mapper.inla_rspde
ibm_amatrix.bru_mapper_inla_rspde <- function(mapper, input, ...) {
  if (is.null(input)) {
    return(Matrix::Matrix(0, 0, inlabru::ibm_n(mapper)))
  }
  if (!is.matrix(input) && !inherits(input, "Spatial")) {
    input <- as.matrix(input)
  }
  model <- mapper[["model"]]
  if(model$est_nu){
    nu <- NULL
  } else{
   nu <- model$nu
  }
  rspde_order <- model$rspde_order
  rSPDE::rspde.make.A(mesh = model$mesh, loc=input,
                                rspde_order = rspde_order,
                                nu=nu)
}

#' @rdname bru_mapper.inla_rspde
bru_get_mapper.inla_rspde <- function(model, ...){
 inlabru::bru_mapper(model)
}
