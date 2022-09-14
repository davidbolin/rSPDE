
#'
#' @title rSPDE inlabru mapper
#' @name bru_mapper.inla.rspde
#' @param model An `inla_rspde` object to use as a mapper
#' @rdname bru_mapper.inla_rspde
#' @rawNamespace if (getRversion() >= "3.6.0") {
#'   S3method(inlabru::bru_mapper, inla_rspde)
#'   # S3method(inlabru::bru_get_mapper, inla_rspde)
#'   S3method(inlabru::ibm_n, bru_mapper_inla_rspde)
#'   S3method(inlabru::ibm_values, bru_mapper_inla_rspde)
#'   S3method(inlabru::ibm_amatrix, bru_mapper_inla_rspde)
#' }
bru_mapper.inla_rspde <- function(model,...) {
  mapper <- list(model = model)
  # Note 1: From inlabru > 2.5.3, use bru_mapper_define instead.
  # Note 2: bru_mapper.default is not exported from inlabru, so
  # must call the generic bru_mapper()
  inlabru::bru_mapper(mapper, new_class = "bru_mapper_inla_rspde")
}

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

#' @param x The model to be passed to obtain the mapper.
#' @rdname bru_mapper.inla_rspde
bru_get_mapper.inla_rspde <- function(model, ...){
 inlabru::bru_mapper(model)
}
