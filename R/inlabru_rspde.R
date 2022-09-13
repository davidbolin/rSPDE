
#' 
#' @title rSPDE inlabru mapper
#' @name bru_mapper.inla.rspde
#' @param model An `inla_rspde` object to use as a mapper
#' @rdname bru_mapper.inla_rspde

bru_mapper.inla_rspde <- function(model,...) {
  mapper <- list(model = model)
  class(mapper) <- c("bru_mapper_inla_rspde", "list")
  bru_mapper(mapper)
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
  model <- mapper[["model"]]
  integer_nu <- model$integer.nu
  rspde_order <- model$rspde_order
  if(integer_nu){
            factor_rspde <- 1
  } else{
            factor_rspde <- rspde_order + 1
  }
  seq_len(factor_rspde*model$n.spde)
}
#' @param input The values for which to produce a mapping matrix
#' @rdname bru_mapper.inla_rspde
ibm_amatrix.bru_mapper_inla_rspde <- function(mapper, input, ...) {
  if (is.null(input)) {
    return(Matrix::Matrix(0, 0, ibm_n(mapper)))
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
 bru_mapper(model) 
}
