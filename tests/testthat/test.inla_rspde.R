library(INLA)

test_that("Checking equality of optimized and non-opt. Prec. matrices for non-integer case", {
set.seed(1)
n <- 10

coords <- cbind(long=sample(1:n), lat=sample(1:n))

prmesh <- inla.mesh.2d(coords, max.edge = c(20, 40))

rspde_order_m = 2

rspde_model <- create_rspde_model(inla_mesh = prmesh, rspde_order = rspde_order_m,
                                  optimize=FALSE, debug=FALSE, 
                                  nu_upper_bound = 2)

prec_m <- rspde_model$f$rgeneric$definition(cmd="Q", theta=c(1,1,1))

rspde_model_opt <- create_rspde_model(inla_mesh = prmesh, rspde_order = rspde_order_m,
                                  optimize=TRUE, sharp = FALSE,debug=FALSE, 
                                  nu_upper_bound = 2)


prec_opt_values <- rspde_model_opt$f$rgeneric$definition(cmd="Q", theta=c(1,1,1))
prec_opt_graph <- rspde_model_opt$f$rgeneric$definition(cmd="graph")
prec_opt <- build_sparse_matrix_rspde(prec_opt_values,prec_opt_graph)

rspde_model_opt_2 <- create_rspde_model(inla_mesh = prmesh, rspde_order = rspde_order_m,
                                      optimize=TRUE, sharp = TRUE,debug=FALSE, 
                                      nu_upper_bound = 2)
prec_opt_values_2 <- rspde_model_opt_2$f$rgeneric$definition(cmd="Q", theta=c(1,1,1))
prec_opt_graph_2 <- rspde_model_opt_2$f$rgeneric$definition(cmd="graph")
prec_opt_2 <- build_sparse_matrix_rspde(prec_opt_values_2,prec_opt_graph_2)

expect_true(all(prec_m==prec_opt))
expect_true(all(prec_m==prec_opt_2))
})

test_that("Checking equality of optimized and non-opt. Prec. matrices for integer case", {

  set.seed(1)
  n <- 10
  
  coords <- cbind(long=sample(1:n), lat=sample(1:n))
  
  prmesh <- inla.mesh.2d(coords, max.edge = c(20, 40))
  rspde_model_int <- create_rspde_model(inla_mesh = prmesh, rspde_order = rspde_order_m,
                                  optimize=FALSE, debug=FALSE, 
                                  nu = 1)

  rspde_model_int_opt <- create_rspde_model(inla_mesh = prmesh, rspde_order = rspde_order_m,
                                      optimize=TRUE, debug=FALSE, 
                                      nu = 1)

  prec_int <- rspde_model_int$f$rgeneric$definition(cmd="Q", theta=c(1,1))
  prec_int_opt_values <- rspde_model_int_opt$f$rgeneric$definition(cmd="Q", theta=c(1,1))
  prec_int_opt_graph <- rspde_model_int_opt$f$rgeneric$definition(cmd="graph")
  prec_int_opt <- build_sparse_matrix_rspde(prec_int_opt_values,prec_int_opt_graph)

  expect_true(all(prec_int==prec_int_opt))
})