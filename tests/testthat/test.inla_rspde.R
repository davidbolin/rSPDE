test_that("Checking equality of optimized and non-opt. Prec. matrices for non-integer case", {
  testthat::skip_on_cran()
  old_threads <- INLA::inla.getOption("num.threads")
  INLA::inla.setOption(num.threads = "1:1")
set.seed(1)
n <- 10

coords <- cbind(long=sample(1:n), lat=sample(1:n))

mesh <- INLA::inla.mesh.2d(coords, max.edge = c(20, 40))

rspde_order_m = 2

rspde_model <- rspde.matern(mesh = mesh, rspde_order = rspde_order_m,
                                  optimize=FALSE, debug=FALSE, 
                                  nu_upper_bound = 2)

prec_m <- rspde_model$f$rgeneric$definition(cmd="Q", theta=c(1,1,1))

rspde_model_opt <- rspde.matern(mesh = mesh, rspde_order = rspde_order_m,
                                  optimize=TRUE, sharp = FALSE,debug=FALSE, 
                                  nu_upper_bound = 2)


prec_opt_values <- rspde_model_opt$f$rgeneric$definition(cmd="Q", theta=c(1,1,1))
prec_opt_graph <- rspde_model_opt$f$rgeneric$definition(cmd="graph")
prec_opt <- build_sparse_matrix_rspde(prec_opt_values,prec_opt_graph)

rspde_model_opt_2 <- rspde.matern(mesh = mesh, rspde_order = rspde_order_m,
                                      optimize=TRUE, sharp = TRUE,debug=FALSE, 
                                      nu_upper_bound = 2)
prec_opt_values_2 <- rspde_model_opt_2$f$rgeneric$definition(cmd="Q", theta=c(1,1,1))
prec_opt_graph_2 <- rspde_model_opt_2$f$rgeneric$definition(cmd="graph")
prec_opt_2 <- build_sparse_matrix_rspde(prec_opt_values_2,prec_opt_graph_2)

expect_true(all(prec_m==prec_opt))
expect_true(all(prec_m==prec_opt_2))
INLA::inla.setOption(num.threads = old_threads)
})

test_that("Checking equality of optimized and non-opt. Prec. matrices for integer case", {
  testthat::skip_on_cran()
  old_threads <- INLA::inla.getOption("num.threads")
  INLA::inla.setOption(num.threads = "1:1")
  set.seed(1)
  n <- 10
  
  coords <- cbind(long=sample(1:n), lat=sample(1:n))
  
  mesh <- INLA::inla.mesh.2d(coords, max.edge = c(20, 40))
  rspde_model_int <- rspde.matern(mesh = mesh, rspde_order = rspde_order_m,
                                  optimize=FALSE, debug=FALSE, 
                                  nu = 1)

  rspde_model_int_opt <- rspde.matern(mesh = mesh, rspde_order = rspde_order_m,
                                      optimize=TRUE, debug=FALSE, 
                                      nu = 1)

  prec_int <- rspde_model_int$f$rgeneric$definition(cmd="Q", theta=c(1,1))
  prec_int_opt_values <- rspde_model_int_opt$f$rgeneric$definition(cmd="Q", theta=c(1,1))
  prec_int_opt_graph <- rspde_model_int_opt$f$rgeneric$definition(cmd="graph")
  prec_int_opt <- build_sparse_matrix_rspde(prec_int_opt_values,prec_int_opt_graph)
  
  expect_true(all(prec_int==prec_int_opt))
  
  spde <- INLA::inla.spde2.matern(mesh,alpha=2)
  prec_temp <- INLA::inla.spde.precision(spde,theta=c(1,3))
  prec_int <- rspde.precision(rspde_model_int, theta=c(1,3))
  expect_true(sum((prec_temp-prec_int)^2) < 10^(-10))
  INLA::inla.setOption(num.threads = old_threads)
})

test_that("Checking equality of optimized and non-opt Prec. matrices for d=1",{
  
  testthat::skip_on_cran()
  
  old_threads <- INLA::inla.getOption("num.threads")
  INLA::inla.setOption(num.threads = "1:1")
  kappa <- 20
  sigma <- 1
  nu <- 0.1
  nu_upper_bound = 2
  
  #create mass and stiffness matrices for a FEM discretization
  nobs = 101
  x <- seq(from = 0, to = 1, length.out = 101)
  mesh <- INLA::inla.mesh.1d(x)
  fem <- rSPDE.fem1d(x)
  
  #compute rational approximation of covariance function at 0.5
  tau <- sqrt(gamma(nu) / (sigma^2 * kappa^(2*nu) * (4*pi)^(1/2) * gamma(nu+1/2)))

  op_cov <- matern.operators(C=fem$C, G=fem$G,nu=nu,
                                     kappa=kappa,sigma=sigma,d=1,m=2)
  
  rspde_model.opt <- rspde.matern(mesh = mesh, rspde_order = 2,
                                  optimize=TRUE,
                              nu_upper_bound = nu_upper_bound)
  
  #Compute the precision matrix
  Q <- rspde.matern.precision(kappa=kappa,nu=nu,tau=tau,
                              rspde_order=2,d=1,fem_mesh_matrices = op_cov$fem_mesh_matrices)

  nu_upper_bound = rspde_model.opt$nu_upper_bound
  
  Q.opt <- rspde.precision(rspde_model.opt, theta=c(log(tau),log(kappa),log(nu/(nu_upper_bound-nu))))
  
  expect_true(sum((Q-Q.opt)^2) < 10^(-10))
  
  rspde_model <- rspde.matern(mesh = mesh, rspde_order = 2,
                                  optimize=FALSE,
                                  nu_upper_bound = nu_upper_bound)
  
  Q2 <- rspde.precision(rspde_model, theta=c(log(tau),log(kappa),log(nu/(nu_upper_bound-nu))))
  
  expect_true(sum((Q2-Q.opt)^2) < 10^(-10))
  
  INLA::inla.setOption(num.threads = old_threads)
  
})