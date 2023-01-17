context("inla_rspde")

test_that("testing cgeneric_integer", {

  testthat::skip_on_cran()
  inla_installed <- "INLA" %in% rownames(installed.packages())
  if(!inla_installed){
    testthat::skip("INLA not installed")
  }
  
  old_threads <- INLA::inla.getOption("num.threads")
  INLA::inla.setOption(num.threads = "1:1")


data(PRprec, package = "INLA")

Y <- rowMeans(PRprec[, 3 + 1:31])

ind <- !is.na(Y)
Y <- Y[ind]
coords <- as.matrix(PRprec[ind, 1:2])

prdomain <- INLA::inla.nonconvex.hull(coords, -0.03, -0.05, resolution = c(50, 50))
prmesh <- INLA::inla.mesh.2d(boundary = prdomain, max.edge = c(0.45, 1), cutoff = 0.5)


rspde_model <- rspde.matern(mesh = prmesh, parameterization = "spde", 
nu = 1)

Q_tmp <- INLA::inla.cgeneric.q(rspde_model)

inla_model <- INLA::inla.spde2.matern(
    mesh = prmesh, alpha = 2
)

Q_1 <- INLA::inla.spde.precision(
    inla_model, theta = Q_tmp$theta
)

testthat::expect_equal(sum( (Q_1 - Q_tmp$Q)^2), 0)

  INLA::inla.setOption(num.threads = old_threads)
})


test_that("testing cgeneric_parsimonious_fixed", {

  testthat::skip_on_cran()
  inla_installed <- "INLA" %in% rownames(installed.packages())
  if(!inla_installed){
    testthat::skip("INLA not installed")
  }
  
  old_threads <- INLA::inla.getOption("num.threads")
  INLA::inla.setOption(num.threads = "1:1")



data(PRprec, package = "INLA")

Y <- rowMeans(PRprec[, 3 + 1:31])

ind <- !is.na(Y)
Y <- Y[ind]
coords <- as.matrix(PRprec[ind, 1:2])

prdomain <- INLA::inla.nonconvex.hull(coords, -0.03, -0.05, resolution = c(50, 50))
prmesh <- INLA::inla.mesh.2d(boundary = prdomain, max.edge = c(0.45, 1), cutoff = 0.5)


rspde_model <- rspde.matern(mesh = prmesh, parameterization = "spde", 
nu = 0.4,
rspde.order = 0)

Q_tmp <- INLA::inla.cgeneric.q(rspde_model)

inla_model <- INLA::inla.spde2.matern(
    mesh = prmesh, alpha = 1.4
)

Q_1 <- INLA::inla.spde.precision(
    inla_model, theta = Q_tmp$theta
)

testthat::expect_equal(sum( (Q_1 - Q_tmp$Q)^2), 0)

  INLA::inla.setOption(num.threads = old_threads)
})


test_that("testing cgeneric_parsimonious_gen", {

  testthat::skip_on_cran()
  inla_installed <- "INLA" %in% rownames(installed.packages())
  if(!inla_installed){
    testthat::skip("INLA not installed")
  }
  
  old_threads <- INLA::inla.getOption("num.threads")
  INLA::inla.setOption(num.threads = "1:1")


data(PRprec, package = "INLA")

Y <- rowMeans(PRprec[, 3 + 1:31])

ind <- !is.na(Y)
Y <- Y[ind]
coords <- as.matrix(PRprec[ind, 1:2])

prdomain <- INLA::inla.nonconvex.hull(coords, -0.03, -0.05, resolution = c(50, 50))
prmesh <- INLA::inla.mesh.2d(boundary = prdomain, max.edge = c(0.45, 1), cutoff = 0.5)


rspde_model <- rspde.matern(mesh = prmesh, parameterization = "spde", 
start.nu = 0.4,
rspde.order = 0)

Q_tmp <- INLA::inla.cgeneric.q(rspde_model)

inla_model <- INLA::inla.spde2.matern(
    mesh = prmesh, alpha = 1.4
)

Q_1 <- INLA::inla.spde.precision(
    inla_model, theta = Q_tmp$theta[1:2]
)

testthat::expect_equal(sum( (Q_1 - Q_tmp$Q)^2), 0)

  INLA::inla.setOption(num.threads = old_threads)
})


test_that("testing cgeneric_rspde_fixed_gen", {

  testthat::skip_on_cran()
  inla_installed <- "INLA" %in% rownames(installed.packages())
  if(!inla_installed){
    testthat::skip("INLA not installed")
  }
  
  old_threads <- INLA::inla.getOption("num.threads")
  INLA::inla.setOption(num.threads = "1:1")


data(PRprec, package = "INLA")

Y <- rowMeans(PRprec[, 3 + 1:31])

ind <- !is.na(Y)
Y <- Y[ind]
coords <- as.matrix(PRprec[ind, 1:2])

prdomain <- INLA::inla.nonconvex.hull(coords, -0.03, -0.05, resolution = c(50, 50))
prmesh <- INLA::inla.mesh.2d(boundary = prdomain, max.edge = c(0.45, 1), cutoff = 0.5)


rspde_model <- rspde.matern(mesh = prmesh, parameterization = "spde", 
start.nu = 0.4,
rspde.order = 2)

Q_tmp <- INLA::inla.cgeneric.q(rspde_model)

rspde_model_fixed <- rspde.matern(mesh = prmesh, parameterization = "spde", 
nu = 0.4,
rspde.order = 2)

Q_tmp2 <- INLA::inla.cgeneric.q(rspde_model_fixed)

testthat::expect_equal(sum( (Q_tmp2$Q - Q_tmp$Q)^2), 0)

  INLA::inla.setOption(num.threads = old_threads)
})

test_that("testing cgeneric_rspde_gen", {

  testthat::skip_on_cran()
  inla_installed <- "INLA" %in% rownames(installed.packages())
  if(!inla_installed){
    testthat::skip("INLA not installed")
  }
  
  old_threads <- INLA::inla.getOption("num.threads")
  INLA::inla.setOption(num.threads = "1:1")


data(PRprec, package = "INLA")

Y <- rowMeans(PRprec[, 3 + 1:31])

ind <- !is.na(Y)
Y <- Y[ind]
coords <- as.matrix(PRprec[ind, 1:2])

prdomain <- INLA::inla.nonconvex.hull(coords, -0.03, -0.05, resolution = c(50, 50))
prmesh <- INLA::inla.mesh.2d(boundary = prdomain, max.edge = c(0.45, 1), cutoff = 0.5)


rspde_model <- rspde.matern(mesh = prmesh, parameterization = "spde", 
start.nu = 0.4,
rspde.order = 2)

Q_tmp <- INLA::inla.cgeneric.q(rspde_model)

C <- rspde_model$fem_mesh[["c0"]]
G <- rspde_model$fem_mesh[["g1"]]

op <- matern.operators(kappa = exp(Q_tmp$theta[2]),
                        nu = 0.4,
                        tau = exp(Q_tmp$theta[1]),
                        C = C, G = G, d = 2, m = 2)

Q_1 <- precision(op)

testthat::expect_equal(sum( (Q_1 - Q_tmp$Q)^2), 0)

  INLA::inla.setOption(num.threads = old_threads)
})


test_that("testing cgeneric_nonstat_gen", {

  testthat::skip_on_cran()
  inla_installed <- "INLA" %in% rownames(installed.packages())
  if(!inla_installed){
    testthat::skip("INLA not installed")
  }
  
  old_threads <- INLA::inla.getOption("num.threads")
  INLA::inla.setOption(num.threads = "1:1")

data(PRprec, package = "INLA")

Y <- rowMeans(PRprec[, 3 + 1:31])

ind <- !is.na(Y)
Y <- Y[ind]
coords <- as.matrix(PRprec[ind, 1:2])

prdomain <- INLA::inla.nonconvex.hull(coords, -0.03, -0.05, resolution = c(50, 50))
prmesh <- INLA::inla.mesh.2d(boundary = prdomain, max.edge = c(0.45, 1), cutoff = 0.5)


rspde_model <- rspde.matern(mesh = prmesh, parameterization = "spde",
                       B.tau = cbind(0,1,rep(0,prmesh$n)))

stopifnot(!rspde_model$stationary)

Q_tmp <- INLA::inla.cgeneric.q(rspde_model)

rspde_stat_model <- rspde.matern(mesh = prmesh,
                            parameterization = "spde")

stopifnot(rspde_stat_model$stationary)

Q_1 <- INLA::inla.cgeneric.q(rspde_stat_model)

testthat::expect_equal(sum( (Q_1$Q - Q_tmp$Q)^2), 0)

  INLA::inla.setOption(num.threads = old_threads)
})

test_that("testing cgeneric_nonstat_fixed", {

  testthat::skip_on_cran()
  inla_installed <- "INLA" %in% rownames(installed.packages())
  if(!inla_installed){
    testthat::skip("INLA not installed")
  }
  
  old_threads <- INLA::inla.getOption("num.threads")
  INLA::inla.setOption(num.threads = "1:1")

data(PRprec, package = "INLA")

Y <- rowMeans(PRprec[, 3 + 1:31])

ind <- !is.na(Y)
Y <- Y[ind]
coords <- as.matrix(PRprec[ind, 1:2])

prdomain <- INLA::inla.nonconvex.hull(coords, -0.03, -0.05, resolution = c(50, 50))
prmesh <- INLA::inla.mesh.2d(boundary = prdomain, max.edge = c(0.45, 1), cutoff = 0.5)


rspde_model <- rspde.matern(mesh = prmesh, nu = 0.7, parameterization = "spde",
                       B.tau = cbind(0,1,rep(0,prmesh$n)))

stopifnot(!rspde_model$stationary)

Q_tmp <- INLA::inla.cgeneric.q(rspde_model)

rspde_stat_model <- rspde.matern(mesh = prmesh, nu = 0.7,
                            parameterization = "spde")

stopifnot(rspde_stat_model$stationary)

Q_1 <- INLA::inla.cgeneric.q(rspde_stat_model)

testthat::expect_equal(sum( (Q_1$Q - Q_tmp$Q)^2), 0)

  INLA::inla.setOption(num.threads = old_threads)
})

test_that("testing cgeneric_nonstat_integer", {

  testthat::skip_on_cran()
  inla_installed <- "INLA" %in% rownames(installed.packages())
  if(!inla_installed){
    testthat::skip("INLA not installed")
  }
  
  old_threads <- INLA::inla.getOption("num.threads")
  INLA::inla.setOption(num.threads = "1:1")

data(PRprec, package = "INLA")

Y <- rowMeans(PRprec[, 3 + 1:31])

ind <- !is.na(Y)
Y <- Y[ind]
coords <- as.matrix(PRprec[ind, 1:2])

prdomain <- INLA::inla.nonconvex.hull(coords, -0.03, -0.05, resolution = c(50, 50))
prmesh <- INLA::inla.mesh.2d(boundary = prdomain, max.edge = c(0.45, 1), cutoff = 0.5)


rspde_model <- rspde.matern(mesh = prmesh, parameterization = "spde", nu = 1,
                       B.tau = cbind(0,1,rep(0,prmesh$n)))

stopifnot(!rspde_model$stationary)

start_theta <- rspde_model$start.theta

tau <- rep(exp(start_theta[1]), prmesh$n)

kappa <- rep(exp(start_theta[2]), prmesh$n)

ns_op <- spde.matern.operators(kappa = kappa, tau = tau, nu = 1,
                                mesh = prmesh, m = 2, type = "covariance")

Q_tmp <- INLA::inla.cgeneric.q(rspde_model)

testthat::expect_equal(sum( (Q_tmp$Q - ns_op$Q)^2), 0)

rspde_stat_model <- rspde.matern(mesh = prmesh, nu = 1,
                            parameterization = "spde")

stopifnot(rspde_stat_model$stationary)

Q_1 <- INLA::inla.cgeneric.q(rspde_stat_model)

testthat::expect_equal(sum( (Q_1$Q - Q_tmp$Q)^2), 0)

  INLA::inla.setOption(num.threads = old_threads)
})
