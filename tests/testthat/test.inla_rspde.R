context("inla_rspde")

test_that("testing cgeneric_integer", {

  testthat::skip_on_cran()
if (!requireNamespace("INLA", quietly=TRUE))
    testthat::skip(message = 'INLA package is not installed. (see www.r-inla.org/download-install)')

  old_threads <- INLA::inla.getOption("num.threads")
  INLA::inla.setOption(num.threads = "1:1")


data(PRprec, package = "INLA")

Y <- rowMeans(PRprec[, 3 + 1:31])

ind <- !is.na(Y)
Y <- Y[ind]
coords <- as.matrix(PRprec[ind, 1:2])

prdomain <- fmesher::fm_nonconvex_hull_inla(coords, -0.03, -0.05, resolution = c(50, 50))
prmesh <- fmesher::fm_mesh_2d_inla(boundary = prdomain, max.edge = c(0.45, 1), cutoff = 0.5)


rspde_model <- rspde.matern(mesh = prmesh, parameterization = "spde",
nu = 1)

Q_tmp <- INLA::inla.cgeneric.q(rspde_model)

Q_tmp2 <- precision(rspde_model)

testthat::expect_equal(sum( (Q_tmp2 - Q_tmp$Q)^2), 0)

inla_model <- INLA::inla.spde2.matern(
    mesh = prmesh, alpha = 2
)

Q_1 <- INLA::inla.spde.precision(
    inla_model, theta = Q_tmp$theta
)

testthat::expect_equal(sum( (Q_1 - Q_tmp$Q)^2), 0)

## Testing matern parameterization

rspde_model <- rspde.matern(mesh = prmesh, parameterization = "matern",
nu = 1)

Q_tmp <- INLA::inla.cgeneric.q(rspde_model)

Q_tmp2 <- precision(rspde_model)

testthat::expect_equal(sum( (Q_tmp2 - Q_tmp$Q)^2), 0)

  INLA::inla.setOption(num.threads = old_threads)
})


# test_that("testing cgeneric_parsimonious_fixed", {

#   testthat::skip_on_cran()
# if (!requireNamespace("INLA", quietly=TRUE))
#     testthat::skip(message = 'INLA package is not installed. (see www.r-inla.org/download-install)')

#   old_threads <- INLA::inla.getOption("num.threads")
#   INLA::inla.setOption(num.threads = "1:1")



# data(PRprec, package = "INLA")

# Y <- rowMeans(PRprec[, 3 + 1:31])

# ind <- !is.na(Y)
# Y <- Y[ind]
# coords <- as.matrix(PRprec[ind, 1:2])

# prdomain <- fmesher::fm_nonconvex_hull_inla(coords, -0.03, -0.05, resolution = c(50, 50))
# prmesh <- fmesher::fm_mesh_2d(boundary = prdomain, max.edge = c(0.45, 1), cutoff = 0.5)


# rspde_model <- rspde.matern(mesh = prmesh, parameterization = "spde",
# nu = 0.4,
# rspde.order = 0)

# Q_tmp <- INLA::inla.cgeneric.q(rspde_model)

# Q_tmp2 <- precision(rspde_model)

# testthat::expect_equal(sum((Q_tmp$Q - Q_tmp2)^2), 0)

# inla_model <- INLA::inla.spde2.matern(
#     mesh = prmesh, alpha = 1.4
# )

# Q_1 <- INLA::inla.spde.precision(
#     inla_model, theta = Q_tmp$theta
# )

# testthat::expect_equal(sum( (Q_1 - Q_tmp$Q)^2), 0)

# ## Testing matern parameterization

# rspde_model <- rspde.matern(mesh = prmesh, parameterization = "matern",
# nu = 0.4,
# rspde.order = 0)

# Q_tmp <- INLA::inla.cgeneric.q(rspde_model)

# Q_tmp2 <- precision(rspde_model)

# testthat::expect_equal(sum((Q_tmp$Q - Q_tmp2)^2), 0)

#   INLA::inla.setOption(num.threads = old_threads)
# })


# test_that("testing cgeneric_parsimonious_gen", {

#   testthat::skip_on_cran()
# if (!requireNamespace("INLA", quietly=TRUE))
#     testthat::skip(message = 'INLA package is not installed. (see www.r-inla.org/download-install)')

#   old_threads <- INLA::inla.getOption("num.threads")
#   INLA::inla.setOption(num.threads = "1:1")


# data(PRprec, package = "INLA")

# Y <- rowMeans(PRprec[, 3 + 1:31])

# ind <- !is.na(Y)
# Y <- Y[ind]
# coords <- as.matrix(PRprec[ind, 1:2])

# prdomain <- fmesher::fm_nonconvex_hull_inla(coords, -0.03, -0.05, resolution = c(50, 50))
# prmesh <- fmesher::fm_mesh_2d(boundary = prdomain, max.edge = c(0.45, 1), cutoff = 0.5)


# rspde_model <- rspde.matern(mesh = prmesh, parameterization = "spde",
# start.nu = 0.4,
# rspde.order = 0)

# Q_tmp <- INLA::inla.cgeneric.q(rspde_model)

# Q_tmp2 <- precision(rspde_model)

# testthat::expect_equal(sum((Q_tmp$Q - Q_tmp2)^2), 0)

# inla_model <- INLA::inla.spde2.matern(
#     mesh = prmesh, alpha = 1.4
# )

# Q_1 <- INLA::inla.spde.precision(
#     inla_model, theta = Q_tmp$theta[1:2]
# )

# testthat::expect_equal(sum( (Q_1 - Q_tmp$Q)^2), 0)

# ## Testing matern parameterization

# rspde_model <- rspde.matern(mesh = prmesh, parameterization = "matern",
# start.nu = 0.4,
# rspde.order = 0)

# Q_tmp <- INLA::inla.cgeneric.q(rspde_model)

# Q_tmp2 <- precision(rspde_model)

# testthat::expect_equal(sum((Q_tmp$Q - Q_tmp2)^2), 0)

#   INLA::inla.setOption(num.threads = old_threads)
# })


test_that("testing cgeneric_rspde_fixed_gen", {

  testthat::skip_on_cran()
if (!requireNamespace("INLA", quietly=TRUE))
    testthat::skip(message = 'INLA package is not installed. (see www.r-inla.org/download-install)')

  old_threads <- INLA::inla.getOption("num.threads")
  INLA::inla.setOption(num.threads = "1:1")


data(PRprec, package = "INLA")

Y <- rowMeans(PRprec[, 3 + 1:31])

ind <- !is.na(Y)
Y <- Y[ind]
coords <- as.matrix(PRprec[ind, 1:2])

prdomain <- fmesher::fm_nonconvex_hull_inla(coords, -0.03, -0.05, resolution = c(50, 50))
prmesh <- fmesher::fm_mesh_2d(boundary = prdomain, max.edge = c(0.45, 1), cutoff = 0.5)


rspde_model <- rspde.matern(mesh = prmesh, parameterization = "spde",
start.nu = 0.4,
rspde.order = 2)

Q_tmp <- INLA::inla.cgeneric.q(rspde_model)

Q_tmp2 <- precision(rspde_model)

testthat::expect_equal(sum((Q_tmp$Q - Q_tmp2)^2), 0)

rspde_model_fixed <- rspde.matern(mesh = prmesh, parameterization = "spde",
nu = 0.4,
rspde.order = 2)

Q_tmp2 <- INLA::inla.cgeneric.q(rspde_model_fixed)

testthat::expect_equal(sum( (Q_tmp2$Q - Q_tmp$Q)^2), 0)

## Testing matern parameterization

rspde_model <- rspde.matern(mesh = prmesh, parameterization = "matern",
start.nu = 0.4,
rspde.order = 2)

Q_tmp <- INLA::inla.cgeneric.q(rspde_model)

Q_tmp2 <- precision(rspde_model)

testthat::expect_equal(sum((Q_tmp$Q - Q_tmp2)^2), 0)

  INLA::inla.setOption(num.threads = old_threads)
})

test_that("testing cgeneric_rspde_gen", {

  testthat::skip_on_cran()
if (!requireNamespace("INLA", quietly=TRUE))
    testthat::skip(message = 'INLA package is not installed. (see www.r-inla.org/download-install)')

  old_threads <- INLA::inla.getOption("num.threads")
  INLA::inla.setOption(num.threads = "1:1")


data(PRprec, package = "INLA")

Y <- rowMeans(PRprec[, 3 + 1:31])

ind <- !is.na(Y)
Y <- Y[ind]
coords <- as.matrix(PRprec[ind, 1:2])

prdomain <- fmesher::fm_nonconvex_hull_inla(coords, -0.03, -0.05, resolution = c(50, 50))
prmesh <- fmesher::fm_mesh_2d(boundary = prdomain, max.edge = c(0.45, 1), cutoff = 0.5)


rspde_model <- rspde.matern(mesh = prmesh, parameterization = "spde",
start.nu = 0.4,
rspde.order = 2)

Q_tmp <- INLA::inla.cgeneric.q(rspde_model)

Q_1 <- precision(rspde_model)

testthat::expect_equal(sum( (Q_1 - Q_tmp$Q)^2), 0)

## Testing matern parameterization

rspde_model <- rspde.matern(mesh = prmesh, parameterization = "matern",
start.nu = 0.4,
rspde.order = 2)

Q_tmp <- INLA::inla.cgeneric.q(rspde_model)

Q_1 <- precision(rspde_model)

testthat::expect_equal(sum( (Q_1 - Q_tmp$Q)^2), 0)

  INLA::inla.setOption(num.threads = old_threads)
})


test_that("testing cgeneric_nonstat_gen", {

  testthat::skip_on_cran()
if (!requireNamespace("INLA", quietly=TRUE))
    testthat::skip(message = 'INLA package is not installed. (see www.r-inla.org/download-install)')

  old_threads <- INLA::inla.getOption("num.threads")
  INLA::inla.setOption(num.threads = "1:1")

data(PRprec, package = "INLA")

Y <- rowMeans(PRprec[, 3 + 1:31])

ind <- !is.na(Y)
Y <- Y[ind]
coords <- as.matrix(PRprec[ind, 1:2])

prdomain <- fmesher::fm_nonconvex_hull_inla(coords, -0.03, -0.05, resolution = c(50, 50))
prmesh <- fmesher::fm_mesh_2d(boundary = prdomain, max.edge = c(0.45, 1), cutoff = 0.5)


rspde_model <- rspde.matern(mesh = prmesh, parameterization = "spde",
                       B.tau = cbind(0,1,rep(0,prmesh$n)))

stopifnot(!rspde_model$stationary)

Q_tmp <- INLA::inla.cgeneric.q(rspde_model)

Q_tmp2 <- precision(rspde_model)

testthat::expect_equal(sum( (Q_tmp2 - Q_tmp$Q)^2), 0)


rspde_stat_model <- rspde.matern(mesh = prmesh,
                            parameterization = "spde")

stopifnot(rspde_stat_model$stationary)

Q_1 <- INLA::inla.cgeneric.q(rspde_stat_model)

testthat::expect_equal(sum( (Q_1$Q - Q_tmp$Q)^2), 0)

## Testing matern parameterization

rspde_model <- rspde.matern(mesh = prmesh, parameterization = "matern",
                       B.sigma = cbind(0,1,rep(0,prmesh$n)))

stopifnot(!rspde_model$stationary)

Q_tmp <- INLA::inla.cgeneric.q(rspde_model)

Q_tmp2 <- precision(rspde_model)

testthat::expect_equal(sum( (Q_tmp2 - Q_tmp$Q)^2), 0)

  INLA::inla.setOption(num.threads = old_threads)
})

test_that("testing cgeneric_nonstat_fixed", {

  testthat::skip_on_cran()
if (!requireNamespace("INLA", quietly=TRUE))
    testthat::skip(message = 'INLA package is not installed. (see www.r-inla.org/download-install)')

  old_threads <- INLA::inla.getOption("num.threads")
  INLA::inla.setOption(num.threads = "1:1")

data(PRprec, package = "INLA")

Y <- rowMeans(PRprec[, 3 + 1:31])

ind <- !is.na(Y)
Y <- Y[ind]
coords <- as.matrix(PRprec[ind, 1:2])

prdomain <- fmesher::fm_nonconvex_hull_inla(coords, -0.03, -0.05, resolution = c(50, 50))
prmesh <- fmesher::fm_mesh_2d(boundary = prdomain, max.edge = c(0.45, 1), cutoff = 0.5)


rspde_model <- rspde.matern(mesh = prmesh, nu = 0.7, parameterization = "spde",
                       B.tau = cbind(0,1,rep(0,prmesh$n)))

stopifnot(!rspde_model$stationary)

Q_tmp <- INLA::inla.cgeneric.q(rspde_model)

Q_tmp2 <- precision(rspde_model)

testthat::expect_equal(sum( (Q_tmp2 - Q_tmp$Q)^2), 0)

rspde_stat_model <- rspde.matern(mesh = prmesh, nu = 0.7,
                            parameterization = "spde")

stopifnot(rspde_stat_model$stationary)

Q_1 <- INLA::inla.cgeneric.q(rspde_stat_model)

testthat::expect_equal(sum( (Q_1$Q - Q_tmp$Q)^2), 0)

## Testing matern parameterization

rspde_model <- rspde.matern(mesh = prmesh, nu = 0.7, parameterization = "matern",
                       B.sigma = cbind(0,1,rep(0,prmesh$n)))

stopifnot(!rspde_model$stationary)

Q_tmp <- INLA::inla.cgeneric.q(rspde_model)

Q_tmp2 <- precision(rspde_model)

testthat::expect_equal(sum( (Q_tmp2 - Q_tmp$Q)^2), 0)

  INLA::inla.setOption(num.threads = old_threads)
})

test_that("testing cgeneric_nonstat_integer", {

  testthat::skip_on_cran()
if (!requireNamespace("INLA", quietly=TRUE))
    testthat::skip(message = 'INLA package is not installed. (see www.r-inla.org/download-install)')

  old_threads <- INLA::inla.getOption("num.threads")
  INLA::inla.setOption(num.threads = "1:1")

data(PRprec, package = "INLA")

Y <- rowMeans(PRprec[, 3 + 1:31])

ind <- !is.na(Y)
Y <- Y[ind]
coords <- as.matrix(PRprec[ind, 1:2])

prdomain <- fmesher::fm_nonconvex_hull_inla(coords, -0.03, -0.05, resolution = c(50, 50))
prmesh <- fmesher::fm_mesh_2d(boundary = prdomain, max.edge = c(0.45, 1), cutoff = 0.5)


rspde_model <- rspde.matern(mesh = prmesh, parameterization = "spde", nu = 1,
                       B.tau = cbind(0,1,rep(0,prmesh$n)))

stopifnot(!rspde_model$stationary)

Q_tmp <- INLA::inla.cgeneric.q(rspde_model)

Q_tmp2 <- precision(rspde_model)

testthat::expect_equal(sum( (Q_tmp$Q - Q_tmp2)^2), 0)

rspde_stat_model <- rspde.matern(mesh = prmesh, nu = 1,
                            parameterization = "spde")

stopifnot(rspde_stat_model$stationary)

Q_1 <- INLA::inla.cgeneric.q(rspde_stat_model)

testthat::expect_equal(sum( (Q_1$Q - Q_tmp$Q)^2), 0)

## Testing matern parameterization

rspde_model <- rspde.matern(mesh = prmesh, parameterization = "matern", nu = 1,
                       B.sigma = cbind(0,1,rep(0,prmesh$n)))

stopifnot(!rspde_model$stationary)

Q_tmp <- INLA::inla.cgeneric.q(rspde_model)

Q_tmp2 <- precision(rspde_model)

testthat::expect_equal(sum( (Q_tmp$Q - Q_tmp2)^2), 0)

  INLA::inla.setOption(num.threads = old_threads)
})
