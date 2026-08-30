## Minimal demo for rSPDE - five vignette-based INLA examples
options(width = 80)
set.seed(1)
cat("\n=== rSPDE minimal demo ===\n")

## ---- Example 1: rspde.matern (spatial, mesh-based) ------------------------------------------
cat("\n-- Example 1: rspde.matern (spatial) --\n")
if (rspde_safe_inla(quietly = TRUE)) try({
  library(INLA); library(rSPDE); library(fmesher)

  coords   <- matrix(runif(40), ncol = 2)
  bnd      <- fm_nonconvex_hull(coords, -0.03, -0.05, resolution = c(30, 30))
  mesh     <- fm_mesh_2d(boundary = bnd, max.edge = c(0.3, 0.6), cutoff = 0.05)

  rspde_model <- rspde.matern(mesh = mesh)
  Abar        <- rspde.make.A(mesh = mesh, loc = coords)
  idx         <- rspde.make.index(name = "field", mesh = mesh)

  Y  <- as.vector(Abar %*% rnorm(ncol(Abar)) + 0.05 * rnorm(nrow(Abar)))

  stk <- inla.stack(data = list(y = Y), A = list(Abar), effects = list(idx), tag = "est")
  res <- inla(y ~ -1 + f(field, model = rspde_model),
              data = inla.stack.data(stk),
              control.predictor = list(A = inla.stack.A(stk), compute = TRUE))
  print(summary(res))
})

## ---- Example 2: rspde.spacetime -----------------------------------------------------------
cat("\n-- Example 2: rspde.spacetime --\n")
if (rspde_safe_inla(quietly = TRUE)) try({
  library(rSPDE); library(INLA)

  s <- seq(0, 20, length.out = 21)
  t <- seq(0,  10, length.out = 6)

  op    <- spacetime.operators(space_loc = s, time_loc = t,
                               kappa = 1, sigma = 1, alpha = 1, beta = 1,
                               rho = 1, gamma = 0.1)
  model <- rspde.spacetime(space_loc = s, time_loc = t, alpha = 1, beta = 1,
                           prior.kappa = list(mean = 1), prior.sigma = list(mean = 1),
                           prior.gamma = list(mean = 0.1), prior.rho = list(mean = 1),
                           drift = TRUE)

  n   <- 80
  loc <- data.frame(x = max(s) * runif(n), t = max(t) * runif(n))
  A   <- make_A(op, loc = loc$x, time = loc$t)
  Y   <- as.vector(A %*% simulate(op, nsim = 1) + 0.05 * rnorm(n))  # spacetime op dim matches A

  stk <- inla.stack(data = list(y = Y), A = list(A),
                    effects = list(list(field = seq_len(ncol(A)))), tag = "est")
  res <- inla(y ~ -1 + f(field, model = model),
              data = inla.stack.data(stk),
              control.predictor = list(A = inla.stack.A(stk), compute = TRUE))
  print(summary(res))
})

## ---- Example 3: rspde.intrinsic ------------------------------------------------------
cat("\n-- Example 3: rspde.intrinsic --\n")
if (rspde_safe_inla(quietly = TRUE)) try({
  library(rSPDE); library(INLA); library(fmesher)

  bnd  <- fm_segm(rbind(c(0,0), c(1,0), c(1,1), c(0,1)), is.bnd = TRUE)
  mesh <- fm_mesh_2d(boundary = bnd, cutoff = 0.05, max.edge = 0.15)

  op    <- intrinsic.operators(tau = 0.2, beta = 1.8, mesh = mesh, m = 1)
  model <- rspde.intrinsic(mesh = mesh, rspde.order = 1)

  loc <- matrix(runif(80 * 2), 80, 2)
  A   <- make_A(op, loc = loc)
  Y   <- as.vector(A %*% rnorm(ncol(A)) + 0.1 * rnorm(nrow(A)))

  stk <- inla.stack(data = list(y = Y), A = list(A),
                    effects = list(list(field = seq_len(ncol(A)))))
  res <- inla(y ~ -1 + f(field, model = model),
              data = inla.stack.data(stk), family = "gaussian",
              control.predictor = list(A = inla.stack.A(stk)))
  print(summary(res))
})

## ---- Example 4: rspde.anistropic2d --------------------------------------
cat("\n-- Example 4: rspde.anistropic2d (anisotropic) --\n")
if (rspde_safe_inla(quietly = TRUE)) try({
  library(rSPDE); library(fmesher); library(INLA)

  loc  <- matrix(runif(30 * 2), 30, 2)
  mesh <- fm_mesh_2d(loc = loc, cutoff = 0.1, max.edge = c(0.25, 0.5))

  op    <- matern2d.operators(hx = 0.08, hy = 0.08, hxy = 0, nu = 0.5, sigma = 1, mesh = mesh)
  model <- rspde.anistropic2d(mesh = mesh,
                               prior.hx = list(mean = 0.08), prior.hy = list(mean = 0.08),
                               prior.hxy = list(mean = 0), prior.sigma = list(mean = 1))

  obs <- cbind(runif(50), runif(50))
  A   <- make_A(op, loc = obs)
  idx <- rspde.make.index(name = "field", mesh = mesh)
  Y   <- as.vector(A %*% rnorm(ncol(A)) + 0.1 * rnorm(nrow(A)))

  stk <- inla.stack(data = list(y = Y), A = list(A), effects = list(idx))
  res <- inla(y ~ -1 + f(field, model = model),
              data = inla.stack.data(stk),
              control.predictor = list(A = inla.stack.A(stk), compute = TRUE))
  print(summary(res))
})

## ---- Example 5: rspde.matern1d (1-D, no FEM) -----------------------------
cat("\n-- Example 5: rspde.matern1d (1-D rational approximation) --\n")
if (rspde_safe_inla(quietly = TRUE)) try({
  library(rSPDE); library(INLA)

  s     <- seq(0, 1, length.out = 101)
  model <- rspde.matern1d(loc = s, nu = 0.5)

  op_cov <- matern.rational(loc = s, nu = 0.5, range = 0.2, sigma = 1,
                             m = 1, parameterization = "matern")
  u <- as.vector(simulate(op_cov, nsim = 1))

  n.obs   <- 80
  obs.loc <- sort(runif(n.obs))
  idx     <- sapply(obs.loc, function(x) which.min(abs(s - x)))
  Aobs    <- model$A[idx, , drop = FALSE]
  Y       <- as.vector(Aobs %*% u + 0.05 * rnorm(n.obs))

  stk <- inla.stack(data = list(y = Y), A = list(Aobs),
                    effects = list(model$index))
  res <- inla(y ~ -1 + f(field, model = model),
              data = inla.stack.data(stk), family = "gaussian",
              control.predictor = list(A = inla.stack.A(stk)))
  print(summary(res))
})

cat("\n=== rSPDE minimal demo finished ===\n")
