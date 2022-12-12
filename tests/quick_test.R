#devtools::install_github("davidbolin/rSPDE", ref = "devel")
library(INLA)
devtools::load_all()
data(PRprec)
data(PRborder)

Y <- rowMeans(PRprec[, 3 + 1:31])

ind <- !is.na(Y)
Y <- Y[ind]
coords <- as.matrix(PRprec[ind, 1:2])

prdomain <- inla.nonconvex.hull(coords, -0.03, -0.05, resolution = c(100, 100))
prmesh <- inla.mesh.2d(boundary = prdomain, max.edge = c(0.45, 1), cutoff = 0.5)


Abar <- rspde.make.A(mesh = prmesh, loc = coords)

mesh.index <- rspde.make.index(name = "field", mesh = prmesh)

rspde_model <- rspde.matern(mesh = prmesh)

# rspde_stat <- rspde.matern(mesh = prmesh,
#                                 parameterization = "spde",
#                                 prior.nu.dist = "beta",
#                                 start.lkappa = rspde_model$param$theta.prior.mean[2],
#                                 start.ltau = rspde_model$param$theta.prior.mean[1],
#                                 prior.kappa = list(meanlog = rspde_model$param$theta.prior.mean[2]),
#                                 prior.tau = list(meanlog = rspde_model$param$theta.prior.mean[1]))

# Q_nonstat <- inla.cgeneric.q(rspde_model)
# # Q_nonstat <- Q_nonstat$Q
# Q_stat <- inla.cgeneric.q(rspde_stat)
# # Q_stat <- Q_stat$Q




stk.dat <- inla.stack(
  data = list(y = Y), A = Abar,
  effects = c(
      mesh.index,
      list(Intercept = 1)
    )
)

f.ns <- y ~ -1 + Intercept +  f(field, model = rspde_model)

rspde_nonstat <- inla(f.ns,
  family = "Gamma", 
  data = inla.stack.data(stk.dat),
  control.inla = list(int.strategy = "eb"),
  verbose = TRUE,
  control.predictor = list(A = inla.stack.A(stk.dat), compute = TRUE),
  num.threads = "1:1"
)


f.stat <- y ~ -1 + Intercept +  f(field, model = rspde_stat)

rspde_stat <- inla(f.stat,
  family = "Gamma", 
  data = inla.stack.data(stk.dat),
  control.inla = list(int.strategy = "eb"),
  verbose = FALSE,
  control.predictor = list(A = inla.stack.A(stk.dat), compute = TRUE),
  num.threads = "1:1"
)
