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

rspde_model <- rspde.matern2(mesh = prmesh)

stk.dat <- inla.stack(
  data = list(y = Y), A = Abar,
  effects = c(
      mesh.index,
      list(Intercept = 1)
    )
)

f.s <- y ~ -1 + Intercept +  f(field, model = rspde_model)

rspde_fit <- inla(f.s,
  family = "Gamma", 
  data = inla.stack.data(stk.dat),
  control.inla = list(int.strategy = "eb"),
  verbose = TRUE,
  control.predictor = list(A = inla.stack.A(stk.dat), compute = TRUE),
  num.threads = "1:1"
)
