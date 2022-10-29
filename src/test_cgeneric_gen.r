library(INLA)
library(rSPDE)


data(PRprec)
data(PRborder)

Y <- rowMeans(PRprec[, 3 + 1:31])
ind <- !is.na(Y)
Y <- Y[ind]
coords <- as.matrix(PRprec[ind, 1:2])
alt <- PRprec$Altitude[ind]
seaDist <- apply(spDists(coords, PRborder[1034:1078, ],
  longlat = TRUE
), 1, min)


prdomain <- inla.nonconvex.hull(coords, -0.03, -0.05, resolution = c(100, 100))
prmesh <- inla.mesh.2d(boundary = prdomain, max.edge = c(0.45, 1), cutoff = 0.2)

Abar <- rspde.make.A(mesh = prmesh, loc = coords)

mesh.index <- rspde.make.index(name = "field", mesh = prmesh,)

rspde_model <- rspde.matern(mesh = prmesh)

stk.dat <- inla.stack(
  data = list(y = Y), A = list(Abar, 1), tag = "est",
  effects = list(
    c(
      mesh.index,
      list(Intercept = 1)
    ),
    list(
      seaDist = inla.group(seaDist)
    )
  )
)

f.s <- y ~ -1 + Intercept + f(seaDist, model="rw1")+
  f(field, model = rspde_model)

rspde_fit <- inla(f.s,
  family = "Gamma", data = inla.stack.data(stk.dat),
  verbose = TRUE,
  control.inla = list(int.strategy = "eb"),
  control.predictor = list(A = inla.stack.A(stk.dat), compute = TRUE),
            inla.mode = "experimental"
)

