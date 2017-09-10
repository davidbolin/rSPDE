library(INLA)
library(rSPDE)
library(fields)

nu = 0.5
r = 0.2
kappa = sqrt(8*nu)/r
sigma = 4

x <- seq(from = 0, to = 1, length.out = 80)
mesh <- inla.mesh.create(lattice = inla.mesh.lattice(x = x, y = x))
fem = inla.fmesher.smorg(mesh$loc, mesh$graph$tv, fem = 2, output = list("c0", "c1", "g1"))
C <- fem$c1
G <- fem$g1

operators <- matern.operators(kappa=kappa,sigma=sigma,nu=nu,G=G,C=C,m=2)
Q <- operators$Q
L2 <- operators$L2

A <- inla.spde.make.A(mesh = mesh, loc = cbind(0.5,0.5))
c <- L2%*%solve(Q,t(L2)%*%t(A))

proj <- inla.mesh.projector(mesh, dims = c(100,100))
d = sqrt(rowSums((mesh$loc[,1:2]-c(0.5,0.5))^2))
c.true <- matern.covariance(d,kappa,nu,sigma)

par(mfrow=c(2,2))
image.plot(proj$x, proj$y, inla.mesh.project(proj, field = as.vector(c)),main="rational approx")
image.plot(proj$x, proj$y, inla.mesh.project(proj, field = as.vector(c.true)),main="true cov")
image.plot(proj$x, proj$y, inla.mesh.project(proj, field = as.vector(c-c.true)),main="diff")
