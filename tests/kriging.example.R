#This example replicates the results of the Application in the article
# "Extending the SPDE approach to Gaussian random fields with general smoothness"
#rm(list=ls())
compute.optimal = TRUE
compute.tapering = TRUE
##################################
#load needed packages
##################################
library(INLA)
library(maps)
library(spam)
library(splancs)
library(excursions)
library(fields)
library(rSPDE)

##################################
#load data
##################################
data("USprecip")
ind = USprecip[,5]==1;
lon = USprecip[ind,1]
lat = USprecip[ind,2]
Y = USprecip[ind,4]

loc = cbind(lon,lat)


##################################
#compute prediction locations
##################################
lattice = inla.mesh.lattice(x = seq(from = min(lon), to = max(lon), by = 0.25),
                            y = seq(from = min(lat), to = max(lat), by = 0.25))

usa <- map('usa')
border <- as.matrix(cbind(usa$x,usa$y))
border.points <- border[!is.na(border[,1]),]

xy.in <- inout(lattice$loc, border.points)
loc.prd = lattice$loc[xy.in,]
lon.prd = lattice$loc[xy.in,1]
lat.prd = lattice$loc[xy.in,2]

##################################
#Parameters of covariance function
##################################
kappa1 = 1/40.73
kappa2 = 1/523.73
phi1 = sqrt(0.277)
phi2 = sqrt(0.722)

##################################
#setup for optimal prediction
##################################
if(compute.optimal){
  t <- proc.time()
  D2 = rdist.earth(loc)
  Sigmao = matern.covariance(D2,kappa1,phi1,nu=1/2) + matern.covariance(D2,kappa2,phi2,nu=1/2)
  D1 = rdist.earth(loc.prd,loc)
  Sigmap = matern.covariance(D1,kappa1,phi1,nu=1/2) + matern.covariance(D1,kappa2,phi2,nu=1/2)
  t.opt1 <- proc.time()-t
}

##################################
#compute optimal prediction
##################################
if(compute.optimal){
  t <- proc.time()
  Yhat = Sigmap%*%solve(Sigma,Y)
  t.opt2 <- proc.time()-t
}

##################################
#setup for tapering prediction
##################################
if(compute.tapering){
  t <- proc.time()
  exp.mix.cov <- function(distance, ...)
    phi1^2 * exp( -abs(distance)*kappa1) + phi2^2 * exp( -abs(distance)*kappa2)

  tap.cov1 <- stationary.taper.cov(loc,loc, Covariance="exp.mix.cov",
                                   Taper="Wendland", Dist.args = list(method="greatcircle"),Taper.args = list(theta=50.0,k=1))

  tap.cov2 <- stationary.taper.cov(loc.prd,loc, Covariance="exp.mix.cov",
                                   Taper="Wendland", Dist.args=list("greatcircle"),Taper.args = list(theta=50.0,k=1))
  t.tap1 <- proc.time() -t
}

##################################
#compute tapering prediction
##################################
if(compute.tapering){
  t <- proc.time()
  Ytap <- tap.cov2%*%solve(tap.cov1,Y)
  t.tap2 <- proc.time()-t
  sqrt(sum((Yhat-Ytap)^2))
}

##################################
#setup for rational approximation
##################################
nugget = 1e-3
R = 3963.34
loc.xy = cBind(cos(lat*pi/180)*cos(lon*pi/180),
               cos(lat*pi/180)*sin(lon*pi/180),
               sin(lat*pi/180))
loc.prd.xy = cBind(cos(lat.prd*pi/180)*cos(lon.prd*pi/180),
                   cos(lat.prd*pi/180)*sin(lon.prd*pi/180),
                   sin(lat.prd*pi/180))

t <- proc.time()
#compute lon-lat mesh and convert to spherical coordinates
m1 <- inla.mesh.2d(loc=loc,max.edge=c(4,8), cutoff=0.55)
Plon <- m1$loc
P <- cBind(cos(m1$loc[,2]*pi/180)*cos(m1$loc[,1]*pi/180),
      cos(m1$loc[,2]*pi/180)*sin(m1$loc[,1]*pi/180),
      sin(m1$loc[,2]*pi/180))
m1$loc = P
m1$manifold = "sphere"
#compute FEM matrices
fem = inla.fmesher.smorg(m1$loc, m1$graph$tv, fem = 2, output = list("c0", "c1", "g1"))
Aprd <- inla.spde.make.A(m1,loc=loc.prd.xy)
A <- inla.spde.make.A(m1, loc=loc.xy)

#compute rational approximations
op1 <- matern.operators(kappa=kappa1*R,sigma=phi1,nu=1/2,G=fem$g1,C=fem$c1)
op2 <- matern.operators(kappa=kappa2*R,sigma=phi2,nu=1/2,G=fem$g1,C=fem$c1)

#compute posterior precision matrix
A = cBind(A%*%op1$L2, A%*%op2$L2)
Qhat = bdiag(op1$Q,op2$Q) + t(A)%*%A/nugget
t.rational1 <- proc.time() - t

##################################
#compute rational prediction
##################################
t <- proc.time()
tmp = inla.qsolve(Qhat,t(A)%*%Y/nugget,reordering="amd")
#tmp = solve(Qhat,t(A)%*%Y/nugget)
Yhat_r = cBind(Aprd%*%op1$L2, Aprd%*%op2$L2)%*%tmp
t.rational2 <- proc.time() - t

if(0){
  kappa1 = 1/40
  kappa2 = 1/500
  op1 <- matern.operators(kappa=kappa1*R,sigma=1,nu=1/2,G=fem$g1,C=fem$c1)
  op2 <- matern.operators(kappa=kappa2*R,sigma=1,nu=1/2,G=fem$g1,C=fem$c1)
  v = rep(0,dim(op1$Q)[1])
  v[500] = 1
  c = op1$L2%*%solve(op1$Q,t(op1$L2)%*%v)+op2$L2%*%solve(op2$Q,t(op2$L2)%*%v)

  d = rdist.earth(Plon[,1:2],t(Plon[500,1:2]))
  c.true <- exponential.covariance(d,kappa1,1) + exponential.covariance(d,kappa2,1)
  par(mfrow=c(1,2))
  quilt.plot(Plon[,1],Plon[,2],as.vector(c))
  quilt.plot(Plon[,1],Plon[,2],c.true)
}
par(mfrow=c(1,2))
quilt.plot(loc.prd[,1],loc.prd[,2],as.vector(Yhat_r))
quilt.plot(loc[,1],loc[,2],as.vector(Y))


##################################
#display results
##################################
if(1){
  errors = c(0,sqrt(sum((Yhat-Ytap)^2)),sqrt(sum((Yhat-Yhat_r)^2)))
  time1 = c(t.opt1[1]+t.opt1[4],t.tap1[1]+t.tap1[4],t.rational1[1]+t.rational1[4])
  time2 = c(t.opt2[1]+t.opt2[4],t.tap2[1]+t.tap2[4],t.rational2[1]+t.rational2[4])
  results <- data.frame(error=errors,setup.time = time1,kriging.time = time2,total.time = time1+time2,
                        row.names = c("Optimal","Tapering","Rational"))
  print(t(results))
}
