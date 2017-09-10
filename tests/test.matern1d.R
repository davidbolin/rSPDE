kappa = 10
sigma = 1
nu = 0.8

x = seq(from = 0, to = 1, length.out = 101)

#create mass and stiffness matrices for a FEM discretization
n = length(x)
d <- c(Inf,diff(x))
dm1 = c(d[2:n],Inf)
G = -bandSparse(n=n,m=n,k=c(-1,0,1),diagonals=cBind(1/dm1, -(1/dm1 + 1/d), 1/dm1))
C = bandSparse(n=n,m=n,k=c(-1,0,1),diagonals=cBind(dm1/6, (dm1+d)/3, d/6))
C[1,1:2] <- c(d[2],d[2]/2)/3
C[n,(n-1):n] <- c(d[n]/2,d[n])/3
Ci = Diagonal(n,1/rowSums(C))

#compute rational approximation of covariance function at 0.5
op <- fractional.operators(L = G/kappa^2+C,beta=(nu+1/2)/2,C=C,Ci=Ci,scale.factor=kappa^2)
tau2 = gamma(nu)/(sigma^2*kappa^(2*nu)*(4*pi)^(1/2)*gamma(nu+1/2))
v = rep(0,n);v[51] = 1
c.approx = op$L2%*%solve(tau2*t(op$L1)%*%op$CiL1,op$L2%*%v)

#plot the result and compare with the true Matern covariance
plot(x,matern.covariance(abs(x-0.5),kappa,nu,sigma),type="l",ylab = "C(h)",xlab="h",
     main = "Matern covariance and rational approximation")
lines(x,c.approx,col=2)

op <- matern.operators(kappa=kappa,sigma=sigma,nu=nu, G=G,C=C,d=1)

v = rep(0,n);v[51] = 1
c.approx = op$L2%*%solve(op$Q,op$L2%*%v)

#plot the result and compare with the true Matern covariance
plot(x,matern.covariance(abs(x-0.5),kappa,nu,sigma),type="l",ylab = "C(h)",xlab="h",
     main = "Matern covariance and rational approximation")
lines(x,c.approx,col=2)

Y <- rSPDE.sample(op,1)
plot(x,Y,type="l",ylab = "u(x)",xlab="x")
