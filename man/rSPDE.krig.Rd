\name{rSPDE.krig}
\alias{rSPDE.krig}
\title{
Compute kriging predictions for a Gaussian fractional SPDE using rational approximations
}
\description{
\code{rSPDE.krig} is used for computing kriging predictions based on data \eqn{Y_i = u(s_i) + \epsilon_i}, where \eqn{\epsilon}{\epsilon} is mean-zero Gaussian measurement noise and \eqn{u(s)}{u(s)} is defined by a fractional SPDE \eqn{L^\beta u(s) = W}{L^\beta u(s) = W}, where \eqn{W}{W} is Gaussian white noise.
}
\usage{
rSPDE.krig(obj, A, Aprd, Y, nugget)
}
\arguments{
\item{obj}{A rational approximation of the model obtained using \code{fractional.operators} or \code{matern.operators}.}
\item{A}{A matrix linking the measurement locations to the basis of the FEM approximation of the latent model.
}
  \item{Aprd}{
A matrix linking the prediction locations to the basis of the FEM approximation of the latent model.
}
  \item{Y}{A vector with the observed data.
}
  \item{nugget}{The variance of the Gaussian measurement noise.
}
}
\value{A vector with the kriging prediction.
}
\author{
David Bolin \email{davidbolin@gmail.com}
}

\seealso{
\code{\link{fractional.operators}}, \code{\link{matern.operators}}.
}
\examples{
#Sample a Gaussian Matern process on R using a rational approximation
kappa = 10
sigma = 1
nu = 0.8
nugget = 0.1

#create mass and stiffness matrices for a FEM discretization
x = seq(from = 0, to = 1, length.out = 101)
fem <- rSPDE.fem1D(x)

#compute rational approximation
op <- matern.operators(kappa=kappa,sigma=sigma,nu=nu, G=fem$G,C=fem$C,d=1)

#Sample the model
u <- rSPDE.sample(op,1)

#Create some data
obs.loc <- runif(n=10,min=0,max=1)
A <- rSPDE.A1d(x,obs.loc)
Y = as.vector(A\%*\%u + sqrt(nugget)*rnorm(10))

#compute kriging predictions at the FEM grid
A.krig <- rSPDE.A1d(x,x)
u.krig <- rSPDE.krig(op,A=A,Aprd = A.krig,Y=Y,nugget = nugget)

plot(obs.loc,Y,ylab = "u(x)",xlab="x",main="Data and kriging")
lines(x,u.krig,ylab = "u(x)",xlab="x")
}
