\name{rSPDE.sample}
\alias{rSPDE.sample}
\title{
  Sample Gaussian fractional SPDE using rational approximations
}
\description{
  \code{rSPDE.sample} is used for sampling the fractional SPDE \eqn{L^\beta u(s) = W}{L^\beta u(s) = W}, where \eqn{W}{W} is Gaussian white noise, using a rational approximation.
}
\usage{
  rSPDE.sample(obj, n)
}
\arguments{
  \item{obj}{A rational approximation of the model obtained using \code{fractional.operators} or \code{matern.operators}.}
  \item{n}{The number of samples required.}
}
\value{
A matrix with the samples as columns.
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

#create mass and stiffness matrices for a FEM discretization
x = seq(from = 0, to = 1, length.out = 101)
fem <- rSPDE.fem1D(x)

#compute rational approximation
op <- matern.operators(kappa=kappa,sigma=sigma,nu=nu, G=fem$G,C=fem$C,d=1)

#Sample the model and plot the result
set.seed(1)
Y <- rSPDE.sample(op,1)
plot(x,Y,type="l",ylab = "u(x)",xlab="x")
}
