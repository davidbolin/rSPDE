\name{rSPDE.loglike}
\alias{rSPDE.loglike}
\title{
Log-likelihood function for latent Gaussian fractional SPDE model
}
\description{
This function evaluates the log-likelihood function for a fractional SPDE model \eqn{L^\beta u(s) = W}{L^\beta u(s) = W} that is observed under Gaussian measurement noise: \eqn{Y_i = u(s_i) + \epsilon_i}{Y_i = u(s_i) + \epsilon_i}, where \eqn{\epsilon_i}{\epsilon_i} are iid mean-zero Gaussian variables.
}
\usage{
rSPDE.loglike(obj, Y, A, nugget)
}
\arguments{
  \item{obj}{A rational approximation of the model obtained using \code{fractional.operators} or \code{matern.operators}.}
  \item{Y}{The measurements.}
  \item{A}{An observation matrix that links the measurement location to the finite elemen basis.}
  \item{nugget}{The variance of the measurement noise.}
}
\value{
The log-likelihood value.
}
\author{
David Bolin \email{davidbolin@gmail.com}
}
\note{
This example below shows how the function can be used to evaluate the likelihood of a latent Matern model. Se \code{\link{matern.loglike}} for an example of how this can be used for maximum likelihood estimation.
}
\seealso{
\code{\link{matern.loglike}}, \code{\link{fractional.operators}}, \code{\link{fractional.operators}}.
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

#compute log-likelihood of the data
rSPDE.loglike(op, Y, A, nugget)
}
