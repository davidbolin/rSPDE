\name{matern.loglike}
\alias{matern.loglike}
\title{
Log-likelihood function for latent Gaussian Matern model
}
\description{
This function evaluates the log-likelihood function for a Gaussian process with a Matern covariance function, that is observed under Gaussian measurement noise: \eqn{Y_i = u(s_i) + \epsilon_i}{Y_i = u(s_i) + \epsilon_i}, where \eqn{\epsilon_i}{\epsilon_i} are iid mean-zero Gaussian variables. The latent model is approximated using a rational approximation of the fractional SPDE model corresponding to the Gaussian process.

}
\usage{
matern.loglike(kappa, sigma, nu, nugget, Y, G, C, A, d = 2, m = 1)
}
\arguments{
  \item{kappa}{Range parameter of the latent process.}
  \item{sigma}{Standard deviation of the latent process.}
  \item{nu}{Shape parameter of the latent process.}
  \item{nugget}{The variance of the measurement noise}
  \item{Y}{The observations, either a vector or a matrix where the columns correspond to independent replicates of observations.}
  \item{G}{The stiffness matrix of a finite element discretization of the domain.}
  \item{C}{The mass matrix of a finite element discretization of the domain.}
  \item{d}{The dimension of the domain. The default value is 2. }
  \item{m}{The order of the rational approximation, which needs to be a positive integer. The default value is 1.}
  \item{A}{A matrix linking the measurement locations to the basis of the FEM approximation of the latent model.}
  \item{d}{The dimension of the domain. The default value is 2.}
  \item{m}{The order of the rational approximation to be used. The default value is 1}
}
\value{
The log-likelihood value.
}
\author{
David Bolin \email{davidbolin@gmail.com}
}

\seealso{
\code{\link{rSPDE.loglike}}, \code{\link{matern.operators}}.
}
\examples{
#this example illustrates how the function can be used for maximum likelihood estimation

#Sample a Gaussian Matern process on R using a rational approximation
sigma = 1
nu = 0.8
practical.range = 0.2
kappa = sqrt(8*nu)/practical.range
nugget = 0.1
n.rep = 20
n.obs = 100
n.x = 101

#create mass and stiffness matrices for a FEM discretization
x = seq(from = 0, to = 1, length.out = n.x)
fem <- rSPDE.fem1D(x)

#compute rational approximation
op <- matern.operators(kappa=kappa,sigma=sigma,nu=nu, G=fem$G,C=fem$C,d=1)

#Sample the model
u <- rSPDE.sample(op,n.rep)

#Create some data
obs.loc <- runif(n=n.obs,min=0,max=1)
A <- rSPDE.A1d(x,obs.loc)
noise <- rnorm(n.obs*n.rep)
dim(noise) <- c(n.obs,n.rep)
Y = as.matrix(A\%*\%u + sqrt(nugget)*noise)

#define negative likelihood function for optimization using matern.loglike
mlik <- function(theta,Y,G,C,A){
return(-matern.loglike(exp(theta[1]),exp(theta[2]),exp(theta[3]),exp(theta[4]),
                         Y=Y, G=G, C=C, A=A, d=1))
}
#Choose some reasonable starting values depending on the size of the domain
theta0 = log(c(sqrt(8)/0.5,sqrt(var(c(Y))),0.9,0.01))

#run estimation and display the results
theta <- optim(theta0,mlik,Y=Y,G=fem$G,C=fem$C,A=A)
print(data.frame(kappa = c(kappa,exp(theta$par[1])),sigma = c(sigma,exp(theta$par[2])),
                    nu = c(nu,exp(theta$par[3])),nugget = c(nugget,exp(theta$par[4])),
                    row.names = c("Truth","Estimates")))
}
