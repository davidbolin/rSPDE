\name{rSPDE.sample}
\alias{rSPDE.sample}
\title{
  Sample fractional SPDE using rational approximation
}
\description{
  \code{rSPDE.sample} is used for sampling the fractional SPDE \eqn{latex}{L^\beta u(s) = W} using a rational approximation.
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
#Compute rational approximation of a Gaussian process with Matern covariance function on R
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
op <- matern.operators(kappa=kappa,sigma=sigma,nu=nu, G=G,C=C,d=1)

Y <- rSPDE.sample(op,1)
plot(x,Y,type="l",ylab = "u(x)",xlab="x")

}
