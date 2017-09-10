\name{matern.operators}
\alias{matern.operators}
\title{
  Rational approximations of Matern covariances
}
\description{
  \code{matern.operators} is used for an approximation \eqn{u(s) = L_2 x(s)}{u(s) = L_2 x(s)},  \eqn{x(s) ~ N(0,Q^{-1})}{x(s) ~ N(0,Q^-1)} of the fractional SPDE corresponding to a Gaussian field on \eqn{latex}{R^d} with a Matern covariance function
  \deqn{C(h) = \frac{\sigma^2}{2^(\nu-1)\Gamma(\nu)}(\kappa h)^\nu K_\nu(\kappa h)}{C(h) = (\sigma^2/(2^(\nu-1)\Gamma(\nu))(\kappa h)^\nu K_\nu(\kappa h)}
}
\usage{
  matern.operators(kappa,
                   sigma,
                   nu,
                   G,
                   C,
                   d=2,
                   m=1)
                       }
\arguments{
  \item{kappa}{Range parameter.}
  \item{sigma}{Standard deviation.}
  \item{nu}{Shape parameter.}
  \item{G}{The stiffness matrix of a finite element discretization of the domain.}
  \item{C}{The mass matrix of a finite element discretization of the domain.}
  \item{d}{The dimension of the domain. The default value is 2. }
  \item{m}{The order of the rational approximation, which needs to be a positive integer. The default value is 1.}
}
\value{
\code{fractional.operators} returns an object of class "rSPDEobj". This is a list that contains the following arguments:
  \item{L2}{The sparse matrix \eqn{L_2}{L_2}}
  \item{Q}{The sparse precision matrix \eqn{Q}{Q}}
  \item{kappa}{Range parameter.}
  \item{sigma}{Standard deviation.}
  \item{nu}{Shape parameter.}
  \item{m}{The order of the rational approximation, which needs to be a positive integer. The default value is 1.}
  \item{type}{String indicating the type of approximation.}
}
\author{
  David Bolin \email{davidbolin@gmail.com}
}

\references{
Inser reference here.
}

\seealso{
  \code{\link{fractional.operators}}.
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

v = rep(0,n);v[51] = 1
c.approx = op$L2%*%solve(op$Q,op$L2%*%v)

#plot the result and compare with the true Matern covariance
plot(x,matern.covariance(abs(x-0.5),kappa,nu,sigma),type="l",ylab = "C(h)",xlab="h",
     main = "Matern covariance and rational approximation")
lines(x,c.approx,col=2)

}
