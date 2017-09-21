\name{fractional.operators}
\alias{fractional.operators}
\title{
  Rational approximations of fractional operators
}
\description{
  \code{fractional.operators} is used for computing the operators \eqn{latex}{L_1} and \eqn{latex}{L_2} in an approximation \eqn{latex}{L_1 u(s) = L_2 W} of the fractional SPDE \eqn{latex}{L^\beta (\tau u(s)) = W}. The approximation requires that the equation first i rescaled as \eqn{latex}{(cL)^\beta u(s) = c^{\beta} W} where \eqn{latex}{c} is a constant chosen so that the smallest eigenvalue of \eqn{latex}{cL} is one.
}
\usage{
  fractional.operators(L,
                       beta,
                       C,
                       scale.factor,
                       m=1,
                       tau = 1)
                       }
\arguments{
  \item{L}{A finite element discretization of the operator \eqn{latex}{cL}.}
  \item{beta}{The fractional power.}
  \item{C}{The mass matrix of the finite element discretization.}
  \item{scale.factor}{The constant c}
  \item{m}{The order of the rational approximation, which needs to be a positive integer. The default value is 1.}
  \item{tau}{The constant that scales the variance of the solution. The default value is 1.}
}
\value{
\code{fractional.operators} returns an object of class "rSPDEobj". This is a list that contains the following arguments:
  \item{L1}{The operator \eqn{L_1}{L_1}}
  \item{L2}{The operator \eqn{L_2}{L_2}}
  \item{C}{The mass matrix.}
  \item{m}{The order of the rational approximation.}
  \item{beta}{The fractional power.}
  \item{type}{String indicating the type of approximation.}
  \item{Q}{The matrix \code{t(L1)\%*\%Ci\%*\%L1}.}
  \item{type}{String indicating the type of approximation.}
}
\author{
  David Bolin \email{davidbolin@gmail.com}
}

\seealso{
  \code{\link{matern.operators}}.
}

\examples{
#Compute rational approximation of a Gaussian process with Matern covariance function on R
kappa = 10
sigma = 1
nu = 0.8

#create mass and stiffness matrices for a FEM discretization
x = seq(from = 0, to = 1, length.out = 101)
fem <- rSPDE.fem1D(x)

#compute rational approximation of covariance function at 0.5
tau = sqrt(gamma(nu)/(sigma^2*kappa^(2*nu)*(4*pi)^(1/2)*gamma(nu+1/2)))
op <- fractional.operators(L = fem$G/kappa^2+fem$C, beta=(nu+1/2)/2, C=fem$C, scale.factor=kappa^2, tau = tau)

v = rep(0,101);v[51] = 1
c.approx = op$L2\%*\%solve(op$Q,op$L2\%*\%v)

#plot the result and compare with the true Matern covariance
plot(x,matern.covariance(abs(x-0.5),kappa,nu,sigma),type="l",ylab = "C(h)",xlab="h",
     main = "Matern covariance and rational approximations")
lines(x,c.approx,col=2)
}
