\name{fractional.operators}
\alias{fractional.operators}
\title{
  Rational approximations of fractional operators
}
\description{
  \code{fractional.operators} is used for computing the operators \eqn{latex}{L_1} and \eqn{latex}{L_2} in an approximation \eqn{latex}{L_1u(s) = L_2 W} of the fractional SPDE \eqn{latex}{L^\beta u(s) = W}. The approximation requires that the equation first i rescaled as \eqn{latex}{(cL)^\beta u(s) = c^{\beta} W} where \eqn{latex}{c} is a constant chosen so that the smallest eigenvalue of \eqn{latex}{cL} is one.
}
\usage{
  fractional.operators(L,
                       beta,
                       C,
                       Ci,
                       scale.factor,
                       m=1,
                       commutative = TRUE)
                       }
\arguments{
  \item{L}{A finite element discretization of the operator \eqn{latex}{cL}.}
  \item{beta}{The fractional power.}
  \item{C}{The mass matrix of the finite element discretization.}
  \item{Ci}{The inverse of the mass matrix of the finite element discretization.}
  \item{scale.factor}{The constant c}
  \item{m}{The order of the rational approximation, which needs to be a positive integer. The default value is 1.}
  \item{commutative}{TRUE or FALSE depending on if the operators are commutative.}
}
\value{
\code{fractional.operators} returns an object of class "rSPDEobj". This is a list that contains the following arguments:
  \item{L1}{The operator \eqn{L_1}{L_1}}
  \item{L2}{The operator \eqn{L_2}{L_2}}
  \item{C}{The mass matrix.}
  \item{CiL1}{The operator \eqn{L_1}{L_1} times the inverse of the mass matrix.}
  \item{m}{The order of the rational approximation.}
  \item{beta}{The fractional power.}
  \item{type}{String indicating the type of approximation.}
  \item{Q}{Either the matrix \code{t(L1)\%*\%CiL1} or \code{L2\%*\%C\%*\%t(L2)} depending on if the operators are commutative or not}
  \item{type}{String indicating the type of approximation.}
  \item{commutative}{TRUE or FALSE depending on if the operators are commutative.}
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

}
