% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fractional.operators.R
\name{spde.matern.operators}
\alias{spde.matern.operators}
\title{Rational approximations of non-stationary Gaussian SPDE Matern random fields}
\usage{
spde.matern.operators(kappa, tau, nu, G, C, d, m = 1)
}
\arguments{
\item{kappa}{Vector with the, possibly spatially varying, range parameter
evaluated at the locations of the mesh used for the finite element
discretization of the SPDE.}

\item{tau}{Vector with the, possibly spatially varying, precision
parameter evaluated at the locations
of the mesh used for the finite element discretization of the SPDE.}

\item{nu}{Shape parameter of the covariance function, related to
\eqn{\beta} through the equation
\eqn{\beta = (\nu + d/2)/2}.}

\item{G}{The stiffness matrix of a finite element discretization of
the domain of interest.}

\item{C}{The mass matrix of a finite element discretization of the
domain of interest.}

\item{d}{The dimension of the domain.}

\item{m}{The order of the rational approximation, which needs to
be a positive integer.
The default value is 1.}
}
\value{
\code{spde.matern.operators} returns an object of
class "rSPDEobj. This object contains the
quantities listed in the output of \code{\link[=fractional.operators]{fractional.operators()}}
as well as the smoothness parameter \eqn{\nu}.
}
\description{
\code{spde.matern.operators} is used for computing a rational SPDE
approximation of a Gaussian random
fields on \eqn{R^d} defined as a solution to the SPDE
\deqn{(\kappa(s) - \Delta)^\beta (\tau(s)u(s)) = W.}
}
\details{
The approximation is based on a rational approximation of the
fractional operator \eqn{(\kappa(s)^2 -\Delta)^\beta}, where
\eqn{\beta = (\nu + d/2)/2}. This results in an approximate model
on the form \deqn{P_l u(s) = P_r W,} where \eqn{P_j = p_j(L)} are
non-fractional operators defined in terms of polynomials \eqn{p_j} for
\eqn{j=l,r}. The order of \eqn{p_r} is given by \code{m} and the order
of \eqn{p_l} is \eqn{m + m_\beta} where \eqn{m_\beta} is the integer
part of \eqn{\beta} if \eqn{\beta>1} and \eqn{m_\beta = 1} otherwise.

The discrete approximation can be written as \eqn{u = P_r x} where
\eqn{x \sim N(0,Q^{-1})}{x ~ N(0,Q^{-1})}
and \eqn{Q = P_l^T C^{-1} P_l}. Note that the matrices \eqn{P_r} and
\eqn{Q} may be be ill-conditioned for \eqn{m>1}.
In this case, the metehods in \code{\link[=operator.operations]{operator.operations()}}
should be used for operations involving the matrices, since
these methods are more numerically stable.
}
\examples{
# Sample non-stationary Matern field on R
tau <- 1
nu <- 0.8

# create mass and stiffness matrices for a FEM discretization
x <- seq(from = 0, to = 1, length.out = 101)
fem <- rSPDE.fem1d(x)

# define a non-stationary range parameter
kappa <- seq(from = 2, to = 20, length.out = length(x))

# compute rational approximation
op <- spde.matern.operators(
  kappa = kappa, tau = tau, nu = nu,
  G = fem$G, C = fem$C, d = 1
)

# sample the field
u <- simulate(op)

# plot the sample
plot(x, u, type = "l", ylab = "u(s)", xlab = "s")

}
\seealso{
\code{\link[=fractional.operators]{fractional.operators()}},
\code{\link[=spde.matern.operators]{spde.matern.operators()}},
\code{\link[=matern.operators]{matern.operators()}}
}