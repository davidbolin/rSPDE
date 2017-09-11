\name{rSPDE}
\alias{rSPDE}
\title{
  Rational approximations of fractional SPDEs.
}
\description{
  \code{rSPDE} is used for approximating fractional elliptic SPDEs
   \eqn{$L^\beta u(s) = W$}{L^\betau(s) = W}, where \eqn{}
  by rational SPDEs
  \eqn{$L_1u(s) = L_2 W$}{L_1u(s) = L_2 W}
  that can be used for statistical inference.

  The main function for computing the rational operators is \code{fractional.operators}, and a simplified interface for the popular case of Matern covariances is provided through the function \code{matern.operators}. Basic statistical operations such as likelihood evaluations (see \code{rSPDE.loglike}) and kriging predictions (see \code{rSPDE.krig}) using the fractional approximations are also implemented.

  For illustration purposes, the package contains a simple FEM implementation for models on R. For spatial models, the FEM implementation in the \code{R-INLA} package is recommended.
}
\author{
  David Bolin \email{davidbolin@gmail.com}
}
\references{
Inser reference here.
}
