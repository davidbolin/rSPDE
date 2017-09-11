\name{rSPDE.A1d}
\alias{rSPDE.A1d}
\title{
Observation matrix for finite element discretization on R
}
\description{
A finite element discretization on R can be written as \eqn{u(s) = \sum_i^n u_i \varphi_i(s)}{u(s) = \sum_i^n u_i \varphi_i(s)} where \eqn{\varphi_i(s)} is a piecewise linear "hat function" centered at location \eqn{x_i}{x_i}. This function computes an \eqn{m\times n}{m x n} matrix \eqn{A}{A} that links the basis function in the expansion to specified locations \eqn{s = (s_1,\ldots, s_m)} in the domain through \eqn{A_ij = \varphi_j(s_i)}{A_ij = \varphi_j(s_i)}.
}
\usage{
rSPDE.A1d(x, loc)
}
\arguments{
  \item{x}{The locations of the nodes in the FEM discretization.}
  \item{loc}{The locations \eqn{(s_1,\ldots, s_m)}}
}
\value{The sparse matrix A.
}
\author{
David Bolin \email{davidbolin@gmail.com}
}

\seealso{
\code{\link{rSPDE.fem1d}}}

\examples{
#create mass and stiffness matrices for a FEM discretization on [0,1]
x = seq(from = 0, to = 1, length.out = 101)
fem <- rSPDE.fem1D(x)

#create the observation matrix for some locations in the domain
obs.loc <- runif(n=10,min=0,max=1)
A <- rSPDE.A1d(x,obs.loc)
}
