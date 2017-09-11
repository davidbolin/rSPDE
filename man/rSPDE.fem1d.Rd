\name{rSPDE.fem1d}
\alias{rSPDE.fem1d}
\title{
Finite element approximation for models on R
}
\description{
This function computes mass and stiffness matrices for a FEM approximation on R. These matrices are needed when discretizing the operators in rational approximations.
}
\usage{
rSPDE.fem1d(x)
}
\arguments{
  \item{x}{Locations of the nodes in the FEM approximation
}
}
\value{
The function returns a list with the following elements
\item{G }{The stiffness matrix.}
\item{C }{The mass matrix.}
}
\author{
David Bolin \email{davidbolin@gmail.com}
}

\examples{
#create mass and stiffness matrices for a FEM discretization on [0,1]
x = seq(from = 0, to = 1, length.out = 101)
fem <- rSPDE.fem1D(x)
}
