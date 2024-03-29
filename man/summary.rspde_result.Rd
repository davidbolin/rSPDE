% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/inla_rspde.R
\name{summary.rspde_result}
\alias{summary.rspde_result}
\title{Summary for posteriors of field parameters for an \code{inla_rspde}
model from a \code{rspde_result} object}
\usage{
\method{summary}{rspde_result}(object, digits = 6, ...)
}
\arguments{
\item{object}{A \code{rspde_result} object.}

\item{digits}{integer, used for number formatting with signif()}

\item{...}{Currently not used.}
}
\value{
Returns a \code{data.frame}
containing the summary.
}
\description{
Summary for posteriors of rSPDE field parameters in
their original scales.
}
\examples{
\donttest{ #devel version
if (requireNamespace("INLA", quietly = TRUE)){
library(INLA)

set.seed(123)

m <- 100
loc_2d_mesh <- matrix(runif(m * 2), m, 2)
mesh_2d <- inla.mesh.2d(
  loc = loc_2d_mesh,
  cutoff = 0.05,
  max.edge = c(0.1, 0.5)
)
sigma <- 1
range <- 0.2
nu <- 0.8
kappa <- sqrt(8 * nu) / range
op <- matern.operators(
  mesh = mesh_2d, nu = nu,
  range = range, sigma = sigma, m = 2,
  parameterization = "matern"
)
u <- simulate(op)
A <- inla.spde.make.A(
  mesh = mesh_2d,
  loc = loc_2d_mesh
)
sigma.e <- 0.1
y <- A \%*\% u + rnorm(m) * sigma.e
Abar <- rspde.make.A(mesh = mesh_2d, loc = loc_2d_mesh)
mesh.index <- rspde.make.index(name = "field", mesh = mesh_2d)
st.dat <- inla.stack(
  data = list(y = as.vector(y)),
  A = Abar,
  effects = mesh.index
)
rspde_model <- rspde.matern(
  mesh = mesh_2d,
  nu.upper.bound = 2
)
f <- y ~ -1 + f(field, model = rspde_model)
rspde_fit <- inla(f,
  data = inla.stack.data(st.dat),
  family = "gaussian",
  control.predictor =
    list(A = inla.stack.A(st.dat))
)
result <- rspde.result(rspde_fit, "field", rspde_model)
summary(result)
}
#devel.tag
}

}
