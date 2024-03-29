% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fractional.computations.R
\name{predict.CBrSPDEobj}
\alias{predict.CBrSPDEobj}
\title{Prediction of a fractional SPDE using the covariance-based
rational SPDE approximation}
\usage{
\method{predict}{CBrSPDEobj}(
  object,
  A,
  Aprd,
  Y,
  sigma.e,
  mu = 0,
  compute.variances = FALSE,
  posterior_samples = FALSE,
  n_samples = 100,
  only_latent = FALSE,
  ...
)
}
\arguments{
\item{object}{The covariance-based rational SPDE approximation,
computed using \code{\link[=matern.operators]{matern.operators()}}}

\item{A}{A matrix linking the measurement locations to the basis of the FEM
approximation of the latent model.}

\item{Aprd}{A matrix linking the prediction locations to the basis of the
FEM approximation of the latent model.}

\item{Y}{A vector with the observed data, can also be a matrix where the
columns are observations
of independent replicates of \eqn{u}.}

\item{sigma.e}{The standard deviation of the Gaussian measurement noise.
Put to zero if the model does not have measurement noise.}

\item{mu}{Expectation vector of the latent field (default = 0).}

\item{compute.variances}{Set to also TRUE to compute the kriging variances.}

\item{posterior_samples}{If \code{TRUE}, posterior samples will be returned.}

\item{n_samples}{Number of samples to be returned. Will only be used if \code{sampling} is \code{TRUE}.}

\item{only_latent}{Should the posterior samples be only given to the laten model?}

\item{...}{further arguments passed to or from other methods.}
}
\value{
A list with elements
\item{mean }{The kriging predictor (the posterior mean of u|Y).}
\item{variance }{The posterior variances (if computed).}
}
\description{
The function is used for computing kriging predictions based
on data \eqn{Y_i = u(s_i) + \epsilon_i}, where \eqn{\epsilon}{\epsilon}
is mean-zero Gaussian measurement noise and \eqn{u(s)}{u(s)} is defined by
a fractional SPDE \eqn{(\kappa^2 I - \Delta)^{\alpha/2} (\tau u(s)) = W},
where \eqn{W}{W} is Gaussian white noise and \eqn{\alpha = \nu + d/2},
where \eqn{d} is the dimension of the domain.
}
\examples{
set.seed(123)
# Sample a Gaussian Matern process on R using a rational approximation
kappa <- 10
sigma <- 1
nu <- 0.8
sigma.e <- 0.3
range <- 0.2

# create mass and stiffness matrices for a FEM discretization
x <- seq(from = 0, to = 1, length.out = 101)
fem <- rSPDE.fem1d(x)

tau <- sqrt(gamma(nu) / (sigma^2 * kappa^(2 * nu) *
       (4 * pi)^(1 / 2) * gamma(nu + 1 / 2)))

# Compute the covariance-based rational approximation
op_cov <- matern.operators(
  loc_mesh = x, nu = nu,
  range = range, sigma = sigma, d = 1, m = 2,
  parameterization = "matern"
)

# Sample the model
u <- simulate(op_cov)

# Create some data
obs.loc <- runif(n = 10, min = 0, max = 1)
A <- rSPDE.A1d(x, obs.loc)
Y <- as.vector(A \%*\% u + sigma.e * rnorm(10))

# compute kriging predictions at the FEM grid
A.krig <- rSPDE.A1d(x, x)
u.krig <- predict(op_cov,
  A = A, Aprd = A.krig, Y = Y, sigma.e = sigma.e,
  compute.variances = TRUE
)

plot(obs.loc, Y,
  ylab = "u(x)", xlab = "x", main = "Data and prediction",
  ylim = c(
    min(u.krig$mean - 2 * sqrt(u.krig$variance)),
    max(u.krig$mean + 2 * sqrt(u.krig$variance))
  )
)
lines(x, u.krig$mean)
lines(x, u.krig$mean + 2 * sqrt(u.krig$variance), col = 2)
lines(x, u.krig$mean - 2 * sqrt(u.krig$variance), col = 2)
}
