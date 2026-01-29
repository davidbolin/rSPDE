# Object-based log-likelihood function for latent Gaussian fractional SPDE model using the rational approximations

This function evaluates the log-likelihood function for a Gaussian
process with a Matern covariance function, that is observed under
Gaussian measurement noise: \\Y_i = u(s_i) + \epsilon_i\\, where
\\\epsilon_i\\ are iid mean-zero Gaussian variables. The latent model is
approximated using the a rational approximation of the fractional SPDE
model corresponding to the Gaussian process.

## Usage

``` r
rSPDE.matern.loglike(
  object,
  Y,
  A,
  sigma.e,
  mu = 0,
  nu = NULL,
  kappa = NULL,
  sigma = NULL,
  range = NULL,
  tau = NULL,
  m = NULL
)
```

## Arguments

- object:

  The rational SPDE approximation, computed using
  [`matern.operators()`](https://davidbolin.github.io/rSPDE/reference/matern.operators.md)

- Y:

  The observations, either a vector or a matrix where the columns
  correspond to independent replicates of observations.

- A:

  An observation matrix that links the measurement location to the
  finite element basis.

- sigma.e:

  The standard deviation of the measurement noise.

- mu:

  Expectation vector of the latent field (default = 0).

- nu:

  If non-null, update the shape parameter of the covariance function.

- kappa:

  If non-null, update the range parameter of the covariance function.

- sigma:

  If non-null, update the standard deviation of the covariance function.

- range:

  If non-null, update the range parameter of the covariance function.

- tau:

  If non-null, update the parameter tau.

- m:

  If non-null, update the order of the rational approximation, which
  needs to be a positive integer.

## Value

The log-likelihood value.

## See also

[`matern.operators()`](https://davidbolin.github.io/rSPDE/reference/matern.operators.md),
[`predict.CBrSPDEobj()`](https://davidbolin.github.io/rSPDE/reference/predict.CBrSPDEobj.md)

## Examples

``` r
# this example illustrates how the function can be used for maximum likelihood estimation

set.seed(123)
# Sample a Gaussian Matern process on R using a rational approximation
nu <- 0.8
kappa <- 5
sigma <- 1
sigma.e <- 0.1
n.rep <- 10
n.obs <- 100
n.x <- 51
range <- 0.2

# create mass and stiffness matrices for a FEM discretization
x <- seq(from = 0, to = 1, length.out = n.x)
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
u <- simulate(op_cov, n.rep)

# Create some data
obs.loc <- runif(n = n.obs, min = 0, max = 1)
A <- rSPDE.A1d(x, obs.loc)
noise <- rnorm(n.obs * n.rep)
dim(noise) <- c(n.obs, n.rep)
Y <- as.matrix(A %*% u + sigma.e * noise)

# Define the negative likelihood function for optimization
# using CBrSPDE.matern.loglike

# Notice that we are also using sigma instead of tau, so it can be compared
# to matern.loglike()
mlik_cov <- function(theta, Y, A, op_cov) {
  kappa <- exp(theta[1])
  sigma <- exp(theta[2])
  nu <- exp(theta[3])
  return(-rSPDE.matern.loglike(
    object = op_cov, Y = Y,
    A = A, kappa = kappa, sigma = sigma,
    nu = nu, sigma.e = exp(theta[4])
  ))
}

# The parameters can now be estimated by minimizing mlik with optim
# \donttest{
# Choose some reasonable starting values depending on the size of the domain
theta0 <- log(c(sqrt(8), 1 / sqrt(var(c(Y))), 0.9, 0.01))

# run estimation and display the results
theta <- optim(theta0, mlik_cov,
  Y = Y, A = A, op_cov = op_cov,
  method = "L-BFGS-B"
)

print(data.frame(
  range = c(range, exp(theta$par[1])), sigma = c(sigma, exp(theta$par[2])),
  nu = c(nu, exp(theta$par[3])), sigma.e = c(sigma.e, exp(theta$par[4])),
  row.names = c("Truth", "Estimates")
))
#>              range    sigma        nu    sigma.e
#> Truth     0.200000 1.000000 0.8000000 0.10000000
#> Estimates 2.828427 1.025364 0.8792568 0.09830274
# }
```
