# Precision matrix of the covariance-based rational approximation of stationary Gaussian Matern random fields

`rspde.matern.precision` is used for computing the precision matrix of
the covariance-based rational SPDE approximation of a stationary
Gaussian random fields on \\R^d\\ with a Matern covariance function
\$\$C(h) = \frac{\sigma^2}{2^(\nu-1)\Gamma(\nu)}(\kappa h)^\nu
K\_\nu(\kappa h)\$\$

## Usage

``` r
rspde.matern.precision(
  kappa,
  nu,
  tau = NULL,
  sigma = NULL,
  rspde.order,
  dim,
  fem_mesh_matrices,
  only_fractional = FALSE,
  return_block_list = FALSE,
  type_rational_approx = "brasil"
)
```

## Arguments

- kappa:

  Range parameter of the covariance function.

- nu:

  Shape parameter of the covariance function.

- tau:

  Scale parameter of the covariance function. If sigma is not provided,
  tau should be provided.

- sigma:

  Standard deviation of the covariance function. If tau is not provided,
  sigma should be provided.

- rspde.order:

  The order of the rational approximation

- dim:

  The dimension of the domain

- fem_mesh_matrices:

  A list containing the FEM-related matrices. The list should contain
  elements c0, g1, g2, g3, etc.

- only_fractional:

  Logical. Should only the fractional-order part of the precision matrix
  be returned?

- return_block_list:

  Logical. For `type = "covariance"`, should the block parts of the
  precision matrix be returned separately as a list?

- type_rational_approx:

  Which type of rational approximation should be used? The current types
  are "brasil", "chebfun" or "chebfunLB".

## Value

The precision matrix

## Examples

``` r
set.seed(123)
nobs <- 101
x <- seq(from = 0, to = 1, length.out = nobs)
fem <- rSPDE.fem1d(x)
kappa <- 40
sigma <- 1
d <- 1
nu <- 2.6
tau <- sqrt(gamma(nu) / (kappa^(2 * nu) * (4 * pi)^(d / 2) *
  gamma(nu + d / 2)))
range <- sqrt(8 * nu) / kappa
op_cov <- matern.operators(
  loc_mesh = x, nu = nu, range = range, sigma = sigma,
  d = 1, m = 2, compute_higher_order = TRUE,
  parameterization = "matern"
)
v <- t(rSPDE.A1d(x, 0.5))
c.true <- matern.covariance(abs(x - 0.5), kappa, nu, sigma)
Q <- rspde.matern.precision(
  kappa = kappa, nu = nu, tau = tau, rspde.order = 2, d = 1,
  fem_mesh_matrices = op_cov$fem_mesh_matrices
)
A <- Diagonal(nobs)
Abar <- cbind(A, A, A)
w <- rbind(v, v, v)
c.approx_cov <- (Abar) %*% solve(Q, w)

# plot the result and compare with the true Matern covariance
plot(x, matern.covariance(abs(x - 0.5), kappa, nu, sigma),
  type = "l", ylab = "C(h)",
  xlab = "h", main = "Matern covariance and rational approximations"
)
lines(x, c.approx_cov, col = 2)
```
