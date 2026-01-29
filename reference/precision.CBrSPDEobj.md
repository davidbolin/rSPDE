# Get the precision matrix of CBrSPDEobj objects

Function to get the precision matrix of a CBrSPDEobj object

## Usage

``` r
precision(object, ...)

# S3 method for class 'CBrSPDEobj'
precision(
  object,
  nu = NULL,
  kappa = NULL,
  sigma = NULL,
  range = NULL,
  tau = NULL,
  m = NULL,
  ...
)
```

## Arguments

- object:

  The covariance-based rational SPDE approximation, computed using
  [`matern.operators()`](https://davidbolin.github.io/rSPDE/reference/matern.operators.md)

- ...:

  Currently not used.

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

The precision matrix.

## See also

[`simulate.CBrSPDEobj()`](https://davidbolin.github.io/rSPDE/reference/simulate.CBrSPDEobj.md),
[`matern.operators()`](https://davidbolin.github.io/rSPDE/reference/matern.operators.md)

## Examples

``` r
# Compute the covariance-based rational approximation of a
# Gaussian process with a Matern covariance function on R
kappa <- 10
sigma <- 1
nu <- 0.8
range <- 0.2

# create mass and stiffness matrices for a FEM discretization
x <- seq(from = 0, to = 1, length.out = 101)
fem <- rSPDE.fem1d(x)

# compute rational approximation of covariance function at 0.5
tau <- sqrt(gamma(nu) / (sigma^2 * kappa^(2 * nu) *
  (4 * pi)^(1 / 2) * gamma(nu + 1 / 2)))
op_cov <- matern.operators(
  loc_mesh = x, nu = nu,
  range = range, sigma = sigma, d = 1, m = 2,
  parameterization = "matern"
)

# Get the precision matrix:
prec_matrix <- precision(op_cov)
```
