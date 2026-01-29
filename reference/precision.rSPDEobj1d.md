# Get the precision matrix of rSPDEobj1d objects

Function to get the precision matrix of a rSPDEobj1d object

## Usage

``` r
# S3 method for class 'rSPDEobj1d'
precision(
  object,
  loc = NULL,
  nu = NULL,
  kappa = NULL,
  sigma = NULL,
  range = NULL,
  tau = NULL,
  m = NULL,
  ordering = c("field", "location"),
  ldl = FALSE,
  ...
)
```

## Arguments

- object:

  The covariance-based rational SPDE approximation, computed using
  [`matern.rational()`](https://davidbolin.github.io/rSPDE/reference/matern.rational.md)

- loc:

  If non-null, update the locations where to evaluate the model.

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

- ordering:

  Return the matrices ordered by field or by location?

- ldl:

  Directly build the LDL factorization of the precision matrix?

- ...:

  Currently not used.

## Value

A list containing the precision matrix `Q` of the process and its
derivatives if they exist, and a matrix `A` that extracts the elements
corresponding to the process. If `ldl=TRUE`, the LDL factorization is
returned instead of `Q`. If the locations are not ordered, the precision
matrix is given for the ordered locations, but the `A` matrix returns to
the original order.

## See also

[`simulate.rSPDEobj1d()`](https://davidbolin.github.io/rSPDE/reference/simulate.rSPDEobj1d.md),
[`matern.rational()`](https://davidbolin.github.io/rSPDE/reference/matern.rational.md)

## Examples

``` r
# Compute the covariance-based rational approximation of a
# Gaussian process with a Matern covariance function on R
sigma <- 1
nu <- 0.8
range <- 0.2

# create mass and stiffness matrices for a FEM discretization
x <- seq(from = 0, to = 1, length.out = 101)

op_cov <- matern.rational(
  loc = x, nu = nu,
  range = range, sigma = sigma, m = 2,
  parameterization = "matern"
)

# Get the precision matrix:
prec_matrix <- precision(op_cov)
```
