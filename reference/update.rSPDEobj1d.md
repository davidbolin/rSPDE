# Update parameters of rSPDEobj1d objects

Function to change the parameters of a rSPDEobj1d object

## Usage

``` r
# S3 method for class 'rSPDEobj1d'
update(
  object,
  nu = NULL,
  alpha = NULL,
  kappa = NULL,
  tau = NULL,
  sigma = NULL,
  range = NULL,
  theta = NULL,
  m = NULL,
  loc = NULL,
  graph = NULL,
  parameterization = NULL,
  type_rational_approximation = object$type_rational_approximation,
  ...
)
```

## Arguments

- object:

  The covariance-based rational SPDE approximation, computed using
  [`matern.rational()`](https://davidbolin.github.io/rSPDE/reference/matern.rational.md)

- nu:

  If non-null, update the shape parameter of the covariance function.
  Will be used if parameterization is 'matern'.

- alpha:

  If non-null, update the fractional SPDE order parameter. Will be used
  if parameterization is 'spde'.

- kappa:

  If non-null, update the parameter kappa of the SPDE. Will be used if
  parameterization is 'spde'.

- tau:

  If non-null, update the parameter tau of the SPDE. Will be used if
  parameterization is 'spde'.

- sigma:

  If non-null, update the standard deviation of the covariance function.
  Will be used if parameterization is 'matern'.

- range:

  If non-null, update the range parameter of the covariance function.
  Will be used if parameterization is 'matern'.

- theta:

  For non-stationary models. If non-null, update the vector of
  parameters.

- m:

  If non-null, update the order of the rational approximation, which
  needs to be a positive integer.

- loc:

  The locations of interest for evaluating the model.

- graph:

  An optional `metric_graph` object.

- parameterization:

  If non-null, update the parameterization.

- type_rational_approximation:

  Which type of rational approximation should be used? The current types
  are "chebfun", "brasil" or "chebfunLB".

- ...:

  Currently not used.

## Value

It returns an object of class "rSPDEobj1d". This object contains the
same quantities listed in the output of
[`matern.rational()`](https://davidbolin.github.io/rSPDE/reference/matern.rational.md).

## See also

[`simulate.rSPDEobj1d()`](https://davidbolin.github.io/rSPDE/reference/simulate.rSPDEobj1d.md),
[`matern.rational()`](https://davidbolin.github.io/rSPDE/reference/matern.rational.md)

## Examples

``` r
s <- seq(from = 0, to = 1, length.out = 101)
kappa <- 20
sigma <- 2
nu <- 0.8
r <- sqrt(8*nu)/kappa #range parameter
op_cov <- matern.rational(loc = s, nu = nu, range = r, sigma = sigma, m = 2, 
parameterization = "matern")
cov1 <- op_cov$covariance(ind = 1)
op_cov <- update(op_cov, range = 0.2)
cov2 <- op_cov$covariance(ind = 1)
plot(s, cov1, type = "l")
lines(s, cov2, col = 2)
```
