# Rational approximation of the Matern covariance

Computes a rational approximation of the Matern covariance function on
intervals.

## Usage

``` r
matern.rational.cov(
  h,
  order,
  kappa,
  nu,
  sigma,
  type_rational = "brasil",
  type_interp = "linear"
)
```

## Arguments

- h:

  Distances to compute the covariance for

- order:

  The order of the approximation

- kappa:

  Range parameter

- nu:

  Smoothness parameter

- sigma:

  Standard deviation

- type_rational:

  Method used to compute the coefficients of the rational approximation.

- type_interp:

  Interpolation method for the rational coefficients.

## Value

The covariance matrix of the approximation

## Examples

``` r
h <- seq(from = 0, to = 1, length.out = 100)
cov.true <- matern.covariance(h, kappa = 10, sigma = 1, nu = 0.8)
cov.approx <- matern.rational.cov(h, kappa = 10, sigma = 1, nu = 0.8, order = 2)

plot(h, cov.true)
lines(h, cov.approx, col = 2)

 
```
