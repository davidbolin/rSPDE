# The Matern covariance function

`matern.covariance` evaluates the Matern covariance function \$\$C(h) =
\frac{\sigma^2}{2^{\nu-1}\Gamma(\nu)}(\kappa h)^\nu K\_\nu(\kappa
h).\$\$

## Usage

``` r
matern.covariance(h, kappa, nu, sigma)
```

## Arguments

- h:

  Distances to evaluate the covariance function at.

- kappa:

  Range parameter.

- nu:

  Shape parameter.

- sigma:

  Standard deviation.

## Value

A vector with the values C(h).

## Examples

``` r
x <- seq(from = 0, to = 1, length.out = 101)
plot(x, matern.covariance(abs(x - 0.5), kappa = 10, nu = 1 / 5, sigma = 1),
  type = "l", ylab = "C(h)", xlab = "h"
)

```
