# Simulation of a Matern field using a rational SPDE approximation

The function samples a Gaussian random field based on a pre-computed
rational SPDE approximation.

## Usage

``` r
# S3 method for class 'rSPDEobj1d'
simulate(object, nsim = 1, seed = NULL, ...)
```

## Arguments

- object:

  The rational SPDE approximation, computed using
  [`matern.rational()`](https://davidbolin.github.io/rSPDE/reference/matern.rational.md).

- nsim:

  The number of simulations.

- seed:

  an object specifying if and how the random number generator should be
  initialized (‘seeded’).

- ...:

  Currently not used.

## Value

A matrix with the `n` samples as columns.

## See also

[`matern.rational()`](https://davidbolin.github.io/rSPDE/reference/matern.rational.md)

## Examples

``` r
# Sample a Gaussian Matern process on R using a rational approximation
range <- 0.2
sigma <- 1
nu <- 0.8

# compute rational approximation
x <- seq(from = 0, to = 1, length.out = 100)
op <- matern.rational(
  range = range, sigma = sigma,
  nu = nu, loc = x
)

# Sample the model and plot the result
Y <- simulate(op)
plot(x, Y, type = "l", ylab = "u(x)", xlab = "x")

```
