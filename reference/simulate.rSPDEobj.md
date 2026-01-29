# Simulation of a fractional SPDE using a rational SPDE approximation

The function samples a Gaussian random field based on a pre-computed
rational SPDE approximation.

## Usage

``` r
# S3 method for class 'rSPDEobj'
simulate(object, nsim = 1, seed = NULL, ...)
```

## Arguments

- object:

  The rational SPDE approximation, computed using
  [`fractional.operators()`](https://davidbolin.github.io/rSPDE/reference/fractional.operators.md),
  [`matern.operators()`](https://davidbolin.github.io/rSPDE/reference/matern.operators.md),
  or
  [`spde.matern.operators()`](https://davidbolin.github.io/rSPDE/reference/spde.matern.operators.md).

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

[`simulate.CBrSPDEobj()`](https://davidbolin.github.io/rSPDE/reference/simulate.CBrSPDEobj.md)

## Examples

``` r
# Sample a Gaussian Matern process on R using a rational approximation
kappa <- 10
sigma <- 1
nu <- 0.8
range <- sqrt(8 * nu) / kappa

# create mass and stiffness matrices for a FEM discretization
x <- seq(from = 0, to = 1, length.out = 101)
fem <- rSPDE.fem1d(x)

# compute rational approximation
op <- matern.operators(
  range = range, sigma = sigma,
  nu = nu, loc_mesh = x, d = 1,
  parameterization = "matern"
)

# Sample the model and plot the result
Y <- simulate(op)
plot(x, Y, type = "l", ylab = "u(x)", xlab = "x")

```
