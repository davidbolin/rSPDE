# Simulation of a fractional SPDE using the covariance-based rational SPDE approximation

The function samples a Gaussian random field based using the
covariance-based rational SPDE approximation.

## Usage

``` r
# S3 method for class 'CBrSPDEobj2d'
simulate(
  object,
  nsim = 1,
  seed = NULL,
  nu = NULL,
  hx = NULL,
  hy = NULL,
  hxy = NULL,
  sigma = NULL,
  m = NULL,
  ...
)
```

## Arguments

- object:

  The covariance-based rational SPDE approximation, computed using
  [`matern2d.operators()`](https://davidbolin.github.io/rSPDE/reference/matern2d.operators.md)

- nsim:

  The number of simulations.

- seed:

  An object specifying if and how the random number generator should be
  initialized (‘seeded’).

- nu:

  If non-null, update the shape parameter of the covariance function.

- hx:

  If non-null, update the hx parameter.

- hy:

  If non-null, update the hy parameter.

- hxy:

  If non-null, update the hxy parameter.

- sigma:

  If non-null, update the standard deviation of the covariance function.

- m:

  If non-null, update the order of the rational approximation, which
  needs to be a positive integer.

- ...:

  Currently not used.

## Value

A matrix with the `n` samples as columns.

## Examples

``` r
library(fmesher)
n_loc <- 2000
loc_2d_mesh <- matrix(runif(n_loc * 2), n_loc, 2)
mesh_2d <- fm_mesh_2d(loc = loc_2d_mesh, cutoff = 0.03, max.edge = c(0.1, 0.5))
op <- matern2d.operators(mesh = mesh_2d, sigma = 1, nu = 1, hx = 0.1, hy = 0.1, hxy = 0)
u <- simulate(op)
```
