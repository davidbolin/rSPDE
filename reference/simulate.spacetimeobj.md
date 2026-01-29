# Simulation of space-time models

Simulation of space-time models

## Usage

``` r
# S3 method for class 'spacetimeobj'
simulate(
  object,
  nsim = 1,
  seed = NULL,
  kappa = NULL,
  sigma = NULL,
  gamma = NULL,
  rho = NULL,
  ...
)
```

## Arguments

- object:

  Space-time object created by
  [`spacetime.operators()`](https://davidbolin.github.io/rSPDE/reference/spacetime.operators.md)

- nsim:

  The number of simulations.

- seed:

  an object specifying if and how the random number generator should be
  initialized (‘seeded’).

- kappa:

  kappa parameter if it should be updated

- sigma:

  sigma parameter if it should be updated

- gamma:

  gamma parameter if it should be updated

- rho:

  rho parameter if it should be updated

- ...:

  Currently not used.

## Value

A matrix with the simulations as columns.

## Examples

``` r
s <- seq(from = 0, to = 20, length.out = 101)
t <- seq(from = 0, to = 20, length.out = 31)

op_cov <- spacetime.operators(space_loc = s, time_loc = t,
                             kappa = 5, sigma = 10, alpha = 1,
                             beta = 2, rho = 1, gamma = 0.05)
x <- simulate(op_cov, nsim = 1) 
image(matrix(x, nrow = length(s), ncol = length(t)))                              
```
