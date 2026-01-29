# Simulation of a fractional intrinsic SPDE using the covariance-based rational SPDE approximation

The function samples a Gaussian random field based using the
covariance-based rational SPDE approximation.

## Usage

``` r
# S3 method for class 'intrinsicCBrSPDEobj'
simulate(
  object,
  nsim = 1,
  seed = NULL,
  kappa = NULL,
  tau = NULL,
  alpha = NULL,
  beta = NULL,
  integral.constraint = TRUE,
  use_kl = NULL,
  ...
)
```

## Arguments

- object:

  The covariance-based rational SPDE approximation, computed using
  [`intrinsic.matern.operators()`](https://davidbolin.github.io/rSPDE/reference/intrinsic.matern.operators.md)

- nsim:

  The number of simulations.

- seed:

  An object specifying if and how the random number generator should be
  initialized (‘seeded’).

- kappa:

  new value of kappa to use

- tau:

  new value of tau to use

- alpha:

  new value of alpha to use

- beta:

  new value of beta to use

- integral.constraint:

  Should the contraint on the integral be done?

- use_kl:

  Simulate based on a KL expansion?

- ...:

  Currently not used.

## Value

A matrix with the `nsim` samples as columns.
