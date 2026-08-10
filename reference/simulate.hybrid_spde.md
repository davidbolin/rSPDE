# Simulation of a hybrid Whittle-Matern SPDE model

Samples the hybrid SPDE model \$\$L^{\alpha/2}(\tau Y(s)) = \beta X(s) +
W(s)\$\$ by first sampling the mean-zero part using the simulation
method of the underlying
[`spde.matern.operators()`](https://davidbolin.github.io/rSPDE/reference/spde.matern.operators.md)
approximation and then adding the deterministic mean \\\mu = \tau^{-1}
\beta L^{-\alpha/2} X\\.

## Usage

``` r
# S3 method for class 'hybrid_spde'
simulate(object, nsim = 1, seed = NULL, ...)
```

## Arguments

- object:

  A `hybrid_spde` object returned by
  [`hybrid.spde()`](https://davidbolin.github.io/rSPDE/reference/hybrid.spde.md).

- nsim:

  Number of samples to generate.

- seed:

  Optional integer used to initialise the random number generator.

- ...:

  Additional arguments passed to the underlying simulate method (e.g.
  updated parameter values).

## Value

A matrix with `nsim` columns and one row per mesh node.

## See also

[`hybrid.spde()`](https://davidbolin.github.io/rSPDE/reference/hybrid.spde.md)
