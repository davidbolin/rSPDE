# Kriging prediction for a hybrid Whittle-Matern SPDE model

Computes posterior predictions for the hybrid model \\Y_i = (A Y)\_i +
\epsilon_i\\, where \\Y\\ satisfies \\L^{\alpha/2}(\tau Y) = \beta X +
W\\. The posterior of the centred latent field \\Y - \mu\\ is computed
using the underlying SPDE prediction routine; the deterministic mean
\\\mu\\ is then added back at the prediction locations.

## Usage

``` r
# S3 method for class 'hybrid_spde'
predict(
  object,
  A,
  Aprd,
  Y,
  sigma.e,
  compute.variances = FALSE,
  posterior_samples = FALSE,
  n_samples = 100,
  only_latent = FALSE,
  ...
)
```

## Arguments

- object:

  A `hybrid_spde` object.

- A:

  Projection matrix linking the observations to the mesh nodes.

- Aprd:

  Projection matrix linking the prediction locations to the mesh nodes.

- Y:

  Observations.

- sigma.e:

  Standard deviation of the measurement noise.

- compute.variances:

  If TRUE, the kriging variances are returned.

- posterior_samples:

  If TRUE, posterior samples are returned.

- n_samples:

  Number of posterior samples.

- only_latent:

  If TRUE, the posterior samples are not perturbed by the measurement
  noise.

- ...:

  Additional arguments passed to the underlying predict method.

## Value

A list with components `mean`, optionally `variance` and `samples`.

## See also

[`hybrid.spde()`](https://davidbolin.github.io/rSPDE/reference/hybrid.spde.md)
