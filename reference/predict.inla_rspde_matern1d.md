# Predict method for 'inlabru' stationary Matern 1d models

Auxiliar function to obtain predictions of the stationary Matern 1d
models using 'inlabru'.

## Usage

``` r
# S3 method for class 'inla_rspde_matern1d'
predict(
  object,
  cmp,
  bru_fit,
  newdata = NULL,
  formula = NULL,
  n.samples = 100,
  seed = 0L,
  probs = c(0.025, 0.5, 0.975),
  return_original_order = TRUE,
  num.threads = NULL,
  used = NULL,
  include = deprecated(),
  exclude = deprecated(),
  drop = FALSE,
  tolerance = 1e-04,
  ...
)
```

## Arguments

- object:

  An `inla_rspde_matern1d` object built with the
  [`rspde.matern1d()`](https://davidbolin.github.io/rSPDE/reference/rspde.matern1d.md)
  function.

- cmp:

  The 'inlabru' component used to fit the model. Only the `formula`
  input syntax is supported.

- bru_fit:

  A fitted model using 'inlabru' or 'INLA'.

- newdata:

  A data.frame of covariates needed for the prediction.

- formula:

  A formula where the right hand side defines an R expression to
  evaluate for each generated sample. If NULL, the latent and
  hyperparameter states are returned as named list elements. See Details
  for more information.

- n.samples:

  Integer setting the number of samples to draw in order to calculate
  the posterior statistics. The default is rather low but provides a
  quick approximate result.

- seed:

  Random number generator seed passed on to
  [`inla.posterior.sample()`](https://rdrr.io/pkg/INLA/man/posterior.sample.html)

- probs:

  A numeric vector of probabilities with values in the standard unit
  interval to be passed to stats::quantile

- return_original_order:

  Should the predictions be returned in the original order?

- num.threads:

  Specification of desired number of threads for parallel computations.
  Default NULL, leaves it up to 'INLA'. When seed != 0, overridden to
  "1:1"

- used:

  Optional
  [`inlabru::bru_used()`](https://inlabru-org.github.io/inlabru/reference/bru_used.html)
  specification for what components are needed by the predictor
  expression. Normally, autodetection works, and `used` can be left as
  `NULL` (the default).

- include, exclude:

  Deprecated unused alternatives to `used`.

- drop:

  logical; If keep=FALSE, data is a SpatialDataFrame, and the prediciton
  summary has the same number of rows as data, then the output is a
  SpatialDataFrame object. Default FALSE.

- tolerance:

  Tolerance for merging locations.

- ...:

  Additional arguments passed on to
  [`inla.posterior.sample()`](https://rdrr.io/pkg/INLA/man/posterior.sample.html).

## Value

A list with predictions.
