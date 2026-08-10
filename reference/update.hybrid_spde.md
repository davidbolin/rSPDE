# Update parameters of hybrid_spde objects

Updates the model parameters of a `hybrid_spde` object. Delegates the
update of the covariance part to the underlying
[`spde.matern.operators()`](https://davidbolin.github.io/rSPDE/reference/spde.matern.operators.md)
update method and additionally rebuilds the operator used to evaluate
the mean and recomputes the stored mean vector.

## Usage

``` r
# S3 method for class 'hybrid_spde'
update(object, X = NULL, beta_X = NULL, ..., kappa_mu = NULL)
```

## Arguments

- object:

  A `hybrid_spde` object.

- X:

  If non-null, update the covariate matrix.

- beta_X:

  If non-null, update the regression coefficients.

- ...:

  Additional arguments passed to the underlying update method of
  [`spde.matern.operators()`](https://davidbolin.github.io/rSPDE/reference/spde.matern.operators.md).

- kappa_mu:

  If non-null, update the range parameter of the operator used for the
  deterministic mean. Only meaningful when the model was constructed
  with a separate `kappa_mu`; if the model was built with
  `kappa_mu = NULL` (linked to `kappa`), any explicit update is honoured
  but the model remains in linked mode unless `separate_kappa_mu` is
  also set.

## Value

An object of class `"hybrid_spde"` with the updated parameters.

## See also

[`hybrid.spde()`](https://davidbolin.github.io/rSPDE/reference/hybrid.spde.md)
