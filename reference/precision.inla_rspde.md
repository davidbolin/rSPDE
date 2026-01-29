# Get the precision matrix of `inla_rspde` objects

Function to get the precision matrix of an `inla_rspde` object created
with the
[`rspde.matern()`](https://davidbolin.github.io/rSPDE/reference/rspde.matern.md)
function.

## Usage

``` r
# S3 method for class 'inla_rspde'
precision(object, theta = NULL, ...)
```

## Arguments

- object:

  The `inla_rspde` object obtained with the
  [`rspde.matern()`](https://davidbolin.github.io/rSPDE/reference/rspde.matern.md)
  function.

- theta:

  If null, the starting values for theta will be used. Otherwise, it
  must be suplied as a vector. For stationary models, we have
  `theta = c(log(tau), log(kappa), nu)`. For nonstationary models, we
  have `theta = c(theta_1, theta_2, ..., theta_n, nu)`.

- ...:

  Currently not used.

## Value

The precision matrix.

## See also

[`precision.CBrSPDEobj()`](https://davidbolin.github.io/rSPDE/reference/precision.CBrSPDEobj.md),
[`matern.operators()`](https://davidbolin.github.io/rSPDE/reference/matern.operators.md)
