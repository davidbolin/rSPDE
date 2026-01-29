# Update parameters of intrinsicCBrSPDEobj objects

Function to change the parameters of a intrinsicCBrSPDEobj object

## Usage

``` r
# S3 method for class 'intrinsicCBrSPDEobj'
update(object, kappa = NULL, tau = NULL, alpha = NULL, beta = NULL, ...)
```

## Arguments

- object:

  Model object created by
  [`intrinsic.matern.operators()`](https://davidbolin.github.io/rSPDE/reference/intrinsic.matern.operators.md)

- kappa:

  kappa value to be updated.

- tau:

  tau value to be update.

- alpha:

  alpha value to be updated.

- beta:

  beta value to be updated. .

- ...:

  currently not used.

## Value

An object of type intrinsicCBrSPDEobj with updated parameters.

## Examples

``` r
if (requireNamespace("RSpectra", quietly = TRUE)) {
  x <- seq(from = 0, to = 10, length.out = 201)
  beta <- 1
  alpha <- 1
  kappa <- 1
  op <- intrinsic.matern.operators(
    kappa = kappa, tau = 1, alpha = alpha,
    beta = beta, loc_mesh = x, d = 1
  )
op <- update(op, beta = 1.1, alpha = 0.9) 
}
```
