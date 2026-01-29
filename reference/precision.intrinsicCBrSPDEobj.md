# Get the precision matrix of intrinsicCBrSPDEobj objects

Function to get the precision matrix of a intrinsicCBrSPDEobj object

## Usage

``` r
# S3 method for class 'intrinsicCBrSPDEobj'
precision(
  object,
  kappa = NULL,
  tau = NULL,
  alpha = NULL,
  beta = NULL,
  ld = FALSE,
  ...
)
```

## Arguments

- object:

  The model object computed using
  [`intrinsic.matern.operators()`](https://davidbolin.github.io/rSPDE/reference/intrinsic.matern.operators.md)

- kappa:

  If non-null, update the range parameter.

- tau:

  If non-null, update the precision parameter.

- alpha:

  If non-null, update the alpha parameter.

- beta:

  If non-null, update the beta parameter.

- ld:

  If TRUE, return the log determinant of the precision matrix instead of
  the precision matrix. By default FALSE.

- ...:

  Currently not used.

## Value

The precision matrix.

## See also

[`simulate.intrinsicCBrSPDEobj()`](https://davidbolin.github.io/rSPDE/reference/simulate.intrinsicCBrSPDEobj.md),
[`intrinsic.matern.operators()`](https://davidbolin.github.io/rSPDE/reference/intrinsic.matern.operators.md)

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
Q <- precision(op) 
}
```
