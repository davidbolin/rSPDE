# Get the precision matrix of spacetimeobj objects

Function to get the precision matrix of a spacetimeobj object

## Usage

``` r
# S3 method for class 'spacetimeobj'
precision(object, kappa = NULL, sigma = NULL, gamma = NULL, rho = NULL, ...)
```

## Arguments

- object:

  The model object computed using
  [`spacetime.operators()`](https://davidbolin.github.io/rSPDE/reference/spacetime.operators.md)

- kappa:

  If non-null, update the range parameter of the covariance function.

- sigma:

  If non-null, update the standard deviation of the covariance function.

- gamma:

  If non-null, update the temporal range parameter of the covariance
  function.

- rho:

  If non-null, update the drift parameter of the covariance function.

- ...:

  Currently not used.

## Value

The precision matrix.

## See also

[`simulate.spacetimeobj()`](https://davidbolin.github.io/rSPDE/reference/simulate.spacetimeobj.md),
[`spacetime.operators()`](https://davidbolin.github.io/rSPDE/reference/spacetime.operators.md)

## Examples

``` r
s <- seq(from = 0, to = 20, length.out = 101)
t <- seq(from = 0, to = 20, length.out = 31)

op_cov <- spacetime.operators(space_loc = s, time_loc = t,
                             kappa = 5, sigma = 10, alpha = 1,
                             beta = 2, rho = 1, gamma = 0.05)
prec_matrix <- precision(op_cov)
```
