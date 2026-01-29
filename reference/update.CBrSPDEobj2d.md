# Update parameters of CBrSPDEobj2d objects

Function to change the parameters of a CBrSPDEobj object

## Usage

``` r
# S3 method for class 'CBrSPDEobj2d'
update(
  object,
  hx = NULL,
  hy = NULL,
  hxy = NULL,
  sigma = NULL,
  nu = NULL,
  m = NULL,
  ...
)
```

## Arguments

- object:

  The covariance-based rational SPDE approximation, computed using
  [`matern2d.operators()`](https://davidbolin.github.io/rSPDE/reference/matern2d.operators.md)

- hx:

  If non-null, update the hx parameter.

- hy:

  If non-null, update the hy parameter.

- hxy:

  If non-null, update the hxy parameter.

- sigma:

  If non-null, update the standard deviation of the covariance function.

- nu:

  If non-null, update the shape parameter of the covariance function.
  Will be used if parameterization is 'matern'.

- m:

  If non-null, update the order of the rational approximation, which
  needs to be a positive integer.

- ...:

  Currently not used.

## Value

It returns an object of class "CBrSPDEobj2d.

## See also

[`simulate.CBrSPDEobj2d()`](https://davidbolin.github.io/rSPDE/reference/simulate.CBrSPDEobj2d.md),
[`matern2d.operators()`](https://davidbolin.github.io/rSPDE/reference/matern2d.operators.md)

## Examples

``` r
library(fmesher)
n_loc <- 2000
loc_2d_mesh <- matrix(runif(n_loc * 2), n_loc, 2)
mesh_2d <- fm_mesh_2d(loc = loc_2d_mesh, cutoff = 0.03, max.edge = c(0.1, 0.5))
op <- matern2d.operators(mesh = mesh_2d)
op <- update(op, nu = 0.5)
```
