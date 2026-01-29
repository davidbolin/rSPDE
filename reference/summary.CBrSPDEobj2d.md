# Summarise CBrSPDEobj2d objects

Summary method for class "CBrSPDEobj2d"

## Usage

``` r
# S3 method for class 'CBrSPDEobj2d'
summary(object, ...)

# S3 method for class 'summary.CBrSPDEobj2d'
print(x, ...)

# S3 method for class 'CBrSPDEobj2d'
print(x, ...)
```

## Arguments

- object:

  an object of class "CBrSPDEobj2d", usually, a result of a call to
  [`matern2d.operators()`](https://davidbolin.github.io/rSPDE/reference/matern2d.operators.md).

- ...:

  further arguments passed to or from other methods.

- x:

  an object of class "summary.CBrSPDEobj2d", usually, a result of a call
  to `summary.CBrSPDEobj2d()`.

## Examples

``` r
library(fmesher)
n_loc <- 2000
loc_2d_mesh <- matrix(runif(n_loc * 2), n_loc, 2)
mesh_2d <- fm_mesh_2d(loc = loc_2d_mesh, cutoff = 0.03, max.edge = c(0.1, 0.5))
op <- matern2d.operators(mesh = mesh_2d)
op
#> Type of approximation:  Covariance-Based aniostropic Matern SPDE Approximation 
#> Type of rational approximation:  brasil 
#> Parameters of covariance function: sigma =  1 , hx =  0.3157421 , hy =  0.3157421 , hxy =  0 , nu =  1 
#> Order or rational approximation:  1 
#> Size of discrete operators:  1  x  1 
```
