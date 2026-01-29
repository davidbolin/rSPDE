# Summarise CBrSPDE objects

Summary method for class "CBrSPDEobj"

## Usage

``` r
# S3 method for class 'CBrSPDEobj'
summary(object, ...)

# S3 method for class 'summary.CBrSPDEobj'
print(x, ...)

# S3 method for class 'CBrSPDEobj'
print(x, ...)
```

## Arguments

- object:

  an object of class "CBrSPDEobj", usually, a result of a call to
  [`matern.operators()`](https://davidbolin.github.io/rSPDE/reference/matern.operators.md).

- ...:

  further arguments passed to or from other methods.

- x:

  an object of class "summary.CBrSPDEobj", usually, a result of a call
  to `summary.CBrSPDEobj()`.

## Examples

``` r
# Compute the covariance-based rational approximation of a
# Gaussian process with a Matern covariance function on R
kappa <- 10
sigma <- 1
nu <- 0.8
range <- sqrt(8 * nu) / kappa

# create mass and stiffness matrices for a FEM discretization
x <- seq(from = 0, to = 1, length.out = 101)
fem <- rSPDE.fem1d(x)

# compute rational approximation of covariance function at 0.5
tau <- sqrt(gamma(nu) / (sigma^2 * kappa^(2 * nu) *
  (4 * pi)^(1 / 2) * gamma(nu + 1 / 2)))
op_cov <- matern.operators(
  loc_mesh = x, nu = nu,
  range = range, sigma = sigma, d = 1, m = 2,
  parameterization = "matern"
)

op_cov
#> Type of approximation:  Covariance-Based Matern SPDE Approximation 
#> Parameterization:  matern 
#> Type of rational approximation:  brasil 
#> Parameters of covariance function: range =  0.2529822 , sigma =  1 , nu =  0.8 
#> Order or rational approximation:  2 
#> Size of discrete operators:  101  x  101 
#> Stationary Model
```
