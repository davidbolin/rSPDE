# Optimized precision matrix of stationary Gaussian Matern random fields with integer covariance exponent

`rspde.matern.precision.integer.opt` is used for computing the optimized
version of the precision matrix of stationary Gaussian random fields on
\\R^d\\ with a Matern covariance function \$\$C(h) =
\frac{\sigma^2}{2^{\nu-1}\Gamma(\nu)}(\kappa h)^\nu K\_\nu(\kappa
h),\$\$ where \\\alpha = \nu + d/2\\ is a natural number.

## Usage

``` r
rspde.matern.precision.integer.opt(
  kappa,
  nu,
  tau,
  d,
  fem_matrices,
  graph = NULL
)
```

## Arguments

- kappa:

  Range parameter of the covariance function.

- nu:

  Shape parameter of the covariance function.

- tau:

  Scale parameter of the covariance function.

- d:

  The dimension of the domain

- fem_matrices:

  A list containing the FEM-related matrices. The list should contain
  elements C, G, G_2, G_3, etc.

- graph:

  The sparsity graph of the matrices. If NULL, only a vector of the
  elements will be returned, if non-NULL, a sparse matrix will be
  returned.

## Value

The precision matrix
