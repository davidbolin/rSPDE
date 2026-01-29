# Optimized precision matrix of the covariance-based rational approximation

`rspde.matern.precision` is used for computing the optimized version of
the precision matrix of the covariance-based rational SPDE approximation
of a stationary Gaussian random fields on \\R^d\\ with a Matern
covariance function \$\$C(h) =
\frac{\sigma^2}{2^{\nu-1}\Gamma(\nu)}(\kappa h)^\nu K\_\nu(\kappa
h).\$\$

## Usage

``` r
rspde.matern.precision.opt(
  kappa,
  nu,
  tau,
  rspde.order,
  dim,
  fem_matrices,
  graph = NULL,
  sharp,
  type_rational_approx
)
```

## Arguments

- kappa:

  Range parameter of the covariance function.

- nu:

  Shape parameter of the covariance function.

- tau:

  Scale parameter of the covariance function.

- rspde.order:

  The order of the rational approximation

- dim:

  The dimension of the domain

- fem_matrices:

  A list containing the FEM-related matrices. The list should contain
  elements C, G, G_2, G_3, etc.

- graph:

  The sparsity graph of the matrices. If NULL, only a vector of the
  elements will be returned, if non-NULL, a sparse matrix will be
  returned.

- sharp:

  The sparsity graph should have the correct sparsity (costs more to
  perform a sparsity analysis) or an upper bound for the sparsity?

- type_rational_approx:

  Which type of rational approximation should be used? The current types
  are "brasil", "chebfun" or "chebfunLB".

## Value

The precision matrix
