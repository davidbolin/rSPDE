# Rational approximations of stationary anisotropic Gaussian Matern random fields

`matern2d.operators` is used for computing a rational SPDE approximation
of a stationary Gaussian random fields on \\R^d\\ with a Matern
covariance function \$\$C(h) =
\frac{\sigma^2}{2^{\nu-1}\Gamma(\nu)}(\sqrt{h^T H^{-1}h})^\nu
K\_\nu(\sqrt{h^T H^{-1}h})\$\$, based on a SPDE representation of the
form \$\$(I - \nabla\cdot(H\nabla))^{(\nu+1)/2}u = c\sigma W\$\$, where
\$c\>0\$ is a constant. The matrix \\H\\ is defined as
\$\$\begin{bmatrix} h_x^2 & h_xh_yh\_{xy} \\ h_xh_yh\_{xy} & h_y^2
\end{bmatrix}\$\$

## Usage

``` r
matern2d.operators(
  hx = NULL,
  hy = NULL,
  hxy = NULL,
  nu = NULL,
  sigma = NULL,
  mesh = NULL,
  fem = NULL,
  m = 1,
  type_rational_approximation = c("brasil", "chebfun", "chebfunLB"),
  return_fem_matrices = FALSE
)
```

## Arguments

- hx:

  Parameter in the H matrix.

- hy:

  Parameter in the H matrix.

- hxy:

  Parameter in the H matrix.

- nu:

  Smoothness parameter.

- sigma:

  standard deviation parameter.

- mesh:

  An `fmesher` mesh.

- fem:

  Optional precomputed FEM matrices.

- m:

  The order of the rational approximation, which needs to be a positive
  integer. The default value is 1.

- type_rational_approximation:

  Which type of rational approximation should be used? The current types
  are "brasil", "chebfun" or "chebfunLB".

- return_fem_matrices:

  Should the FEM matrices be returned?

## Value

An object of type `CBrSPDEobj2d`

## See also

[`fractional.operators()`](https://davidbolin.github.io/rSPDE/reference/fractional.operators.md),
[`spde.matern.operators()`](https://davidbolin.github.io/rSPDE/reference/spde.matern.operators.md),
[`matern.operators()`](https://davidbolin.github.io/rSPDE/reference/matern.operators.md)

## Examples

``` r
library(fmesher)
n_loc <- 2000
loc_2d_mesh <- matrix(runif(n_loc * 2), n_loc, 2)
mesh_2d <- fm_mesh_2d(loc = loc_2d_mesh, cutoff = 0.03, max.edge = c(0.1, 0.5))
op <- matern2d.operators(mesh = mesh_2d)
```
