# Rational approximations of non-stationary Gaussian SPDE Matern random fields

`spde.matern.operators` is used for computing a rational SPDE
approximation of a Gaussian random fields on \\R^d\\ defined as a
solution to the SPDE \$\$(\kappa(s) - \Delta)^\beta (\tau(s)u(s)) =
W.\$\$

## Usage

``` r
spde.matern.operators(
  kappa = NULL,
  tau = NULL,
  theta = NULL,
  B.tau = matrix(c(0, 1, 0), 1, 3),
  B.kappa = matrix(c(0, 0, 1), 1, 3),
  B.sigma = matrix(c(0, 1, 0), 1, 3),
  B.range = matrix(c(0, 0, 1), 1, 3),
  alpha = NULL,
  nu = NULL,
  parameterization = c("spde", "matern"),
  G = NULL,
  C = NULL,
  d = NULL,
  graph = NULL,
  mesh = NULL,
  range_mesh = NULL,
  loc_mesh = NULL,
  m = 1,
  type = c("covariance", "operator"),
  type_rational_approximation = c("brasil", "chebfun", "chebfunLB"),
  check_stationarity = TRUE
)
```

## Arguments

- kappa:

  Vector with the, possibly spatially varying, range parameter evaluated
  at the locations of the mesh used for the finite element
  discretization of the SPDE.

- tau:

  Vector with the, possibly spatially varying, precision parameter
  evaluated at the locations of the mesh used for the finite element
  discretization of the SPDE.

- theta:

  Theta parameter that connects B.tau and B.kappa to tau and kappa
  through a log-linear regression, in case the parameterization is
  `spde`, and that connects B.sigma and B.range to tau and kappa in case
  the parameterization is `matern`. When both tau and kappa are constant
  after this transformation, `spde.matern.operators()` delegates to
  [`matern.operators()`](https://davidbolin.github.io/rSPDE/reference/matern.operators.md).

- B.tau:

  Matrix with specification of log-linear model for \\\tau\\. Will be
  used if `parameterization = 'spde'`.

- B.kappa:

  Matrix with specification of log-linear model for \\\kappa\\. Will be
  used if `parameterization = 'spde'`.

- B.sigma:

  Matrix with specification of log-linear model for \\\sigma\\. Will be
  used if `parameterization = 'matern'`.

- B.range:

  Matrix with specification of log-linear model for \\\rho\\, which is a
  range-like parameter (it is exactly the range parameter in the
  stationary case). Will be used if `parameterization = 'matern'`.

- alpha:

  smoothness parameter. Will be used if the parameterization is 'spde'.

- nu:

  Shape parameter of the covariance function. Will be used if the
  parameterization is 'matern'.

- parameterization:

  Which parameterization to use? `matern` uses range, std. deviation and
  nu (smoothness). `spde` uses kappa, tau and nu (smoothness). The
  default is `matern`.

- G:

  The stiffness matrix of a finite element discretization of the domain
  of interest.

- C:

  The mass matrix of a finite element discretization of the domain of
  interest.

- d:

  The dimension of the domain. Does not need to be given if `mesh` is
  used.

- graph:

  An optional `metric_graph` object. Replaces `d`, `C` and `G`.

- mesh:

  An optional inla mesh. `d`, `C` and `G` must be given if `mesh` is not
  given.

- range_mesh:

  The range of the mesh. Will be used to provide starting values for the
  parameters. Will be used if `mesh` and `graph` are `NULL`, and if one
  of the parameters (kappa or tau for spde parameterization, or sigma or
  range for matern parameterization) are not provided.

- loc_mesh:

  The mesh locations used to construct the matrices C and G. This option
  should be provided if one wants to use the
  [`rspde_lme()`](https://davidbolin.github.io/rSPDE/reference/rspde_lme.md)
  function and will not provide neither graph nor mesh. Only works for
  1d data. Does not work for metric graphs. For metric graphs you should
  supply the graph using the `graph` argument.

- m:

  The order of the rational approximation, which needs to be a positive
  integer. The default value is 1.

- type:

  The type of the rational approximation. The options are "covariance"
  and "operator". The default is "covariance".

- type_rational_approximation:

  Which type of rational approximation should be used? The current types
  are "brasil", "chebfun" or "chebfunLB".

- check_stationarity:

  Logical; if TRUE, automatically returns a stationary model when
  tau/kappa (or sigma/range) are constant. Set to FALSE to keep a
  non-stationary model even when parameters are constant.

## Value

`spde.matern.operators` returns an object of class "rSPDEobj. This
object contains the quantities listed in the output of
[`fractional.operators()`](https://davidbolin.github.io/rSPDE/reference/fractional.operators.md)
as well as the smoothness parameter \\\nu\\.

## Details

The approximation is based on a rational approximation of the fractional
operator \\(\kappa(s)^2 -\Delta)^\beta\\, where \\\beta = (\nu +
d/2)/2\\. This results in an approximate model on the form \$\$P_l u(s)
= P_r W,\$\$ where \\P_j = p_j(L)\\ are non-fractional operators defined
in terms of polynomials \\p_j\\ for \\j=l,r\\. The order of \\p_r\\ is
given by `m` and the order of \\p_l\\ is \\m + m\_\beta\\ where
\\m\_\beta\\ is the integer part of \\\beta\\ if \\\beta\>1\\ and
\\m\_\beta = 1\\ otherwise.

The discrete approximation can be written as \\u = P_r x\\ where \\x
\sim N(0,Q^{-1})\\ and \\Q = P_l^T C^{-1} P_l\\. Note that the matrices
\\P_r\\ and \\Q\\ may be be ill-conditioned for \\m\>1\\. In this case,
the metehods in
[`operator.operations()`](https://davidbolin.github.io/rSPDE/reference/operator.operations.md)
should be used for operations involving the matrices, since these
methods are more numerically stable.

## See also

[`fractional.operators()`](https://davidbolin.github.io/rSPDE/reference/fractional.operators.md),
`spde.matern.operators()`,
[`matern.operators()`](https://davidbolin.github.io/rSPDE/reference/matern.operators.md)

## Examples

``` r
# Sample non-stationary Matern field on R
tau <- 1
nu <- 0.8

# create mass and stiffness matrices for a FEM discretization
x <- seq(from = 0, to = 1, length.out = 101)
fem <- rSPDE.fem1d(x)

# define a non-stationary range parameter
kappa <- seq(from = 2, to = 20, length.out = length(x))
alpha <- nu + 1 / 2
# compute rational approximation
op <- spde.matern.operators(
  kappa = kappa, tau = tau, alpha = alpha,
  G = fem$G, C = fem$C, d = 1
)

# sample the field
u <- simulate(op)

# plot the sample
plot(x, u, type = "l", ylab = "u(s)", xlab = "s")

```
