# Hybrid non-stationary Whittle-Matern SPDE model with non-zero mean

`hybrid.spde` constructs a finite element/rational SPDE approximation of
the model \$\$Y(s) = \beta L^{-\alpha/2} X(s) + \tau^{-1} L^{-\alpha/2}
W(s),\$\$ where \\L = \kappa(s)^2 - \Delta\\ is the (possibly
non-stationary) Whittle-Matern differential operator, \\\tau\\ is a
positive scale that affects only the random component, \\X\\ is a
(matrix of) covariate field(s) supplied at the mesh nodes, \\\beta\\ is
a vector of regression coefficients and \\W\\ is Gaussian white noise.

## Usage

``` r
hybrid.spde(
  X = NULL,
  beta_X = NULL,
  kappa = NULL,
  tau = NULL,
  kappa_mu = NULL,
  theta = NULL,
  B.tau = NULL,
  B.kappa = NULL,
  B.sigma = NULL,
  B.range = NULL,
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

- X:

  Numeric vector or matrix containing the covariate field(s) evaluated
  at the mesh nodes. If a matrix is supplied, each column is treated as
  a separate covariate field. The number of rows must equal the number
  of mesh nodes of the underlying FEM discretization.

- beta_X:

  Numeric vector of regression coefficients corresponding to the columns
  of `X`. Defaults to a vector of zeros of length `ncol(X)`.

- kappa, tau, theta, B.tau, B.kappa, B.sigma, B.range, alpha, nu,
  parameterization, G, C, d, graph, mesh, range_mesh, loc_mesh, m, type,
  type_rational_approximation, check_stationarity:

  Arguments passed to
  [`spde.matern.operators()`](https://davidbolin.github.io/rSPDE/reference/spde.matern.operators.md).
  The default `check_stationarity = TRUE` allows the underlying
  covariance model to be returned as a stationary
  [`matern.operators()`](https://davidbolin.github.io/rSPDE/reference/matern.operators.md)
  object whenever `tau` and `kappa` are constant after the log-linear
  regression on `theta` is applied.

- kappa_mu:

  Optional. If `NULL` (the default), the operator \\L\_\mu = \kappa^2 -
  \Delta\\ used for the deterministic mean \\\mu = \beta
  L\_\mu^{-\alpha/2} X\\ shares the range parameter `kappa` with the
  random component. If a positive scalar (or vector of length n_mesh for
  non-stationary models) is supplied, \\L\_\mu = \kappa\_\mu^2 -
  \Delta\\ is used instead.

## Value

`hybrid.spde` returns an object of class `"hybrid_spde"` which also
inherits from the class of the object returned by
[`spde.matern.operators()`](https://davidbolin.github.io/rSPDE/reference/spde.matern.operators.md)
(`spde_matern_operator` plus either `rSPDEobj` or `CBrSPDEobj` depending
on `type`). The object stores the covariate matrix `X`, the regression
coefficients `beta_X`, the cached operator object `op_mean` used to
evaluate \\L^{-\alpha/2}\\, and the evaluated mean function `mu` at the
mesh nodes.

## See also

[`spde.matern.operators()`](https://davidbolin.github.io/rSPDE/reference/spde.matern.operators.md),
[`simulate.hybrid_spde()`](https://davidbolin.github.io/rSPDE/reference/simulate.hybrid_spde.md),
[`predict.hybrid_spde()`](https://davidbolin.github.io/rSPDE/reference/predict.hybrid_spde.md)

## Examples

``` r
set.seed(1)
x <- seq(from = 0, to = 1, length.out = 51)
kappa <- 10
tau <- 0.5
alpha <- 1.3
X <- cbind(sin(2 * pi * x), cos(2 * pi * x))
beta_X <- c(1.5, -0.5)
op <- hybrid.spde(
  X = X, beta_X = beta_X,
  kappa = kappa, tau = tau, alpha = alpha,
  loc_mesh = x, d = 1, parameterization = "spde"
)
Y <- simulate(op, nsim = 1)
```
