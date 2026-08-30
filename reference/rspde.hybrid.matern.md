# Hybrid Whittle-Matern SPDE model for INLA / inlabru (alpha = 2)

Creates an INLA cgeneric latent model for the hybrid Whittle-Matern SPDE
\$\$Y(s) = \beta_X^T L^{-1} X(s) + \tau^{-1} L^{-1} W(s),\$\$ where \\L
= \kappa^2 - \Delta\\ and the deterministic mean component \\\beta_X^T
L^{-1} X\\ is part of the latent model itself: the cgeneric stores the
centred field \\u(s) = \tau^{-1} L^{-1} W(s)\\ and returns \\\mu(s) =
\beta_X^T L^{-1} X(s)\\ via the `mu` callback. The coefficients
`beta_X`, together with `tau` and `kappa`, are estimated jointly as
cgeneric hyperparameters – this is the genuine fully Bayesian fit, not a
plug-in approximation.

## Usage

``` r
rspde.hybrid.matern(
  mesh,
  X,
  prior.tau = NULL,
  prior.kappa = NULL,
  prior.beta_x = NULL,
  start.ltau = NULL,
  start.lkappa = NULL,
  start.beta_x = NULL,
  ...,
  separate_kappa_mu = FALSE,
  prior.kappa_mu = NULL,
  start.lkappa_mu = NULL,
  debug = FALSE,
  shared_lib = "detect"
)
```

## Arguments

- mesh:

  An `fmesher` mesh (1d or 2d).

- X:

  Numeric vector or matrix of covariate fields evaluated at the mesh
  nodes; one column per covariate.

- prior.tau:

  A list with `mean` and `prec` (Gaussian prior on `log(tau)`).
  Defaults: mean = 0, prec = 0.1.

- prior.kappa:

  A list with `mean` and `prec` (Gaussian prior on `log(kappa)`).
  Defaults: mean derived from the mesh range, prec = 0.1.

- prior.beta_x:

  A list with `mean` (vector of length `p`) and `prec` (vector of length
  `p`). Independent Gaussian priors on the regression coefficients.
  Defaults: mean = 0, prec = 0.001 (uninformative).

- start.ltau, start.lkappa:

  Starting values for `log(tau)` and `log(kappa)`. Default to the prior
  means.

- start.beta_x:

  Numeric vector of starting values for `beta_X`. Default: zeros.

- ...:

  Additional arguments passed to the underlying model constructor.

- separate_kappa_mu:

  Logical. If `TRUE`, the deterministic mean operator uses its own range
  parameter `kappa_mu`, estimated as an extra hyperparameter, rather
  than sharing `kappa` with the random field. Default: `FALSE`.

- prior.kappa_mu:

  A list with `mean` and `prec` (Gaussian prior on `log(kappa_mu)`).
  Only used when `separate_kappa_mu = TRUE`. Defaults to the same prior
  as `prior.kappa`.

- start.lkappa_mu:

  Starting value for `log(kappa_mu)`. Only used when
  `separate_kappa_mu = TRUE`. Defaults to the prior mean.

- debug:

  Passed to INLA.

- shared_lib:

  Which shared lib to use for the cgeneric. See
  [`rspde.matern()`](https://davidbolin.github.io/rSPDE/reference/rspde.matern.md)
  for details.

## Value

An INLA cgeneric latent model that can be used in INLA formulae as
`f(idx, model = <model>)`.

## Details

Currently the implementation is restricted to the case alpha = 2 (which
is the main case of practical interest: in 2d this is the Matern with
`nu = 1`).

The cgeneric hyperparameters are \$\$\theta = (\log\tau,\\ \log\kappa,\\
\beta\_{X,1}, \ldots, \beta\_{X,p}).\$\$ The model stores the FEM `C`
and `G` matrices and the covariate matrix `X` and (re)computes \\Q =
\tau^2\\ L\\ C^{-1}\\ L\\ and \\\mu = L^{-1}(X \beta_X)\\ at every INLA
call. The mean computation is a single sparse Cholesky solve done with
Eigen.

## See also

[`hybrid.spde()`](https://davidbolin.github.io/rSPDE/reference/hybrid.spde.md),
[`rspde.matern()`](https://davidbolin.github.io/rSPDE/reference/rspde.matern.md)

## Examples

``` r
library(INLA)
x <- seq(0, 1, length.out = 81)
mesh <- fmesher::fm_mesh_1d(x)
X <- cbind(sin(2 * pi * x), cos(2 * pi * x))
hyb <- rspde.hybrid.matern(mesh = mesh, X = X)
# ... use in formula: y ~ -1 + f(idx, model = hyb)
```
