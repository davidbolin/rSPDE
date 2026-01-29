# Rational approximations of stationary anisotropic Gaussian Matern random fields

`rspde.anistropic2d` computes a Finite Element Method (FEM)
approximation of a Gaussian random field defined as the solution to the
stochastic partial differential equation (SPDE): \$\$C(h) =
\frac{\sigma^2}{2^{\nu-1}\Gamma(\nu)}(\sqrt{h^T H^{-1}h})^\nu
K\_\nu(\sqrt{h^T H^{-1}h})\$\$, based on a SPDE representation of the
form \$\$(I - \nabla\cdot(H\nabla))^{(\nu+1)/2}u = c\sigma W\$\$, where
\$c\>0\$ is a constant. The matrix \\H\\ is defined as
\$\$\begin{bmatrix} h_x^2 & h_xh_yh\_{xy} \\ h_xh_yh\_{xy} & h_y^2
\end{bmatrix}\$\$

## Usage

``` r
rspde.anistropic2d(
  mesh,
  nu = NULL,
  nu.upper.bound = 2,
  rspde.order = 1,
  prior.hx = NULL,
  prior.hy = NULL,
  prior.hxy = NULL,
  prior.sigma = NULL,
  prior.precision = NULL,
  prior.nu = NULL,
  prior.nu.dist = "lognormal",
  nu.prec.inc = 0.01,
  type.rational.approx = "brasil",
  shared_lib = "detect",
  debug = FALSE,
  ...
)
```

## Arguments

- mesh:

  Spatial mesh for the FEM approximation.

- nu:

  If nu is set to a parameter, nu will be kept fixed and will not be
  estimated. If nu is `NULL`, it will be estimated.

- nu.upper.bound:

  Upper bound for the smoothness parameter \\\nu\\. If `NULL`, it will
  be set to 2.

- rspde.order:

  The order of the covariance-based rational SPDE approach. The default
  order is 1.

- prior.hx:

  A list specifying the prior for the parameter \\h_x\\ in the matrix
  \\H\\. This list may contain two elements: `mean` and/or `precision`,
  both of which must be numeric scalars. The precision refers to the
  prior on \\\log(h_x)\\. If `NULL`, default values will be used. The
  `mean` value is also used as starting value for hx.

- prior.hy:

  A list specifying the prior for the parameter \\h_y\\ in the matrix
  \\H\\. This list may contain two elements: `mean` and/or `precision`,
  both of which must be numeric scalars. The precision refers to the
  prior on \\\log(h_x)\\. If `NULL`, default values will be used. The
  `mean` value is also used as starting value for hy.

- prior.hxy:

  A list specifying the prior for the parameter \\h_x\\ in the matrix
  \\H\\. This list may contain two elements: `mean` and/or `precision`,
  both of which must be numeric scalars. The precision refers to the
  prior on \\\log((h\_{xy}+1)/(1-h\_{xy}))\\. If `NULL`, default values
  will be used. The `mean` value is also used as starting value for hxy.

- prior.sigma:

  A list specifying the prior for the variance parameter \\\sigma\\.
  This list may contain two elements: `mean` and/or `precision`, both of
  which must be numeric scalars. The precision refers to the prior on
  \\\log(\sigma)\\. If `NULL`, default values will be used. The `mean`
  value is also used as starting value for sigma.

- prior.precision:

  A precision matrix for \\\log(h_x), \log(h_y),
  \log((h\_{xy}+1)/(1-h\_{xy})), \log(\sigma)\\. This matrix replaces
  the precision element from `prior.kappa`, `prior.sigma`,
  `prior.gamma`, and `prior.rho` respectively. For dimension 1
  `prior.precision` must be a 4x4 matrix. For dimension 2, \\\rho\\ is a
  vector of length 2, so in this case `prior.precision` must be a 5x5
  matrix. If `NULL`, a diagonal precision matrix with default values
  will be used.

- prior.nu:

  a list containing the elements `mean` and `prec` for beta
  distribution, or `loglocation` and `logscale` for a truncated
  lognormal distribution. `loglocation` stands for the location
  parameter of the truncated lognormal distribution in the log scale.
  `prec` stands for the precision of a beta distribution. `logscale`
  stands for the scale of the truncated lognormal distribution on the
  log scale. Check details below.

- prior.nu.dist:

  The distribution of the smoothness parameter. The current options are
  "beta" or "lognormal". The default is "lognormal".

- nu.prec.inc:

  Amount to increase the precision in the beta prior distribution. Check
  details below.

- type.rational.approx:

  Which type of rational approximation should be used? The current types
  are "brasil", "chebfun" or "chebfunLB".

- shared_lib:

  String specifying which shared library to use for the Cgeneric
  implementation. Options are "detect", "INLA", or "rSPDE". You may also
  specify the direct path to a .so (or .dll) file.

- debug:

  Logical value indicating whether to enable INLA debug mode.

- ...:

  Additional arguments passed internally for configuration purposes.

## Value

An object of class `inla_rspde_spacetime` representing the FEM
approximation of the space-time Gaussian random field.

## Examples

``` r
library(fmesher)
n_loc <- 2000
loc_2d_mesh <- matrix(runif(n_loc * 2), n_loc, 2)
mesh_2d <- fm_mesh_2d(loc = loc_2d_mesh, cutoff = 0.03, max.edge = c(0.1, 0.5))
aniso_model <- rspde.anistropic2d(mesh = mesh_2d)
```
