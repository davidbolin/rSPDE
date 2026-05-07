# Matern rSPDE model object for metric graphs in INLA

Creates an INLA object for a stationary Matern model on a metric graph
with general smoothness parameter.

## Usage

``` r
rspde.metric_graph(
  graph_obj,
  h = NULL,
  nu.upper.bound = 2,
  rspde.order = 1,
  nu = NULL,
  debug = FALSE,
  B.sigma = matrix(c(0, 1, 0), 1, 3),
  B.range = matrix(c(0, 0, 1), 1, 3),
  parameterization = c("matern", "spde"),
  B.tau = matrix(c(0, 1, 0), 1, 3),
  B.kappa = matrix(c(0, 0, 1), 1, 3),
  start.nu = NULL,
  start.theta = NULL,
  prior.nu = NULL,
  theta.prior.mean = NULL,
  theta.prior.prec = 0.1,
  prior.std.dev.nominal = 1,
  prior.range.nominal = NULL,
  prior.kappa.mean = NULL,
  prior.tau.mean = NULL,
  start.lstd.dev = NULL,
  start.lrange = NULL,
  start.ltau = NULL,
  start.lkappa = NULL,
  prior.theta.param = c("theta", "spde"),
  prior.nu.dist = c("lognormal", "beta"),
  nu.prec.inc = 1,
  type.rational.approx = c("brasil", "chebfun", "chebfunLB"),
  shared_lib = "detect"
)
```

## Arguments

- graph_obj:

  The graph object to build the model. Needs to be of class
  `metric_graph`. It should have a built mesh. If the mesh is not built,
  one will be built using h=0.01 as default.

- h:

  The width of the mesh in case the mesh was not built.

- nu.upper.bound:

  Upper bound for the smoothness parameter.

- rspde.order:

  The order of the covariance-based rational SPDE approach.

- nu:

  If nu is set to a parameter, nu will be kept fixed and will not be
  estimated. If nu is `NULL`, it will be estimated.

- debug:

  INLA debug argument

- B.sigma:

  Matrix with specification of log-linear model for \\\sigma\\. Will be
  used if `parameterization = 'matern'`.

- B.range:

  Matrix with specification of log-linear model for \\\rho\\, which is a
  range-like parameter (it is exactly the range parameter in the
  stationary case). Will be used if `parameterization = 'matern'`.

- parameterization:

  Which parameterization to use? `matern` uses range, std. deviation and
  nu (smoothness). `spde` uses kappa, tau and nu (smoothness). The
  default is `matern`.

- B.tau:

  Matrix with specification of log-linear model for \\\tau\\. Will be
  used if `parameterization = 'spde'`.

- B.kappa:

  Matrix with specification of log-linear model for \\\kappa\\. Will be
  used if `parameterization = 'spde'`.

- start.nu:

  Starting value for nu.

- start.theta:

  Starting values for the model parameters. In the stationary case, if
  `parameterization='matern'`, then `theta[1]` is the std.dev and
  `theta[2]` is the range parameter. If `parameterization = 'spde'`,
  then `theta[1]` is `tau` and `theta[2]` is `kappa`.

- prior.nu:

  a list containing the elements `mean` and `prec` for beta
  distribution, or `loglocation` and `logscale` for a truncated
  lognormal distribution. `loglocation` stands for the location
  parameter of the truncated lognormal distribution in the log scale.
  `prec` stands for the precision of a beta distribution. `logscale`
  stands for the scale of the truncated lognormal distribution on the
  log scale. Check details below.

- theta.prior.mean:

  A vector for the mean priors of `theta`.

- theta.prior.prec:

  A precision matrix for the prior of `theta`.

- prior.std.dev.nominal:

  Prior std. deviation to be used for the priors and for the starting
  values.

- prior.range.nominal:

  Prior range to be used for the priors and for the starting values.

- prior.kappa.mean:

  Prior kappa to be used for the priors and for the starting values.

- prior.tau.mean:

  Prior tau to be used for the priors and for the starting values.

- start.lstd.dev:

  Starting value for log of std. deviation. Will not be used if
  start.ltau is non-null. Will be only used in the stationary case and
  if `parameterization = 'matern'`.

- start.lrange:

  Starting value for log of range. Will not be used if start.lkappa is
  non-null. Will be only used in the stationary case and if
  `parameterization = 'matern'`.

- start.ltau:

  Starting value for log of tau. Will be only used in the stationary
  case and if `parameterization = 'spde'`.

- start.lkappa:

  Starting value for log of kappa. Will be only used in the stationary
  case and if `parameterization = 'spde'`.

- prior.theta.param:

  Should the lognormal prior be on `theta` or on the SPDE parameters
  (`tau` and `kappa` on the stationary case)?

- prior.nu.dist:

  The distribution of the smoothness parameter. The current options are
  "beta" or "lognormal". The default is "beta".

- nu.prec.inc:

  Amount to increase the precision in the beta prior distribution. Check
  details below.

- type.rational.approx:

  Which type of rational approximation should be used? The current types
  are "brasil", "chebfun" or "chebfunLB".

- shared_lib:

  Which shared lib to use for the cgeneric implementation? If "INLA", it
  will use the shared lib from INLA's installation. If 'rSPDE', then it
  will use the local installation (does not work if your installation is
  from CRAN). Otherwise, you can directly supply the path of the .so (or
  .dll) file.

- prior.kappa:

  a `list` containing the elements `meanlog` and `sdlog`, that is, the
  mean and standard deviation on the log scale.

- prior.tau:

  a list containing the elements `meanlog` and `sdlog`, that is, the
  mean and standard deviation on the log scale.

- prior.range:

  a `list` containing the elements `meanlog` and `sdlog`, that is, the
  mean and standard deviation on the log scale. Will not be used if
  prior.kappa is non-null.

- prior.std.dev:

  a `list` containing the elements `meanlog` and `sdlog`, that is, the
  mean and standard deviation on the log scale. Will not be used if
  prior.tau is non-null.

## Value

An INLA model.
