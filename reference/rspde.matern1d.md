# Matern rSPDE model object for INLA

Creates an INLA object for a stationary Matern model with general
smoothness parameter.

## Usage

``` r
rspde.matern1d(
  loc,
  nu.upper.bound = NULL,
  rspde.order = 1,
  nu = NULL,
  parameterization = c("spde", "matern", "matern2"),
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
  prior.nu.dist = c("beta", "lognormal"),
  nu.prec.inc = 1,
  type.rational.approx = c("brasil", "chebfun", "chebfunLB"),
  debug = FALSE,
  shared_lib = "detect",
  ...
)
```

## Arguments

- loc:

  A vector of spatial locations.

- nu.upper.bound:

  Upper bound for the smoothness parameter. If `NULL`, it will be set to
  2.

- rspde.order:

  The order of the covariance-based rational SPDE approach. The default
  order is 1.

- nu:

  If nu is set to a parameter, nu will be kept fixed and will not be
  estimated. If nu is `NULL`, it will be estimated.

- parameterization:

  Which parameterization to use? `matern` uses range, std. deviation and
  nu (smoothness). `spde` uses kappa, tau and nu (smoothness). `matern2`
  uses range-like (1/kappa), variance and nu (smoothness). The default
  is `spde`.

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
  "beta" or "lognormal". The default is "lognormal".

- nu.prec.inc:

  Amount to increase the precision in the beta prior distribution.

- type.rational.approx:

  Which type of rational approximation should be used? The current types
  are "brasil", "chebfun" or "chebfunLB".

- debug:

  INLA debug argument

- shared_lib:

  Which shared lib to use for the cgeneric implementation? If "detect",
  it will check if the shared lib exists locally, in which case it will
  use it. Otherwise it will use INLA's shared library. If "INLA", it
  will use the shared lib from INLA's installation. If 'rSPDE', then it
  will use the local installation (does not work if your installation is
  from CRAN). Otherwise, you can directly supply the path of the .so (or
  .dll) file.

- ...:

  Only being used internally.

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
