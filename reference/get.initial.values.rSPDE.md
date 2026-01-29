# Initial values for log-likelihood optimization in rSPDE models with a latent stationary Gaussian Matern model

Auxiliar function to obtain domain-based initial values for
log-likelihood optimization in rSPDE models with a latent stationary
Gaussian Matern model

## Usage

``` r
get.initial.values.rSPDE(
  mesh = NULL,
  mesh.range = NULL,
  graph.obj = NULL,
  n.spde = 1,
  dim = NULL,
  B.tau = NULL,
  B.kappa = NULL,
  B.sigma = NULL,
  B.range = NULL,
  nu = NULL,
  parameterization = c("matern", "spde"),
  include.nu = TRUE,
  log.scale = TRUE,
  nu.upper.bound = NULL
)
```

## Arguments

- mesh:

  An in INLA mesh

- mesh.range:

  The range of the mesh.

- graph.obj:

  A `metric_graph` object. To be used in case both `mesh` and
  `mesh.range` are `NULL`.

- n.spde:

  The number of basis functions in the mesh model.

- dim:

  The dimension of the domain.

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

- nu:

  The smoothness parameter.

- parameterization:

  Which parameterization to use? `matern` uses range, std. deviation and
  nu (smoothness). `spde` uses kappa, tau and nu (smoothness). The
  default is `matern`.

- include.nu:

  Should we also provide an initial guess for nu?

- log.scale:

  Should the results be provided in log scale?

- nu.upper.bound:

  Should an upper bound for nu be considered?

## Value

A vector of the form (theta_1,theta_2,theta_3) or where theta_1 is the
initial guess for tau, theta_2 is the initial guess for kappa and
theta_3 is the initial guess for nu.
