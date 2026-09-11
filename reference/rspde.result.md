# rSPDE result extraction from INLA estimation results

Extract field and parameter values and distributions for an rspde effect
from an inla result object.

## Usage

``` r
rspde.result(
  inla,
  name,
  rspde,
  compute.summary = TRUE,
  parameterization = "detect",
  n_samples = 5000,
  n_density = 1024
)
```

## Arguments

- inla:

  An `inla` object obtained from a call to
  [`inla()`](https://rdrr.io/pkg/INLA/man/inla.html).

- name:

  A character string with the name of the rSPDE effect in the inla
  formula.

- rspde:

  The `inla_rspde` object used for the effect in the inla formula.

- compute.summary:

  Should the summary be computed?

- parameterization:

  If 'detect', the parameterization from the model will be used.
  Otherwise, the options are 'spde', 'matern' and 'matern2'.

- n_samples:

  The number of samples to be used if parameterization is different from
  the one used to fit the model.

- n_density:

  The number of equally spaced points to estimate the density.

## Value

If the model was fitted with `matern` parameterization (the default), it
returns a list containing:

- marginals.range:

  Marginal densities for the range parameter

- marginals.log.range:

  Marginal densities for log(range)

- marginals.std.dev:

  Marginal densities for std. deviation

- marginals.log.std.dev:

  Marginal densities for log(std. deviation)

- marginals.values:

  Marginal densities for the field values

- summary.log.range:

  Summary statistics for log(range)

- summary.log.std.dev:

  Summary statistics for log(std. deviation)

- summary.values:

  Summary statistics for the field values

If `compute.summary` is `TRUE`, then the list will also contain

- summary.kappa:

  Summary statistics for kappa

- summary.tau:

  Summary statistics for tau

If the model was fitted with the `spde` parameterization, it returns a
list containing:

- marginals.kappa:

  Marginal densities for kappa

- marginals.log.kappa:

  Marginal densities for log(kappa)

- marginals.log.tau:

  Marginal densities for log(tau)

- marginals.tau:

  Marginal densities for tau

- marginals.values:

  Marginal densities for the field values

- summary.log.kappa:

  Summary statistics for log(kappa)

- summary.log.tau:

  Summary statistics for log(tau)

- summary.values:

  Summary statistics for the field values

If `compute.summary` is `TRUE`, then the list will also contain

- summary.kappa:

  Summary statistics for kappa

- summary.tau:

  Summary statistics for tau

For both cases, if nu was estimated, then the list will also contain

- marginals.nu:

  Marginal densities for nu

If nu was estimated and a beta prior was used, then the list will also
contain

- marginals.logit.nu:

  Marginal densities for logit(nu)

- summary.logit.nu:

  Marginal densities for logit(nu)

If nu was estimated and a truncated lognormal prior was used, then the
list will also contain

- marginals.log.nu:

  Marginal densities for log(nu)

- summary.log.nu:

  Marginal densities for log(nu)

If nu was estimated and `compute.summary` is `TRUE`, then the list will
also contain

- summary.nu:

  Summary statistics for nu

## Examples

``` r
# \donttest{
if (rspde_safe_inla()) {
  library(INLA)

  set.seed(123)

  m <- 50
  loc_2d_mesh <- matrix(runif(m * 2), m, 2)
  mesh_2d <- inla.mesh.2d(
    loc = loc_2d_mesh,
    cutoff = 0.1,
    max.edge = c(0.5, 0.5)
  )
  sigma <- 1
  range <- 0.2
  nu <- 0.8
  kappa <- sqrt(8 * nu) / range
  op <- matern.operators(
    mesh = mesh_2d, nu = nu,
    range = range, sigma = sigma, m = 1,
    parameterization = "matern"
  )
  u <- simulate(op)
  A <- inla.spde.make.A(
    mesh = mesh_2d,
    loc = loc_2d_mesh
  )
  sigma.e <- 0.1
  y <- A %*% u + rnorm(m) * sigma.e
  Abar <- rspde.make.A(mesh = mesh_2d, loc = loc_2d_mesh)
  mesh.index <- rspde.make.index(name = "field", mesh = mesh_2d)
  st.dat <- inla.stack(
    data = list(y = as.vector(y)),
    A = Abar,
    effects = mesh.index
  )
  rspde_model <- rspde.matern(
    mesh = mesh_2d,
    nu.upper.bound = 1
  )
  f <- y ~ -1 + f(field, model = rspde_model)
  rspde_fit <- inla(f,
    data = inla.stack.data(st.dat),
    family = "gaussian",
    control.predictor =
      list(A = inla.stack.A(st.dat))
  )
  result <- rspde.result(rspde_fit, "field", rspde_model)
  summary(result)
}
#> Warning: the mean or mode of nu is very close to nu.upper.bound, please consider increasing nu.upper.bound, and refitting the model.
#>            mean        sd 0.025quant   0.5quant 0.975quant        mode
#> tau    0.146132  0.309758 0.00009511  0.0341516   0.963156 1.06359e-06
#> kappa 14.892600 13.024600 4.01405000 10.7904000  50.476100 6.69931e+00
#> nu     0.736905  0.294736 0.09246090  0.8774770   0.999912 9.99990e-01
# }
```
