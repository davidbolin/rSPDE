# Summary for posteriors of field parameters for an `inla_rspde` model from a `rspde_result` object

Summary for posteriors of rSPDE field parameters in their original
scales.

## Usage

``` r
# S3 method for class 'rspde_result'
summary(object, digits = 6, ...)
```

## Arguments

- object:

  A `rspde_result` object.

- digits:

  integer, used for number formatting with signif()

- ...:

  Currently not used.

## Value

Returns a `data.frame` containing the summary.

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
