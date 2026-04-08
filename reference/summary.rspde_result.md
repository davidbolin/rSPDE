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
# devel version
if (requireNamespace("INLA", quietly = TRUE)) {
  library(INLA)

  set.seed(123)

  m <- 100
  loc_2d_mesh <- matrix(runif(m * 2), m, 2)
  mesh_2d <- inla.mesh.2d(
    loc = loc_2d_mesh,
    cutoff = 0.05,
    max.edge = c(0.1, 0.5)
  )
  sigma <- 1
  range <- 0.2
  nu <- 0.8
  kappa <- sqrt(8 * nu) / range
  op <- matern.operators(
    mesh = mesh_2d, nu = nu,
    range = range, sigma = sigma, m = 2,
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
    nu.upper.bound = 2
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
#>             mean        sd 0.025quant   0.5quant 0.975quant       mode
#> tau    0.0332515 0.0157152  0.0125518  0.0299566  0.0730695  0.0243231
#> kappa 16.5447000 3.2988400 10.9905000 16.2312000 23.8983000 15.6021000
#> nu     0.8724800 0.1225610  0.6366930  0.8713660  1.1154300  0.8701090
# devel.tag
# }
```
