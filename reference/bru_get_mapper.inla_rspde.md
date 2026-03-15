# rSPDE inlabru mapper

rSPDE inlabru mapper

## Usage

``` r
# S3 method for class 'bru_mapper_inla_rspde_fintrinsic'
ibm_n(mapper, ...)

# S3 method for class 'bru_mapper_inla_rspde_fintrinsic'
ibm_values(mapper, ...)

# S3 method for class 'bru_mapper_inla_rspde_fintrinsic'
ibm_jacobian(mapper, input, ...)

# S3 method for class 'inla_rspde'
bru_get_mapper(model, ...)

# S3 method for class 'bru_mapper_inla_rspde'
ibm_n(mapper, ...)

# S3 method for class 'bru_mapper_inla_rspde'
ibm_values(mapper, ...)

# S3 method for class 'bru_mapper_inla_rspde'
ibm_jacobian(mapper, input, ...)
```

## Arguments

- mapper:

  A `bru_mapper_inla_rspde` object

- ...:

  Arguments passed on to other methods

- input:

  The values for which to produce a mapping matrix

- model:

  An `inla_rspde` object for which to construct or extract a mapper

## Examples

``` r
# \donttest{
# devel version
if (requireNamespace("INLA", quietly = TRUE) &&
  requireNamespace("inlabru", quietly = TRUE)) {
  library(INLA)
  library(inlabru)

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
  y <- as.vector(y)

  data_df <- data.frame(
    y = y, x1 = loc_2d_mesh[, 1],
    x2 = loc_2d_mesh[, 2]
  )
  rspde_model <- rspde.matern(
    mesh = mesh_2d,
    nu_upper_bound = 2
  )

  cmp <- y ~ Intercept(1) +
    field(cbind(x1,x2), model = rspde_model)


  rspde_fit <- bru(cmp, data = data_df)
  summary(rspde_fit)
}
#> 
#> Loading required package: fmesher
#> Warning: `inla.mesh.2d()` was deprecated in INLA 23.08.18.
#> ℹ Please use `fmesher::fm_mesh_2d_inla()` instead.
#> ℹ For more information, see
#>   https://inlabru-org.github.io/fmesher/articles/inla_conversion.html
#> ℹ To silence these deprecation messages in old legacy code, set
#>   `inla.setOption(fmesher.evolution.warn = FALSE)`.
#> ℹ To ensure visibility of these messages in package tests, also set
#>   `inla.setOption(fmesher.evolution.verbosity = 'warn')`.
#> inlabru version: 2.14.0 
#> INLA version: 26.02.06 
#> Latent components:
#> Intercept: main = linear(1)
#> field: main = cgeneric(cbind(x1, x2))
#> Observation models:
#>   Model tag: <No tag>
#>     Family: 'gaussian'
#>     Data class: 'data.frame'
#>     Response class: 'numeric'
#>     Predictor: y ~ Intercept + field
#>     Additive/Linear/Rowwise: TRUE/TRUE/TRUE
#>     Used components: effect[Intercept, field], latent[] 
#> Time used:
#>     Pre = 0.343, Running = 1.44, Post = 0.108, Total = 1.89 
#> Fixed effects:
#>            mean   sd 0.025quant 0.5quant 0.975quant  mode kld
#> Intercept 0.306 0.17     -0.031    0.306      0.645 0.306   0
#> 
#> Random effects:
#>   Name     Model
#>     field CGeneric
#> 
#> Model hyperparameters:
#>                                            mean      sd 0.025quant 0.5quant
#> Precision for the Gaussian observations 174.007 119.933     37.647  143.961
#> Theta1 for field                         -3.609   0.564     -4.751   -3.599
#> Theta2 for field                          2.870   0.212      2.454    2.870
#> Theta3 for field                         -0.222   0.316     -0.824   -0.228
#>                                         0.975quant   mode
#> Precision for the Gaussian observations    486.027 94.389
#> Theta1 for field                            -2.530 -3.553
#> Theta2 for field                             3.289  2.867
#> Theta3 for field                             0.418 -0.256
#> 
#> Marginal log-Likelihood:  -108.14 
#>  is computed 
#> Posterior summaries for the linear predictor and the fitted values are computed
#> (Posterior marginals needs also 'control.compute=list(return.marginals.predictor=TRUE)')
#> 
# devel.tag
# }
```
