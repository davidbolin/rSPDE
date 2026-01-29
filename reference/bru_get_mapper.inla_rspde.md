# rSPDE inlabru mapper

rSPDE inlabru mapper

## Usage

``` r
ibm_n.bru_mapper_inla_rspde_fintrinsic(mapper, ...)

ibm_values.bru_mapper_inla_rspde_fintrinsic(mapper, ...)

ibm_jacobian.bru_mapper_inla_rspde_fintrinsic(mapper, input, ...)

bru_get_mapper.inla_rspde(model, ...)

ibm_n.bru_mapper_inla_rspde(mapper, ...)

ibm_values.bru_mapper_inla_rspde(mapper, ...)

ibm_jacobian.bru_mapper_inla_rspde(mapper, input, ...)
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
#> inlabru version: 2.13.0
#> INLA version: 26.01.26-1
#> Components:
#> Intercept: main = linear(1), group = exchangeable(1L), replicate = iid(1L), NULL
#> field: main = cgeneric(cbind(x1, x2)), group = exchangeable(1L), replicate = iid(1L), NULL
#> Observation models:
#>   Family: 'gaussian'
#>     Tag: <No tag>
#>     Data class: 'data.frame'
#>     Response class: 'numeric'
#>     Predictor: y ~ .
#>     Additive/Linear: TRUE/TRUE
#>     Used components: effects[Intercept, field], latent[]
#> Time used:
#>     Pre = 0.563, Running = 1.35, Post = 0.107, Total = 2.02 
#> Fixed effects:
#>            mean    sd 0.025quant 0.5quant 0.975quant  mode kld
#> Intercept 0.305 0.169     -0.028    0.305      0.641 0.305   0
#> 
#> Random effects:
#>   Name     Model
#>     field CGeneric
#> 
#> Model hyperparameters:
#>                                           mean     sd 0.025quant 0.5quant
#> Precision for the Gaussian observations 155.13 76.326     49.368  140.761
#> Theta1 for field                         -3.78  0.498     -4.567   -3.833
#> Theta2 for field                          2.94  0.202      2.540    2.943
#> Theta3 for field                         -0.12  0.268     -0.721   -0.093
#>                                         0.975quant    mode
#> Precision for the Gaussian observations    342.064 113.162
#> Theta1 for field                            -2.652  -4.094
#> Theta2 for field                             3.336   2.947
#> Theta3 for field                             0.314   0.037
#> 
#> Marginal log-Likelihood:  -108.72 
#>  is computed 
#> Posterior summaries for the linear predictor and the fitted values are computed
#> (Posterior marginals needs also 'control.compute=list(return.marginals.predictor=TRUE)')
#> 
# devel.tag
# }
```
