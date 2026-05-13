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
#> Warning: `inla.mesh.2d()` was deprecated in INLA 23.08.18.
#> ℹ Please use `fmesher::fm_mesh_2d_inla()` instead.
#> ℹ For more information, see
#>   https://inlabru-org.github.io/fmesher/articles/inla_conversion.html
#> ℹ To silence these deprecation messages in old legacy code, set
#>   `inla.setOption(fmesher.evolution.warn = FALSE)`.
#> ℹ To ensure visibility of these messages in package tests, also set
#>   `inla.setOption(fmesher.evolution.verbosity = 'warn')`.
#> inlabru version: 2.14.1 
#> INLA version: 26.05.10 
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
#>     Pre = 0.365, Running = 1.39, Post = 0.424, Total = 2.18 
#> Fixed effects:
#>            mean    sd 0.025quant 0.5quant 0.975quant  mode kld
#> Intercept 0.306 0.169      -0.03    0.305      0.643 0.305   0
#> 
#> Random effects:
#>   Name     Model
#>     field CGeneric
#> 
#> Model hyperparameters:
#>                                           mean     sd 0.025quant 0.5quant
#> Precision for the Gaussian observations 150.25 84.062     36.062   133.03
#> Theta1 for field                         -8.54  4.525    -19.272    -7.79
#> Theta2 for field                          3.29  0.463      2.621     3.23
#> Theta3 for field                          2.43  2.444     -0.872     2.02
#>                                         0.975quant  mode
#> Precision for the Gaussian observations     354.33 96.34
#> Theta1 for field                             -2.43 -3.88
#> Theta2 for field                              4.37  2.89
#> Theta3 for field                              8.23 -0.08
#> 
#> Marginal log-Likelihood:  -106.76 
#>  is computed 
#> Posterior summaries for the linear predictor and the fitted values are computed
#> (Posterior marginals needs also 'control.compute=list(return.marginals.predictor=TRUE)')
#> 
# devel.tag
# }
```
