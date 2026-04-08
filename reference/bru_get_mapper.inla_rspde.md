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
#> INLA version: 26.03.19 
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
#>     Pre = 0.326, Running = 1.41, Post = 0.102, Total = 1.84 
#> Fixed effects:
#>            mean    sd 0.025quant 0.5quant 0.975quant  mode kld
#> Intercept 0.305 0.169      -0.03    0.305      0.642 0.305   0
#> 
#> Random effects:
#>   Name     Model
#>     field CGeneric
#> 
#> Model hyperparameters:
#>                                            mean     sd 0.025quant 0.5quant
#> Precision for the Gaussian observations 173.230 86.584     56.928  156.033
#> Theta1 for field                         -3.729  0.542     -4.546   -3.799
#> Theta2 for field                          2.941  0.202      2.538    2.943
#> Theta3 for field                         -0.131  0.279     -0.768   -0.099
#>                                         0.975quant    mode
#> Precision for the Gaussian observations    388.718 124.749
#> Theta1 for field                            -2.481  -4.150
#> Theta2 for field                             3.332   2.952
#> Theta3 for field                             0.305   0.059
#> 
#> Marginal log-Likelihood:  -108.69 
#>  is computed 
#> Posterior summaries for the linear predictor and the fitted values are computed
#> (Posterior marginals needs also 'control.compute=list(return.marginals.predictor=TRUE)')
#> 
# devel.tag
# }
```
