# Prediction of an anisotropic Whittle-Matern field

The function is used for computing kriging predictions based on data
\\Y_i = u(s_i) + \epsilon_i\\, where \\\epsilon\\ is mean-zero Gaussian
measurement noise and \\u(s)\\ is defined by a SPDE as described in
[`matern2d.operators()`](https://davidbolin.github.io/rSPDE/reference/matern2d.operators.md).

## Usage

``` r
# S3 method for class 'CBrSPDEobj2d'
predict(
  object,
  A,
  Aprd,
  Y,
  sigma.e,
  mu = 0,
  compute.variances = FALSE,
  posterior_samples = FALSE,
  n_samples = 100,
  only_latent = FALSE,
  ...
)
```

## Arguments

- object:

  The covariance-based rational SPDE approximation, computed using
  [`matern2d.operators()`](https://davidbolin.github.io/rSPDE/reference/matern2d.operators.md)

- A:

  A matrix linking the measurement locations to the basis of the FEM
  approximation of the latent model.

- Aprd:

  A matrix linking the prediction locations to the basis of the FEM
  approximation of the latent model.

- Y:

  A vector with the observed data, can also be a matrix where the
  columns are observations of independent replicates of \\u\\.

- sigma.e:

  The standard deviation of the Gaussian measurement noise. Put to zero
  if the model does not have measurement noise.

- mu:

  Expectation vector of the latent field (default = 0).

- compute.variances:

  Set to also TRUE to compute the kriging variances.

- posterior_samples:

  If `TRUE`, posterior samples will be returned.

- n_samples:

  Number of samples to be returned. Will only be used if `sampling` is
  `TRUE`.

- only_latent:

  Should the posterior samples be only given to the laten model?

- ...:

  further arguments passed to or from other methods.

## Value

A list with elements

- mean :

  The kriging predictor (the posterior mean of u\|Y).

- variance :

  The posterior variances (if computed).

## Examples

``` r
 library(fmesher)
 n_loc <- 2000
 loc_2d_mesh <- matrix(runif(n_loc * 2), n_loc, 2)
 mesh_2d <- fm_mesh_2d(loc = loc_2d_mesh, cutoff = 0.01, max.edge = c(0.1, 0.5))
 op <- matern2d.operators(hx = 0.08, hy = 0.08, hxy = 0.5, nu = 0.5, 
 sigma = 1, mesh = mesh_2d)
 u <- simulate(op)
 n.obs <- 2000
 obs.loc <- cbind(runif(n.obs),runif(n.obs))
 A <- fm_basis(mesh_2d,obs.loc)
 sigma.e <- 0.1
 Y <- as.vector(A%*%u + sigma.e*rnorm(n.obs))
 A <- make_A(op, obs.loc)
 proj <- fm_evaluator(mesh_2d, dims = c(100, 100),
             xlim = c(0,1), ylim = c(0,1))
 Aprd <- make_A(op, proj$lattice$loc)
 u.krig <- predict(op, A = A, Aprd = Aprd, Y = Y, sigma.e = sigma.e)
```
