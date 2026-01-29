# Prediction of a space-time SPDE

The function is used for computing kriging predictions based on data
\\Y_i = u(s_i,t_i) + \epsilon_i\\, where \\\epsilon\\ is mean-zero
Gaussian measurement noise and \\u(s,t)\\ is defined by a
spatio-temporal SPDE as described in
[`spacetime.operators()`](https://davidbolin.github.io/rSPDE/reference/spacetime.operators.md).

## Usage

``` r
# S3 method for class 'spacetimeobj'
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
  [`spacetime.operators()`](https://davidbolin.github.io/rSPDE/reference/spacetime.operators.md)

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
s <- seq(from = 0, to = 20, length.out = 101)
t <- seq(from = 0, to = 20, length.out = 31)

op_cov <- spacetime.operators(space_loc = s, time_loc = t,
                             kappa = 5, sigma = 10, alpha = 1,
                             beta = 2, rho = 1, gamma = 0.05)
# generate data
sigma.e <- 0.01
n.obs <- 500
obs.loc <- data.frame(x = max(s)*runif(n.obs), 
                     t = max(t)*runif(n.obs))
A <- rSPDE.Ast(space_loc = s, time_loc = t, obs.s = obs.loc$x, obs.t = obs.loc$t)
Aprd <- Diagonal(dim(A)[2])
x <- simulate(op_cov, nsim = 1) 
Y <- A%*%x + sigma.e*rnorm(n.obs)
u.krig <- predict(op_cov, A, Aprd, Y, sigma.e)
```
