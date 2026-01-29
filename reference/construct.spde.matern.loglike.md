# Constructor of Matern loglikelihood functions for non-stationary models.

This function evaluates the log-likelihood function for observations of
a non-stationary Gaussian process defined as the solution to the SPDE
\$\$(\kappa(s) - \Delta)^\beta (\tau(s)u(s)) = W.\$\$

The observations are assumed to be generated as \\Y_i = u(s_i) +
\epsilon_i\\, where \\\epsilon_i\\ are iid mean-zero Gaussian variables.
The latent model is approximated using a rational approximation of the
fractional SPDE model.

## Usage

``` r
construct.spde.matern.loglike(
  object,
  Y,
  A,
  sigma.e = NULL,
  mu = 0,
  nu = NULL,
  m = NULL,
  log_scale = TRUE,
  return_negative_likelihood = TRUE
)
```

## Arguments

- object:

  The rational SPDE approximation, computed using
  [`matern.operators()`](https://davidbolin.github.io/rSPDE/reference/matern.operators.md)

- Y:

  The observations, either a vector or a matrix where the columns
  correspond to independent replicates of observations.

- A:

  An observation matrix that links the measurement location to the
  finite element basis.

- sigma.e:

  IF non-null, the standard deviation of the measurement noise will be
  kept fixed in the returned likelihood.

- mu:

  Expectation vector of the latent field (default = 0).

- nu:

  If non-null, the shape parameter will be kept fixed in the returned
  likelihood.

- m:

  If non-null, update the order of the rational approximation, which
  needs to be a positive integer.

- log_scale:

  Should the parameters be evaluated in log-scale?

- return_negative_likelihood:

  Return minus the likelihood to turn the maximization into a
  minimization?

## Value

The log-likelihood function. The parameters of the returned function are
given in the order theta, nu, sigma.e, whenever they are available.

## See also

[`matern.operators()`](https://davidbolin.github.io/rSPDE/reference/matern.operators.md),
[`predict.CBrSPDEobj()`](https://davidbolin.github.io/rSPDE/reference/predict.CBrSPDEobj.md)

## Examples

``` r
# this example illustrates how the function can be used for maximum
# likelihood estimation
# Sample a Gaussian Matern process on R using a rational approximation
set.seed(123)
sigma.e <- 0.1
n.rep <- 10
n.obs <- 100
n.x <- 51
# create mass and stiffness matrices for a FEM discretization
x <- seq(from = 0, to = 1, length.out = n.x)
fem <- rSPDE.fem1d(x)
tau <- rep(0.5, n.x)
nu <- 0.8
alpha <- nu + 0.5
kappa <- rep(1, n.x)
# Matern parameterization
# compute rational approximation
op <- spde.matern.operators(
  loc_mesh = x,
  kappa = kappa, tau = tau, alpha = alpha,
  parameterization = "spde", d = 1
)
# Sample the model
u <- simulate(op, n.rep)
# Create some data
obs.loc <- runif(n = n.obs, min = 0, max = 1)
A <- rSPDE.A1d(x, obs.loc)
noise <- rnorm(n.obs * n.rep)
dim(noise) <- c(n.obs, n.rep)
Y <- as.matrix(A %*% u + sigma.e * noise)
# define negative likelihood function for optimization using matern.loglike
mlik <- construct.spde.matern.loglike(op, Y, A)
#' #The parameters can now be estimated by minimizing mlik with optim
# \donttest{
# Choose some reasonable starting values depending on the size of the domain
theta0 <- log(c(1 / sqrt(var(c(Y))), sqrt(8), 0.9, 0.01))
# run estimation and display the results
theta <- optim(theta0, mlik)
print(data.frame(
  tau = c(tau[1], exp(theta$par[1])), kappa = c(kappa[1], exp(theta$par[2])),
  nu = c(nu, exp(theta$par[3])), sigma.e = c(sigma.e, exp(theta$par[4])),
  row.names = c("Truth", "Estimates")
))
#>                 tau    kappa       nu    sigma.e
#> Truth     0.5000000 1.000000 0.800000 0.10000000
#> Estimates 0.3850559 1.059108 1.535358 0.09917289
# }

# SPDE parameterization
# compute rational approximation
op <- spde.matern.operators(
  kappa = kappa, tau = tau, alpha = alpha,
  loc_mesh = x, d = 1,
  parameterization = "spde"
)
# Sample the model
u <- simulate(op, n.rep)
# Create some data
obs.loc <- runif(n = n.obs, min = 0, max = 1)
A <- rSPDE.A1d(x, obs.loc)
noise <- rnorm(n.obs * n.rep)
dim(noise) <- c(n.obs, n.rep)
Y <- as.matrix(A %*% u + sigma.e * noise)
# define negative likelihood function for optimization using matern.loglike
mlik <- construct.spde.matern.loglike(op, Y, A)
#' #The parameters can now be estimated by minimizing mlik with optim
# \donttest{
# Choose some reasonable starting values depending on the size of the domain
theta0 <- log(c(1 / sqrt(var(c(Y))), sqrt(8), 0.9, 0.01))
# run estimation and display the results
theta <- optim(theta0, mlik)
print(data.frame(
  tau = c(tau[1], exp(theta$par[1])), kappa = c(kappa[1], exp(theta$par[2])),
  nu = c(nu, exp(theta$par[3])), sigma.e = c(sigma.e, exp(theta$par[4])),
  row.names = c("Truth", "Estimates")
))
#>                 tau    kappa       nu    sigma.e
#> Truth     0.5000000 1.000000 0.800000 0.10000000
#> Estimates 0.5511011 1.230306 1.251033 0.09867291
# }
```
