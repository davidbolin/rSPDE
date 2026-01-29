# Parameter-based log-likelihood for a latent Gaussian Matern SPDE model using a rational SPDE approximation

This function evaluates the log-likelihood function for observations of
a Gaussian process defined as the solution to the SPDE \$\$(\kappa(s) -
\Delta)^\beta (\tau(s)u(s)) = W.\$\$

## Usage

``` r
spde.matern.loglike(
  object,
  Y,
  A,
  sigma.e,
  mu = 0,
  nu = NULL,
  kappa = NULL,
  tau = NULL,
  theta = NULL,
  m = NULL
)
```

## Arguments

- object:

  The rational SPDE approximation, computed using
  [`spde.matern.operators()`](https://davidbolin.github.io/rSPDE/reference/spde.matern.operators.md)

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

- kappa:

  If non-null, updates the range parameter.

- tau:

  If non-null, updates the parameter tau.

- theta:

  If non-null, updates the parameter theta (that connects tau and kappa
  to the model matrices in `object`).

- m:

  If non-null, update the order of the rational approximation, which
  needs to be a positive integer.

## Value

The log-likelihood value.

## Details

The observations are assumed to be generated as \\Y_i = u(s_i) +
\epsilon_i\\, where \\\epsilon_i\\ are iid mean-zero Gaussian variables.
The latent model is approximated using a rational approximation of the
fractional SPDE model.

## See also

[`rSPDE.loglike()`](https://davidbolin.github.io/rSPDE/reference/rSPDE.loglike.md).

## Examples

``` r
# this example illustrates how the function can be used for maximum
# likelihood estimation
# Sample a Gaussian Matern process on R using a rational approximation
sigma.e <- 0.1
n.rep <- 10
n.obs <- 100
n.x <- 51
# create mass and stiffness matrices for a FEM discretization
x <- seq(from = 0, to = 1, length.out = n.x)
fem <- rSPDE.fem1d(x)
tau <- rep(0.5, n.x)
nu <- 0.8
alpha <- nu + 1 / 2
kappa <- rep(1, n.x)
# compute rational approximation
op <- spde.matern.operators(
  kappa = kappa, tau = tau, alpha = alpha,
  parameterization = "spde", d = 1,
  loc_mesh = x
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
mlik <- function(theta) {
  return(-spde.matern.loglike(op, Y, A,
    sigma.e = exp(theta[4]),
    nu = exp(theta[3]),
    kappa = exp(theta[2]),
    tau = exp(theta[1])
  ))
}
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
#>                 tau     kappa       nu   sigma.e
#> Truth     0.5000000 1.0000000 0.800000 0.1000000
#> Estimates 0.5481522 0.7711699 2.050953 0.1015812
# }
```
