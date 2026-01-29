# Object-based log-likelihood function for latent Gaussian fractional SPDE model

This function evaluates the log-likelihood function for a fractional
SPDE model \\L^\beta u(s) = W\\ that is observed under Gaussian
measurement noise: \\Y_i = u(s_i) + \epsilon_i\\, where \\\epsilon_i\\
are iid mean-zero Gaussian variables and \\x(s) = \mu(s) + u(s)\\, where
\\\mu(s)\\ is the expectation vector of the latent field.

## Usage

``` r
rSPDE.loglike(obj, Y, A, sigma.e, mu = 0)
```

## Arguments

- obj:

  The rational SPDE approximation, computed using
  [`fractional.operators()`](https://davidbolin.github.io/rSPDE/reference/fractional.operators.md),
  [`matern.operators()`](https://davidbolin.github.io/rSPDE/reference/matern.operators.md),
  or
  [`spde.matern.operators()`](https://davidbolin.github.io/rSPDE/reference/spde.matern.operators.md).

- Y:

  The observations, either a vector or a matrix where the columns
  correspond to independent replicates of observations.

- A:

  An observation matrix that links the measurement location to the
  finite element basis.

- sigma.e:

  The standard deviation of the measurement noise.

- mu:

  Expectation vector of the latent field (default = 0).

## Value

The log-likelihood value.

## Note

This example below shows how the function can be used to evaluate the
likelihood of a latent Matern model.

## See also

[`spde.matern.loglike()`](https://davidbolin.github.io/rSPDE/reference/spde.matern.loglike.md)

## Examples

``` r
# Sample a Gaussian Matern process on R using a rational approximation
kappa <- 10
sigma <- 1
nu <- 0.8
sigma.e <- 0.3
range <- 0.2

# create mass and stiffness matrices for a FEM discretization
x <- seq(from = 0, to = 1, length.out = 101)
fem <- rSPDE.fem1d(x)

# compute rational approximation
op <- matern.operators(
  range = range, sigma = sigma, nu = nu,
  loc_mesh = x, d = 1,
  type = "operator", parameterization = "matern"
)

# Sample the model
u <- simulate(op)

# Create some data
obs.loc <- runif(n = 10, min = 0, max = 1)
A <- rSPDE.A1d(x, obs.loc)
Y <- as.vector(A %*% u + sigma.e * rnorm(10))

# compute log-likelihood of the data
lik1 <- rSPDE.loglike(op, Y, A, sigma.e)
cat(lik1)
#> -8.21687
```
