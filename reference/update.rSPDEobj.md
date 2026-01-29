# Update parameters of rSPDEobj objects

Function to change the parameters of a rSPDEobj object

## Usage

``` r
# S3 method for class 'rSPDEobj'
update(
  object,
  nu = NULL,
  alpha = NULL,
  kappa = NULL,
  sigma = NULL,
  range = NULL,
  tau = NULL,
  theta = NULL,
  m = NULL,
  mesh = NULL,
  loc_mesh = NULL,
  graph = NULL,
  range_mesh = NULL,
  parameterization = NULL,
  check_stationarity = TRUE,
  ...
)
```

## Arguments

- object:

  The operator-based rational SPDE approximation, computed using
  [`matern.operators()`](https://davidbolin.github.io/rSPDE/reference/matern.operators.md)
  with `type="operator"`

- nu:

  If non-null, update the shape parameter of the covariance function.

- alpha:

  If non-null, update the fractional order.

- kappa:

  If non-null, update the range parameter of the covariance function.

- sigma:

  If non-null, update the standard deviation of the covariance function.

- range:

  If non-null, update the range parameter of the covariance function.

- tau:

  If non-null, update the parameter tau.

- theta:

  If non-null, update the parameter theta, that connects tau and kappa
  to the model matrices.

- m:

  If non-null, update the order of the rational approximation, which
  needs to be a positive integer.

- mesh:

  An optional inla mesh. Replaces `d`, `C` and `G`.

- loc_mesh:

  The mesh locations used to construct the matrices C and G. This option
  should be provided if one wants to use the
  [`rspde_lme()`](https://davidbolin.github.io/rSPDE/reference/rspde_lme.md)
  function and will not provide neither graph nor mesh. Only works for
  1d data. Does not work for metric graphs. For metric graphs you should
  supply the graph using the `graph` argument.

- graph:

  An optional `metric_graph` object. Replaces `d`, `C` and `G`.

- range_mesh:

  The range of the mesh. Will be used to provide starting values for the
  parameters. Will be used if `mesh` and `graph` are `NULL`, and if one
  of the parameters (kappa or tau for spde parameterization, or sigma or
  range for matern parameterization) are not provided.

- parameterization:

  If non-null, update the parameterization. Only works for stationary
  models.

- check_stationarity:

  Logical; if TRUE, automatically returns a stationary model when
  tau/kappa (or sigma/range) are constant. Set to FALSE to keep a
  non-stationary model even when parameters are constant.

- ...:

  Currently not used.

## Value

It returns an object of class "rSPDEobj. This object contains the same
quantities listed in the output of
[`matern.operators()`](https://davidbolin.github.io/rSPDE/reference/matern.operators.md).

## See also

[`simulate.rSPDEobj()`](https://davidbolin.github.io/rSPDE/reference/simulate.rSPDEobj.md),
[`matern.operators()`](https://davidbolin.github.io/rSPDE/reference/matern.operators.md)

## Examples

``` r
# Compute the operator-based rational approximation of a
# Gaussian process with a Matern covariance function on R
kappa <- 10
sigma <- 1
nu <- 0.8
range <- sqrt(8 * nu) / kappa

# create mass and stiffness matrices for a FEM discretization
x <- seq(from = 0, to = 1, length.out = 101)
fem <- rSPDE.fem1d(x)

# compute rational approximation of covariance function at 0.5
op <- matern.operators(
  loc_mesh = x, nu = nu,
  range = range, sigma = sigma, d = 1, m = 2, type = "operator",
  parameterization = "matern"
)
op
#> Type of approximation:  Matern approximation 
#> Parameterization:  matern 
#> Parameters of covariance function: range =  0.2529822 , sigma =  1 , nu =  0.8 
#> Order or rational approximation:  2 
#> Size of discrete operators:  101  x  101 
#> Stationary Model

# Update the range parameter of the model:
op <- update(op, kappa = 20)
op
#> Type of approximation:  Matern approximation 
#> Parameterization:  matern 
#> Parameters of covariance function: range =  0.2529822 , sigma =  1 , nu =  0.8 
#> Order or rational approximation:  2 
#> Size of discrete operators:  101  x  101 
#> Stationary Model
```
