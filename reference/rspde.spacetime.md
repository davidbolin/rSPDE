# Space-Time Random Fields via SPDE Approximation

`rspde.spacetime` computes a Finite Element Method (FEM) approximation
of a Gaussian random field defined as the solution to the stochastic
partial differential equation (SPDE): \$\$d u + \gamma(\kappa^2 +
\kappa^{d/2}\rho\cdot\nabla - \Delta)^\alpha u = \sigma dW_C\$\$ where
\\C\\ is a Whittle-Matérn covariance operator with smoothness parameter
\\\beta\\ and range parameter \\\kappa\\. This function is designed to
handle space-time random fields using either 1D spatial models or
higher-dimensional FEM-based approaches.

## Usage

``` r
rspde.spacetime(
  mesh_space = NULL,
  mesh_time = NULL,
  space_loc = NULL,
  time_loc = NULL,
  drift = TRUE,
  alpha,
  beta,
  prior.kappa = NULL,
  prior.sigma = NULL,
  prior.rho = NULL,
  prior.gamma = NULL,
  prior.precision = NULL,
  graph_dirichlet = TRUE,
  bounded_rho = TRUE,
  bound_rho = NULL,
  shared_lib = "detect",
  debug = FALSE,
  ...
)
```

## Arguments

- mesh_space:

  Spatial mesh for the FEM approximation, or a `metric_graph` object for
  handling models on metric graphs.

- mesh_time:

  Temporal mesh for the FEM approximation.

- space_loc:

  A vector of spatial locations for mesh nodes in 1D spatial models.
  This should be provided when `mesh_space` is not specified.

- time_loc:

  A vector of temporal locations for mesh nodes. This should be provided
  when `mesh_time` is not specified.

- drift:

  Logical value indicating whether the drift term should be included. If
  `FALSE`, the drift coefficient \\\rho\\ is set to zero.

- alpha:

  Integer smoothness parameter \\\alpha\\.

- beta:

  Integer smoothness parameter \\\beta\\.

- prior.kappa:

  A list specifying the prior for the range parameter \\\kappa\\. This
  list may contain two elements: `mean` and/or `precision`, both of
  which must be numeric scalars (numeric vectors of length 1). The
  precision refers to the prior on \\\log(\kappa)\\. If `NULL`, default
  values will be used. The `mean` value is also used as starting value
  for kappa.

- prior.sigma:

  A list specifying the prior for the variance parameter \\\sigma\\.
  This list may contain two elements: `mean` and/or `precision`, both of
  which must be numeric scalars. The precision refers to the prior on
  \\\log(\sigma)\\. If `NULL`, default values will be used. The `mean`
  value is also used as starting value for sigma.

- prior.rho:

  A list specifying the prior for the drift coefficient \\\rho\\. This
  list may contain two elements: `mean` and/or `precision`, both of
  which must be numeric scalars if dimension is one, and numeric vectors
  of length 2 if dimension is 2. The precision applies directly to
  \\\rho\\ without log transformation. If `NULL`, default values will be
  used. Will not be used if `drift = FALSE`. The `mean` value is also
  used as starting value for rho.

- prior.gamma:

  A list specifying the prior for the weight \\\gamma\\ in the SPDE
  operator. This list may contain two elements: `mean` and/or
  `precision`, both of which must be numeric scalars. The precision
  refers to the prior on \\\log(\gamma)\\. If `NULL`, default values
  will be used. The `mean` value is also used as starting value for
  gamma.

- prior.precision:

  A precision matrix for \\\log(\kappa), \log(\sigma), \log(\gamma),
  \rho\\. This matrix replaces the precision element from `prior.kappa`,
  `prior.sigma`, `prior.gamma`, and `prior.rho` respectively. For
  dimension 1 `prior.precision` must be a 4x4 matrix. For dimension 2,
  \\\rho\\ is a vector of length 2, so in this case `prior.precision`
  must be a 5x5 matrix. If `NULL`, a diagonal precision matrix with
  default values will be used.

- graph_dirichlet:

  For models on metric graphs, use Dirichlet vertex conditions at
  vertices of degree 1?

- bounded_rho:

  Logical. Should `rho` be bounded to ensure the existence, uniqueness,
  and well-posedness of the solution? Defaults to `TRUE`. Note that this
  bounding is not a strict condition; there may exist values of rho
  beyond the upper bound that still satisfy these properties. When
  `bounded_rho = TRUE`, the `rspde_lme` models enforce bounded `rho` for
  consistency. If the estimated value of `rho` approaches the upper
  bound too closely, we recommend refitting the model with
  `bounded_rho = FALSE`. However, this should be done with caution, as
  it may lead to instability in some cases, though it can also result in
  a better model fit. The actual bound used for `rho` can be accessed
  from the `bound_rho` element of the returned object.

- bound_rho:

  A positive number specifying the bound for `rho`. If `NULL`, the
  default bound will be used.

- shared_lib:

  String specifying which shared library to use for the Cgeneric
  implementation. Options are "detect", "INLA", or "rSPDE". You may also
  specify the direct path to a .so (or .dll) file.

- debug:

  Logical value indicating whether to enable INLA debug mode.

- ...:

  Additional arguments passed internally for configuration purposes.

## Value

An object of class `inla_rspde_spacetime` representing the FEM
approximation of the space-time Gaussian random field.

## Examples

``` r
library(INLA)
library(MetricGraph)
#> This is MetricGraph 1.6.0
#> - See https://davidbolin.github.io/MetricGraph for vignettes and manuals.
#> 
#> Attaching package: ‘MetricGraph’
#> The following object is masked from ‘package:rSPDE’:
#> 
#>     cross_validation
#> The following object is masked from ‘package:stats’:
#> 
#>     filter
graph <- metric_graph$new()
#> Starting graph creation...
#> LongLat is set to FALSE
#> The current tolerances are:
#>   Vertex-Vertex 0.001
#>   Vertex-Edge 0.001
#>   Edge-Edge 0
#> Creating edges...
#> Setting edge weights...
#> Computing bounding box...
#> Setup edges and merge close vertices
#> Computing the relative positions of the edges...
#> Snap vertices to close edges
#> Total construction time: 1.55 secs
#> Creating and updating vertices...
#> Storing the initial graph...
#> Computing the relative positions of the edges...
graph$build_mesh(h = 0.1)
graph$compute_fem()

# Define the time locations
time_loc <- seq(from = 0, to = 10, length.out = 11)

# Create the model
model <- rspde.spacetime(mesh_space = graph,
                         time_loc = time_loc,
                         alpha = 2,
                         beta = 1)
```
