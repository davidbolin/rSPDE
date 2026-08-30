# rSPDE model index vector generation

Generates a list of named index vectors for an rSPDE model.

## Usage

``` r
rspde.make.index(
  name,
  n.spde = NULL,
  n.group = 1,
  n.repl = 1,
  mesh = NULL,
  rspde.order = 1,
  nu = NULL,
  dim = NULL
)
```

## Arguments

- name:

  A character string with the base name of the effect.

- n.spde:

  The number of basis functions in the mesh model.

- n.group:

  The size of the group model.

- n.repl:

  The total number of replicates.

- mesh:

  An `inla.mesh`, an `inla.mesh.1d` object or a `metric_graph` object.

- rspde.order:

  The order of the rational approximation

- nu:

  If `NULL`, then the model will assume that nu will be estimated. If nu
  is fixed, you should provide the value of nu.

- dim:

  the dimension of the domain. Should only be provided if `mesh` is not
  provided.

## Value

A list of named index vectors.

- name:

  Indices into the vector of latent variables

- name.group:

  'group' indices

- name.repl:

  Indices for replicates

## Examples

``` r
# \donttest{
if (rspde_safe_inla()) {
  library(INLA)

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
  Abar <- rspde.make.A(mesh = mesh_2d, loc = loc_2d_mesh)
  mesh.index <- rspde.make.index(name = "field", mesh = mesh_2d)
  st.dat <- inla.stack(
    data = list(y = as.vector(y)),
    A = Abar,
    effects = mesh.index
  )
}
# }
```
