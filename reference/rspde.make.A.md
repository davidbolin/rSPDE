# Observation/prediction matrices for rSPDE models.

Constructs observation/prediction weight matrices for rSPDE models based
on `inla.mesh` or `inla.mesh.1d` objects.

## Usage

``` r
rspde.make.A(
  mesh = NULL,
  loc = NULL,
  A = NULL,
  dim = NULL,
  rspde.order = 1,
  nu = NULL,
  index = NULL,
  group = NULL,
  repl = 1L,
  n.group = NULL,
  n.repl = NULL
)
```

## Arguments

- mesh:

  An `inla.mesh`, an `inla.mesh.1d` object or a `metric_graph` object.

- loc:

  Locations, needed if an INLA mesh is provided

- A:

  The A matrix from the standard SPDE approach, such as the matrix
  returned by `inla.spde.make.A`. Should only be provided if `mesh` is
  not provided.

- dim:

  the dimension. Should only be provided if an `mesh` is not provided.

- rspde.order:

  The order of the covariance-based rational SPDE approach.

- nu:

  If `NULL`, then the model will assume that nu will be estimated. If nu
  is fixed, you should provide the value of nu.

- index:

  For each observation/prediction value, an index into loc. Default is
  seq_len(nrow(A.loc)).

- group:

  For each observation/prediction value, an index into the group model.

- repl:

  For each observation/prediction value, the replicate index.

- n.group:

  The size of the group model.

- n.repl:

  The total number of replicates.

## Value

The \\A\\ matrix for rSPDE models.

## Examples

``` r
# \donttest{
# devel version
if (requireNamespace("INLA", quietly = TRUE)) {
  library(INLA)

  set.seed(123)
  loc <- matrix(runif(100 * 2) * 100, 100, 2)
  mesh <- inla.mesh.2d(
    loc = loc,
    cutoff = 50,
    max.edge = c(50, 500)
  )
  A <- rspde.make.A(mesh, loc = loc, rspde.order = 3)
}
# devel.tag
# }
```
