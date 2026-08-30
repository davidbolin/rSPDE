# Observation/prediction matrices for rSPDE models with integer smoothness.

Constructs observation/prediction weight matrices for rSPDE models with
integer smoothness based on `inla.mesh` or `inla.mesh.1d` objects.

## Usage

``` r
spde.make.A(
  mesh = NULL,
  loc = NULL,
  A = NULL,
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
if (requireNamespace("fmesher", quietly = TRUE)) {
  library(fmesher)

  set.seed(123)
  loc <- matrix(runif(100 * 2) * 100, 100, 2)
  mesh <- fm_mesh_2d(
    loc = loc,
    cutoff = 50,
    max.edge = c(50, 500)
  )
  A <- spde.make.A(mesh, loc = loc)
}
# }
```
