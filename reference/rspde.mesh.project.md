# Calculate a lattice projection to/from an `inla.mesh` for rSPDE objects

Calculate a lattice projection to/from an `inla.mesh` for rSPDE objects

## Usage

``` r
rspde.mesh.project(...)

rspde.mesh.projector(
  mesh,
  nu = NULL,
  rspde.order = 1,
  loc = NULL,
  lattice = NULL,
  xlim = NULL,
  ylim = NULL,
  dims = c(100, 100),
  projection = NULL,
  ...
)

# S3 method for class 'inla.mesh'
rspde.mesh.project(
  mesh,
  loc = NULL,
  field = NULL,
  rspde.order = 1,
  nu = NULL,
  ...
)

# S3 method for class 'rspde.mesh.projector'
rspde.mesh.project(projector, field, ...)

# S3 method for class 'inla.mesh.1d'
rspde.mesh.project(mesh, loc, field = NULL, rspde.order = 1, nu = NULL, ...)
```

## Arguments

- ...:

  Additional parameters.

- mesh:

  An `inla.mesh` or `inla.mesh.1d` object.

- nu:

  The smoothness parameter. If `NULL`, it will be assumed that nu was
  estimated.

- rspde.order:

  The order of the rational approximation.

- loc:

  Projection locations. Can be a matrix or a SpatialPoints or a
  SpatialPointsDataFrame object.

- lattice:

  An `inla.mesh.lattice` object.

- xlim:

  X-axis limits for a lattice. For R2 meshes, defaults to covering the
  domain.

- ylim:

  Y-axis limits for a lattice. For R2 meshes, defaults to covering the
  domain.

- dims:

  Lattice dimensions.

- projection:

  One of c("default", "longlat", "longsinlat", "mollweide").

- field:

  Basis function weights, one per mesh basis function, describing the
  function to be evaluated at the projection locations.

- projector:

  A `rspde.mesh.projector` object.

## Value

A list with projection information for rspde.mesh.project. For
rspde.mesh.projector(mesh, ...), a rspde.mesh.projector object. For
rspde.mesh.project(projector, field, ...), a field projected from the
mesh onto the locations given by the projector object.

## Details

This function is built upon the inla.mesh.project and
inla.mesh.projector functions from INLA.
