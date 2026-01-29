# Finite element calculations for problems in 2D

This function computes mass and stiffness matrices for a mesh in 2D,
assuming Neumann boundary conditions.

## Usage

``` r
rSPDE.fem2d(FV, P)
```

## Arguments

- FV:

  Matrix where each row defines a triangle

- P:

  Locations of the nodes in the mesh.

## Value

The function returns a list with the following elements

- G :

  The stiffness matrix with elements \\(\nabla \phi_i, \nabla \phi_j)\\.

- C :

  The mass matrix with elements \\(\phi_i, \phi_j)\\.

- Cd :

  The mass lumped matrix with diagonal elements \\(\phi_i, 1)\\.

- Hxx :

  Matrix with elements \\(\partial_x \phi_i, \partial_x \phi_j)\\.

- Hyy :

  Matrix with elements \\(\partial_y \phi_i, \partial_y \phi_j)\\.

- Hxy :

  Matrix with elements \\(\partial_x \phi_i, \partial_y \phi_j)\\.

- Hyx :

  Matrix with elements \\(\partial_y \phi_i, \partial_x \phi_j)\\.

- Bx :

  Matrix with elements \\(\partial_x \phi_i, \phi_j)\\.

- By :

  Matrix with elements \\(\partial_y \phi_i, \phi_j)\\.

## See also

[`rSPDE.fem1d()`](https://davidbolin.github.io/rSPDE/reference/rSPDE.fem1d.md)

## Author

David Bolin <davidbolin@gmail.com>

## Examples

``` r
P <- rbind(c(0, 0), c(1, 0), c(1, 1), c(0, 1))
FV <- rbind(c(1, 2, 3), c(2, 3, 4))
fem <- rSPDE.fem2d(FV, P)
```
