# Finite element calculations for problems on R

This function computes mass and stiffness matrices for a FEM
approximation on R, assuming Neumann boundary conditions. These matrices
are needed when discretizing the operators in rational approximations.

## Usage

``` r
rSPDE.fem1d(x)
```

## Arguments

- x:

  Locations of the nodes in the FEM approximation.

## Value

The function returns a list with the following elements

- G :

  The stiffness matrix with elements \\(\nabla \phi_i, \nabla \phi_j)\\.

- C :

  The mass matrix with elements \\(\phi_i, \phi_j)\\.

- Cd :

  Mass lumped mass matrix.

- B :

  Matrix with elements \\(\nabla \phi_i, \phi_j)\\.

## See also

[`rSPDE.A1d()`](https://davidbolin.github.io/rSPDE/reference/rSPDE.A1d.md)

## Author

David Bolin <davidbolin@gmail.com>

## Examples

``` r
# create mass and stiffness matrices for a FEM discretization on [0,1]
x <- seq(from = 0, to = 1, length.out = 101)
fem <- rSPDE.fem1d(x)
```
