# Observation matrix for finite element discretization on R

A finite element discretization on R can be written as \\u(s) = \sum_i^n
u_i \varphi_i(s)\\ where \\\varphi_i(s)\\ is a piecewise linear "hat
function" centered at location \\x_i\\. This function computes an
\\m\times n\\ matrix \\A\\ that links the basis function in the
expansion to specified locations \\s = (s_1,\ldots, s_m)\\ in the domain
through \\A_ij = \varphi_j(s_i)\\.

## Usage

``` r
rSPDE.A1d(x, loc)
```

## Arguments

- x:

  The locations of the nodes in the FEM discretization.

- loc:

  The locations \\(s_1,\ldots, s_m)\\

## Value

The sparse matrix `A`.

## See also

[`rSPDE.fem1d()`](https://davidbolin.github.io/rSPDE/reference/rSPDE.fem1d.md)

## Author

David Bolin <davidbolin@gmail.com>

## Examples

``` r
# create mass and stiffness matrices for a FEM discretization on [0,1]
x <- seq(from = 0, to = 1, length.out = 101)
fem <- rSPDE.fem1d(x)

# create the observation matrix for some locations in the domain
obs.loc <- runif(n = 10, min = 0, max = 1)
A <- rSPDE.A1d(x, obs.loc)
```
