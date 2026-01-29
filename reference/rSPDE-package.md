# Rational approximations of fractional SPDEs.

`rSPDE` is used for approximating fractional elliptic SPDEs \$\$L^\beta
(\tau u(s)) = W,\$\$ where \\L\\ is a differential operator and
\\\beta\>0\\ is a general fractional power.

## Details

The approximation is based on a rational approximation of the fractional
operator, and allows for computationally efficient inference and
simulation.

The main functions for computing rational approximation objects are:

- [`fractional.operators()`](https://davidbolin.github.io/rSPDE/reference/fractional.operators.md):

  works for general rational operators

- [`matern.operators()`](https://davidbolin.github.io/rSPDE/reference/matern.operators.md):

  works for random fields with stationary Matern covariance functions

- [`spde.matern.operators()`](https://davidbolin.github.io/rSPDE/reference/spde.matern.operators.md):

  works for random fields with defined as solutions to a possibly
  non-stationary Matern-type SPDE model.

- [`rspde.matern()`](https://davidbolin.github.io/rSPDE/reference/rspde.matern.md):

  R-INLA implementation of the covariance-based rational approximation
  for random fields with stationary Matern covariance functions

Basic statistical operations such as likelihood evaluations (see
`[rSPDE.loglike], [rSPDE.matern.loglike]`) and kriging predictions (see
`[predict.rSPDEobj], [predict.CBrSPDEobj]`) using the rational
approximations are also implemented.

For illustration purposes, the package contains a simple FEM
implementation for models on R. For spatial models, the FEM
implementation in the `R-INLA` package is recommended.

For a more detailed introduction to the package, see the rSPDE
Vignettes.

## See also

Useful links:

- <https://davidbolin.github.io/rSPDE/>

- Report bugs at <https://github.com/davidbolin/rSPDE/issues>

## Author

**Maintainer**: David Bolin <davidbolin@gmail.com>

Authors:

- Alexandre Simas <alexandre.impa@gmail.com>

Other contributors:

- Finn Lindgren <finn.lindgren@ed.ac.uk> \[contributor\]
