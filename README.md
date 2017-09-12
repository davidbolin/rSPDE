# Description #
rSPDE is an R package used for approximating fractional elliptic SPDEs 

L^\beta u(s) = W,

by rational SPDEs 

L_1u(s) = L_2 W 

where L_1 and L_2 are non-fractional operators. These rational SPDE can then be used for computatially efficient statistical inference.

Basic statistical operations such as likelihood evaluations and kriging predictions using the fractional approximations are also implemented.

For illustration purposes, the package contains a simple FEM implementation for models on R. For spatial models, the FEM implementation in the R-INLA package is recommended.

# Installation instructions #
The package can be installed directly from CRAN, or by using the command
```
#!r

devtools::install_bitbucket("rSPDE","davidbolin",ref="default")
```
in R. 
