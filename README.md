# Description #
rSPDE is an R package used for approximating fractional SPDEs 

L^\beta u(s) = W,

by rational SPDEs 

P_l u(s) = P_r W 

where P_l and P_r are non-fractional operators. These rational SPDE can then be used for computatially efficient statistical inference.

Basic statistical operations such as likelihood evaluations and kriging predictions using the fractional approximations are also implemented.

For illustration purposes, the package contains a simple FEM implementation for models on R. For spatial models, the FEM implementation in the R-INLA package is recommended.

# Reference #
D. Bolin and K. Kichner, [The rational SPDE approach for Gaussian random fields with general smoothness][ref]. Preprint, arXiv:1711.04333

# Installation instructions #
The package can be installed in R using the command
```
#!r

devtools::install_bitbucket("davidbolin/rspde",ref="default")
```

[ref]: https://arxiv.org/abs/1711.04333  "The SPDE approach for Gaussian random fields with general smoothness"