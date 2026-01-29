# Intrinsic Matern rSPDE model object for INLA

Creates an INLA object for a stationary intrinsic Matern model.
Currently, alpha is fixed to 2 and beta is fixed to 1.

## Usage

``` r
rspde.intrinsic.matern(
  mesh,
  alpha = 2,
  mean.correction = FALSE,
  prior.lkappa.mean = NULL,
  prior.ltau.mean = 1,
  prior.lkappa.prec = 0.1,
  prior.ltau.prec = 0.1,
  start.ltau = NULL,
  start.lkappa = NULL,
  true.scaling = TRUE,
  diagonal = 0,
  debug = FALSE,
  shared_lib = "detect",
  ...
)
```

## Arguments

- mesh:

  The mesh to build the model. It can be an `inla.mesh` or an
  `inla.mesh.1d` object. Otherwise, should be a list containing elements
  d, the dimension, C, the mass matrix, and G, the stiffness matrix.

- alpha:

  Smoothness parameter, need to be 1 or 2.

- mean.correction:

  Add mean correction for extreme value models?

- prior.lkappa.mean:

  Prior on log kappa to be used for the priors and for the starting
  values.

- prior.ltau.mean:

  Prior on log tau to be used for the priors and for the starting
  values.

- prior.lkappa.prec:

  Precision to be used on the prior on log kappa to be used for the
  priors and for the starting values.

- prior.ltau.prec:

  Precision to be used on the prior on log tau to be used for the priors
  and for the starting values.

- start.ltau:

  Starting value for log of tau.

- start.lkappa:

  Starting value for log of kappa.

- true.scaling:

  Compute the true normalizing constant manually? Default `TRUE`. The
  alternative is to set this to `FALSE` and set the `diagonal` argument
  to some small positive value. In the latter case, the model is
  approximated by a non-intrinsic model with a precision matrix that has
  the `diagonal` value added to the diagonal.

- diagonal:

  Value of diagonal correction for INLA stability. Default 0.

- debug:

  INLA debug argument

- shared_lib:

  Which shared lib to use for the cgeneric implementation? If "detect",
  it will check if the shared lib exists locally, in which case it will
  use it. Otherwise it will use INLA's shared library. If "INLA", it
  will use the shared lib from INLA's installation. If 'rSPDE', then it
  will use the local installation (does not work if your installation is
  from CRAN). Otherwise, you can directly supply the path of the .so (or
  .dll) file.

- ...:

  Only being used internally.

## Value

An INLA model.
