# rSPDE stationary inlabru mapper

rSPDE stationary inlabru mapper

## Usage

``` r
bru_get_mapper.inla_rspde_matern1d(model, ...)

ibm_n.bru_mapper_inla_rspde_matern1d(mapper, ...)

ibm_values.bru_mapper_inla_rspde_matern1d(mapper, ...)

ibm_jacobian.bru_mapper_inla_rspde_matern1d(mapper, input, ...)
```

## Arguments

- model:

  An `inla_rspde_matern1d` object for which to construct or extract a
  mapper

- ...:

  Arguments passed on to other methods

- mapper:

  A `bru_mapper_inla_rspde_matern1d` object

- input:

  The values for which to produce a mapping matrix
