# rSPDE stationary inlabru mapper

rSPDE stationary inlabru mapper

## Usage

``` r
# S3 method for class 'inla_rspde_matern1d'
bru_get_mapper(model, ...)

# S3 method for class 'bru_mapper_inla_rspde_matern1d'
ibm_n(mapper, ...)

# S3 method for class 'bru_mapper_inla_rspde_matern1d'
ibm_values(mapper, ...)

# S3 method for class 'bru_mapper_inla_rspde_matern1d'
ibm_jacobian(mapper, input, ...)
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
