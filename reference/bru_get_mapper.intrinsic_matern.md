# rSPDE inlabru mapper

rSPDE inlabru mapper

## Usage

``` r
# S3 method for class 'intrinsic_matern'
bru_get_mapper(model, ...)

# S3 method for class 'bru_mapper_intrinsic_matern'
ibm_n(mapper, ...)

# S3 method for class 'bru_mapper_intrinsic_matern'
ibm_values(mapper, ...)

# S3 method for class 'bru_mapper_intrinsic_matern'
ibm_jacobian(mapper, input, ...)
```

## Arguments

- model:

  An `intrinsic_matern` object for which to construct or extract a
  mapper

- ...:

  Arguments passed on to other methods

- mapper:

  A `bru_mapper_intrinsic_matern` object

- input:

  The values for which to produce a mapping matrix
