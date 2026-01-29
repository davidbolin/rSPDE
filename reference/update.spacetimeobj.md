# Update parameters of spacetimeobj objects

Function to change the parameters of a spacetimeobj object

## Usage

``` r
# S3 method for class 'spacetimeobj'
update(object, kappa = NULL, sigma = NULL, gamma = NULL, rho = NULL, ...)
```

## Arguments

- object:

  Space-time object created by
  [`spacetime.operators()`](https://davidbolin.github.io/rSPDE/reference/spacetime.operators.md)

- kappa:

  kappa value to be updated.

- sigma:

  sigma value to be updated.

- gamma:

  gamma value to be updated.

- rho:

  rho value to be updated.

- ...:

  currently not used.

## Value

An object of type spacetimeobj with updated parameters.

## Examples

``` r
s <- seq(from = 0, to = 20, length.out = 101)
t <- seq(from = 0, to = 20, length.out = 31)

op_cov <- spacetime.operators(space_loc = s, time_loc = t,
                             kappa = 5, sigma = 10, alpha = 1,
                             beta = 2, rho = 1, gamma = 0.05)
op_cov <- update(op_cov, kappa = 4, 
                             sigma = 2, gamma = 0.1)                              
```
