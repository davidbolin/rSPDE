# Transform Anisotropic SPDE Model Parameters to Original Scale

This function takes a vector of transformed parameters and applies the
appropriate transformations to return them in the original scale for use
in anisotropic SPDE models.

## Usage

``` r
transform_parameters_anisotropic(theta, nu_upper_bound = NULL)
```

## Arguments

- theta:

  A numeric vector of length 4 or 5, containing the transformed
  parameters in this order:

  lhx

  :   The logarithmic representation of hx.

  lhy

  :   The logarithmic representation of hy.

  logit_hxy

  :   The logit-transformed representation of hxy.

  lsigma

  :   The logarithmic representation of sigma.

  lnu (optional)

  :   The logarithmic representation of nu. If not provided, nu is not
      returned.

- nu_upper_bound:

  (optional) A numeric value representing the upper bound for the
  smoothness parameter nu. This is only used, and must be provided, if
  `lnu` is provided.

## Value

A named list with the parameters in the original scale:

- hx:

  The original scale for hx (exponential of lhx).

- hy:

  The original scale for hy (exponential of lhy).

- hxy:

  The original scale for hxy (inverse logit transformation of
  logit_hxy).

- sigma:

  The original scale for sigma (exponential of lsigma).

- nu (optional):

  The original scale for nu (using the forward_nu transformation). Only
  included if `lnu` is provided.

## Examples

``` r
# With lnu
theta <- c(log(0.1), log(0.2), log((0.3 + 1) / (1 - 0.3)), log(0.5), log(1))
nu_upper_bound <- 2
transform_parameters_anisotropic(theta, nu_upper_bound)
#> $hx
#> [1] 0.1
#> 
#> $hy
#> [1] 0.2
#> 
#> $hxy
#> [1] 0.3
#> 
#> $sigma
#> [1] 0.5
#> 
#> $nu
#> [1] 1
#> 

# Without lnu
theta <- c(log(0.1), log(0.2), log((0.3 + 1) / (1 - 0.3)), log(0.5))
transform_parameters_anisotropic(theta)
#> $hx
#> [1] 0.1
#> 
#> $hy
#> [1] 0.2
#> 
#> $hxy
#> [1] 0.3
#> 
#> $sigma
#> [1] 0.5
#> 
```
