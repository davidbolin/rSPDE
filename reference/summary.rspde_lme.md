# Summary Method for `rspde_lme` Objects.

Function providing a summary of results related to mixed effects
regression models with Whittle-Matern latent models.

## Usage

``` r
# S3 method for class 'rspde_lme'
summary(object, all_times = FALSE, ...)
```

## Arguments

- object:

  an object of class "rspde_lme" containing results from the fitted
  model.

- all_times:

  Show all computed times.

- ...:

  not used.

## Value

An object of class `summary_rspde_lme` containing several informations
of a *rspde_lme* object.
