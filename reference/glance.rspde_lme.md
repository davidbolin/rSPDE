# Glance at an `rspde_lme` object

Glance accepts a `rspde_lme` object and returns a
[`tidyr::tibble()`](https://tidyr.tidyverse.org/reference/reexports.html)
with exactly one row of model summaries. The summaries are the square
root of the estimated variance of the measurement error, residual
degrees of freedom, AIC, BIC, log-likelihood, the type of latent model
used in the fit and the total number of observations.

## Usage

``` r
# S3 method for class 'rspde_lme'
glance(x, ...)
```

## Arguments

- x:

  An `rspde_lme` object.

- ...:

  Currently not used.

## Value

A
[`tidyr::tibble()`](https://tidyr.tidyverse.org/reference/reexports.html)
with exactly one row and columns:

- `nobs` Number of observations used.

- `sigma` the square root of the estimated residual variance

- `logLik` The log-likelihood of the model.

- `AIC` Akaike's Information Criterion for the model.

- `BIC` Bayesian Information Criterion for the model.

- `deviance` Deviance of the model.

- `df.residual` Residual degrees of freedom.

- `model.type` Type of latent model fitted.

## See also

[augment.rspde_lme](https://davidbolin.github.io/rSPDE/reference/augment.rspde_lme.md)
