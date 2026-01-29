# Data frame for rspde_result objects to be used in ggplot2

Returns a ggplot-friendly data-frame with the marginal posterior
densities.

## Usage

``` r
# S3 method for class 'rspde_result'
gg_df(
  result,
  parameter = result$params,
  transform = TRUE,
  restrict_x_axis = NULL,
  restrict_quantiles = NULL,
  ...
)
```

## Arguments

- result:

  An rspde_result object.

- parameter:

  Vector. Which parameters to get the posterior density in the
  data.frame? The options are `std.dev`, `range`, `tau`, `kappa` and
  `nu`.

- transform:

  Should the posterior density be given in the original scale?

- restrict_x_axis:

  Variables to restrict the range of x axis based on quantiles.

- restrict_quantiles:

  Named list of quantiles to restrict x axis. It should contain the name
  of the parameter along with a vector with two elements specifying the
  lower and upper quantiles. The names should be match the ones in
  result\$params. For example, if we want to restrict nu to the 0.05 and
  0.95 quantiles we do `restrict_quantiles = c(0.05, 0.95)`.

- ...:

  currently not used.

## Value

A data frame containing the posterior densities.
