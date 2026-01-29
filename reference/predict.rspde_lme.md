# Prediction of a mixed effects regression model on a metric graph.

Prediction of a mixed effects regression model on a metric graph.

## Usage

``` r
# S3 method for class 'rspde_lme'
predict(
  object,
  newdata = NULL,
  loc = NULL,
  time = NULL,
  mesh = FALSE,
  which_repl = NULL,
  compute_variances = FALSE,
  posterior_samples = FALSE,
  n_samples = 100,
  sample_latent = FALSE,
  return_as_list = FALSE,
  return_original_order = TRUE,
  ...,
  data = lifecycle::deprecated()
)
```

## Arguments

- object:

  The fitted object with the
  [`rspde_lme()`](https://davidbolin.github.io/rSPDE/reference/rspde_lme.md)
  function

- newdata:

  A `data.frame` or a `list` containing the covariates, the edge number
  and the distance on edge for the locations to obtain the prediction.

- loc:

  Prediction locations. Can either be a `data.frame`, a `matrix` or a
  character vector, that contains the names of the columns of the
  coordinates of the locations. For models using `metric_graph` objects,
  plase use `edge_number` and `distance_on_edge` instead.

- time:

  Prediction times for spatio-temporal models.

- mesh:

  Obtain predictions for mesh nodes? The graph must have a mesh, and
  either `only_latent` is set to TRUE or the model does not have
  covariates.

- which_repl:

  Which replicates to use? If `NULL` all replicates will be used.

- compute_variances:

  Set to also TRUE to compute the kriging variances.

- posterior_samples:

  If `TRUE`, posterior samples will be returned.

- n_samples:

  Number of samples to be returned. Will only be used if `sampling` is
  `TRUE`.

- sample_latent:

  Do posterior samples only for the random effects?

- return_as_list:

  Should the means of the predictions and the posterior samples be
  returned as a list, with each replicate being an element?

- return_original_order:

  Should the results be return in the original (input) order or in the
  order inside the graph?

- ...:

  Additional arguments. Expert use only.

- data:

  **\[deprecated\]** Use `newdata` instead.

## Value

A list with elements `mean`, which contains the means of the
predictions, `fe_mean`, which is the prediction for the fixed effects,
`re_mean`, which is the prediction for the random effects, `variance`
(if `compute_variance` is `TRUE`), which contains the variances of the
predictions, `samples` (if `posterior_samples` is `TRUE`), which
contains the posterior samples.
