# Augment data with information from a `rspde_lme` object

Augment accepts a model object and a dataset and adds information about
each observation in the dataset. It includes predicted values in the
`.fitted` column, residuals in the `.resid` column, and standard errors
for the fitted values in a `.se.fit` column. It also contains the New
columns always begin with a . prefix to avoid overwriting columns in the
original dataset.

## Usage

``` r
# S3 method for class 'rspde_lme'
augment(
  x,
  newdata = NULL,
  loc = NULL,
  mesh = FALSE,
  which_repl = NULL,
  se_fit = FALSE,
  conf_int = FALSE,
  pred_int = FALSE,
  level = 0.95,
  n_samples = 100,
  ...
)
```

## Arguments

- x:

  A `rspde_lme` object.

- newdata:

  A `data.frame` or a `list` containing the covariates, the edge number
  and the distance on edge for the locations to obtain the prediction.
  If `NULL`, the fitted values will be given for the original locations
  where the model was fitted.

- loc:

  Prediction locations. Can either be a `data.frame`, a `matrix` or a
  character vector, that contains the names of the columns of the
  coordinates of the locations. For models using `metric_graph` objects,
  plase use `edge_number` and `distance_on_edge` instead.

- mesh:

  Obtain predictions for mesh nodes? The graph must have a mesh, and
  either `only_latent` is set to TRUE or the model does not have
  covariates.

- which_repl:

  Which replicates to obtain the prediction. If `NULL` predictions will
  be obtained for all replicates. Default is `NULL`.

- se_fit:

  Logical indicating whether or not a .se.fit column should be added to
  the augmented output. If TRUE, it only returns a non-NA value if type
  of prediction is 'link'.

- conf_int:

  Logical indicating whether or not confidence intervals for the fitted
  variable should be built.

- pred_int:

  Logical indicating whether or not prediction intervals for future
  observations should be built.

- level:

  Level of confidence and prediction intervals if they are constructed.

- n_samples:

  Number of samples when computing prediction intervals.

- ...:

  Additional arguments. Expert use only.

## Value

A
[`tidyr::tibble()`](https://tidyr.tidyverse.org/reference/reexports.html)
with columns:

- `.fitted` Fitted or predicted value.

- `.fittedlwrconf` Lower bound of the confidence interval, if conf_int =
  TRUE

- `.fitteduprconf` Upper bound of the confidence interval, if conf_int =
  TRUE

- `.fittedlwrpred` Lower bound of the prediction interval, if pred_int =
  TRUE

- `.fitteduprpred` Upper bound of the prediction interval, if pred_int =
  TRUE

- `.fixed` Prediction of the fixed effects.

- `.random` Prediction of the random effects.

- `.resid` The ordinary residuals, that is, the difference between
  observed and fitted values.

- `.se_fit` Standard errors of fitted values, if se_fit = TRUE.

## See also

[glance.rspde_lme](https://davidbolin.github.io/rSPDE/reference/glance.rspde_lme.md)
