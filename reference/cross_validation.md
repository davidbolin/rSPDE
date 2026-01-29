# Perform cross-validation on a list of fitted models.

Obtain several scores for a list of fitted models according to a folding
scheme.

## Usage

``` r
cross_validation(
  models,
  model_names = NULL,
  scores = c("mae", "mse", "crps", "scrps", "dss"),
  cv_type = c("k-fold", "loo", "lpo"),
  weight_thr = NULL,
  k = 5,
  percentage = 20,
  number_folds = 10,
  n_samples = 1000,
  return_scores_folds = FALSE,
  orientation_results = c("negative", "positive"),
  include_best = TRUE,
  train_test_indexes = NULL,
  return_train_test = FALSE,
  return_post_samples = FALSE,
  return_true_test_values = FALSE,
  parallelize_RP = FALSE,
  n_cores_RP = parallel::detectCores() - 1,
  true_CV = TRUE,
  save_settings = FALSE,
  model_options_bru = list(),
  print = TRUE,
  fit_verbose = FALSE
)
```

## Arguments

- models:

  A fitted model obtained from calling the
  [`bru()`](https://inlabru-org.github.io/inlabru/reference/bru.html)
  function or a list of fitted models. All models in the list must have
  the same number of likelihoods and must be fitted to identical
  datasets.

- model_names:

  A vector containing the names of the models to appear in the returned
  `data.frame`. If `NULL`, the names will be of the form `Model 1`,
  `Model 2`, and so on. By default, it will try to obtain the name from
  the models list.

- scores:

  A vector containing the scores to be computed. The options are "mse",
  "crps", "scrps", "dss", "wcrps" and "swcrps". By default, all scores
  are computed.

- cv_type:

  The type of the folding to be carried out. The options are `k-fold`
  for `k`-fold cross-validation, in which case the parameter `k` should
  be provided, `loo`, for leave-one-out and `lpo` for
  leave-percentage-out, in this case, the parameter `percentage` should
  be given, and also the `number_folds` with the number of folds to be
  done. The default is `k-fold`.

- weight_thr:

  When computing "wcrps" or "swcrps", the threshold to be used to
  compute the weights. Must be supplied if any of these scores are
  requested. No default value is provided.

- k:

  The number of folds to be used in `k`-fold cross-validation. Will only
  be used if `cv_type` is `k-fold`.

- percentage:

  The percentage (from 1 to 99) of the data to be used to train the
  model. Will only be used if `cv_type` is `lpo`.

- number_folds:

  Number of folds to be done if `cv_type` is `lpo`.

- n_samples:

  Number of samples to compute the posterior statistics to be used to
  compute the scores.

- return_scores_folds:

  If `TRUE`, the scores for each fold will also be returned.

- orientation_results:

  character vector. The options are "negative" and "positive". If
  "negative", the smaller the scores the better. If "positive", the
  larger the scores the better.

- include_best:

  Should a row indicating which model was the best for each score be
  included?

- train_test_indexes:

  A list where each element corresponds to a fold. Each fold contains:

  - `train`: A list of training index vectors, one for each likelihood.

  - `test`: A list of test index vectors, one for each likelihood, with
    the same length as `train`. This list is typically obtained by
    setting the argument `return_train_test` to `TRUE`. When supplying
    `train_test_indexes`, the `cv_type`, `k`, `percentage` and
    `number_folds` arguments are ignored.

- return_train_test:

  Logical. Should the training and test indexes be returned? If 'TRUE'
  the train and test indexes will the 'train_test' element of the
  returned list.

- return_post_samples:

  If `TRUE` the posterior samples will be included in the returned list.

- return_true_test_values:

  If `TRUE` the true test values will be included in the returned list.

- parallelize_RP:

  Logical. Should the computation of CRPS and SCRPS (and for some cases,
  DSS) be parallelized?

- n_cores_RP:

  Number of cores to be used if `parallelize_rp` is `TRUE`.

- true_CV:

  Should a `TRUE` cross-validation be performed? If `TRUE` the models
  will be fitted on the training dataset. If `FALSE`, the parameters
  will be kept fixed at the ones obtained in the result object.

- save_settings:

  Logical. If `TRUE`, the settings used in the cross-validation will
  also be returned.

- model_options_bru:

  A list of options to be passed to inlabru.

- print:

  Should partial results be printed throughout the computation?

- fit_verbose:

  Should INLA's run during cross-validation be verbose?

## Value

A data.frame with the fitted models and the corresponding scores.
