# Posterior cross-validation for `rspde_lme` models

Performs cross-validation on objects fitted with
[`rspde_lme`](https://davidbolin.github.io/rSPDE/reference/rspde_lme.md)
and computes common predictive scores. The interface mirrors
[`MetricGraph::posterior_crossvalidation`](https://davidbolin.github.io/MetricGraph/reference/posterior_crossvalidation.html)
so that the two can be used interchangeably for graph-based and
non-graph rSPDE fits. Pure OLS fits (`rspde_lme(formula, data)` with no
`model`) are also supported: predictions reduce to evaluating the
covariates at the test locations, with predictive variance
\\\sigma\_\epsilon^2 \\ (1 + x^\top (X^\top X)^{-1} x)\\.

For `true_CV = FALSE` (default, "pseudo" cross-validation) the model
parameters are kept fixed at the estimates from the full fit and only
the held-out points are masked at prediction time. When
`use_precomputed = TRUE` the parameter-dependent structures (updated
rSPDE operator and the precision matrix Q) are computed once and reused
across folds, mirroring the `advanced_options$precompute_data` /
`na_test_idx` path of `predict.rspde_lme`.

For `true_CV = TRUE` the model is refit on each training fold (with the
held-out response values set to `NA`). The previous fit is forwarded via
`previous_fit` so that the optimisation starts from the full-data
estimates.

## Usage

``` r
posterior_crossvalidation(
  object,
  scores = c("logscore", "crps", "scrps", "mae", "rmse"),
  mode = "k-fold",
  k = 10,
  percentage = 20,
  number_folds = 10,
  train_test_indices = NULL,
  true_CV = FALSE,
  factor = 1,
  tibble = TRUE,
  parallel_folds = FALSE,
  parallel_fitting = FALSE,
  n_cores = parallel::detectCores() - 1,
  print = FALSE,
  seed = NULL,
  return_indices = FALSE,
  use_precomputed = TRUE,
  data = NULL,
  nelder_mead_init = FALSE
)
```

## Arguments

- object:

  A fitted object of class `rspde_lme`, or a (preferably named) list of
  such objects. When a list is supplied, the function is applied to each
  element and the scores are returned in a single `data.frame` /
  `tibble`.

- scores:

  Character vector of scores to compute. Possible values are
  `"logscore"`, `"crps"`, `"scrps"`, `"mae"` and `"rmse"`.

- mode:

  One of `"k-fold"`, `"loo"` or `"lpo"`.

- k:

  Number of folds for k-fold cross-validation. Default is 10.

- percentage:

  Percentage (1-99) of observations used for training in
  leave-percentage-out (`"lpo"`) cross-validation. Default is 20.

- number_folds:

  Number of folds for `"lpo"`. Default is 10.

- train_test_indices:

  Optional pre-specified list of folds. Each element must be a list with
  integer vectors `train` and `test` giving positions into the fitted
  observations (length `object$nobs`). When supplied, `mode`, `k`,
  `percentage` and `number_folds` are ignored.

- true_CV:

  If `TRUE` the model is refit on every training fold; if `FALSE`
  (default) the original fitted parameters are reused (pseudo
  cross-validation).

- factor:

  Multiplier applied to the mean of each score (default 1).

- tibble:

  If `TRUE` (default), the returned scores are coerced to a `tibble`
  (requires the tidyr package).

- parallel_folds:

  Process folds in parallel. Default `FALSE`.

- parallel_fitting:

  Run the per-fold model refit in parallel (only relevant when
  `true_CV = TRUE`). Default `FALSE`. `parallel_folds` and
  `parallel_fitting` cannot both be `TRUE`.

- n_cores:

  Number of cores for parallel computation. Default
  `parallel::detectCores() - 1`.

- print:

  Print progress messages.

- seed:

  Random seed for fold creation reproducibility.

- return_indices:

  If `TRUE` the train/test indices used in each fold are also returned.

- use_precomputed:

  If `TRUE` (default) the parameter-dependent structures used by
  `predict.rspde_lme` are precomputed once and reused across folds. Only
  relevant when `true_CV = FALSE`.

- data:

  Optional `data.frame` containing the original data used to fit the
  model. Required for non-graph fits with non-trivial covariates (e.g.
  factors or transformations), where the design matrix stored in
  `object$model_matrix` is insufficient to rebuild the covariates at the
  test locations. Ignored for graph-based fits, where the data is taken
  from the metric graph.

- nelder_mead_init:

  Logical, forwarded to
  [`rspde_lme`](https://davidbolin.github.io/rSPDE/reference/rspde_lme.md)
  when `true_CV = TRUE`. Default is `FALSE`: per-fold refits warm-start
  from the full-fit parameters via `previous_fit`, so the Nelder-Mead
  pre-pass adds little except cost. Set to `TRUE` to re-enable it for
  robustness if the full fit is itself suspected of being a poor local
  optimum.

## Value

A list with elements

- mu:

  vector of posterior means, one per observation

- var:

  vector of posterior variances (including measurement error), one per
  observation

- scores:

  a `data.frame` (or `tibble`) with the requested aggregated scores

- indices:

  list of train/test indices used (if `return_indices = TRUE`)
