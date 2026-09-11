# Posterior cross-validation for fitted models

Generic function for posterior cross-validation. rSPDE provides methods
for
[`rspde_lme`](https://davidbolin.github.io/rSPDE/reference/rspde_lme.md)
fits (see
[`posterior_crossvalidation.rspde_lme`](https://davidbolin.github.io/rSPDE/reference/posterior_crossvalidation.rspde_lme.md))
and for lists of fitted models, and the MetricGraph package provides a
method for `graph_lme` fits. Because both packages use this generic,
loading them in either order does not mask one implementation with the
other.

## Usage

``` r
posterior_crossvalidation(object, ...)

# S3 method for class 'list'
posterior_crossvalidation(object, ..., tibble = TRUE, return_indices = FALSE)

# Default S3 method
posterior_crossvalidation(object, ...)
```

## Arguments

- object:

  A fitted model, or a (preferably named) list of fitted models. The
  elements of a list can be of any class with a
  `posterior_crossvalidation` method, so for instance `rspde_lme` and
  `graph_lme` fits can be compared in one call.

- ...:

  Arguments passed on to the methods. For a list of models, only the
  arguments that are supplied are passed on, so each method uses its own
  defaults for the rest.

- tibble:

  If `TRUE` (default), the scores for a list of models are returned as a
  `tibble` with a `Model` column.

- return_indices:

  If `TRUE`, the train/test indices used for the first model are also
  returned.

## Value

A list with elements `mu`, `var` and `scores`, and `indices` if
requested. For a list of models, `mu` and `var` are lists with one
element per model, and `scores` has one row per model.
