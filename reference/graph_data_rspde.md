# Data extraction from metric graphs for 'rSPDE' models

Extracts data from metric graphs to be used by 'INLA' and 'inlabru'.

## Usage

``` r
graph_data_rspde(
  graph_rspde,
  name = "field",
  repl = NULL,
  repl_col = NULL,
  group = NULL,
  group_col = NULL,
  only_pred = FALSE,
  time = NULL,
  bru = FALSE,
  tibble = FALSE,
  drop_na = FALSE,
  drop_all_na = TRUE
)
```

## Arguments

- graph_rspde:

  An `inla_metric_graph_spde` or `inla_rspde_spacetime` object built
  with the
  [`rspde.metric_graph()`](https://davidbolin.github.io/rSPDE/reference/rspde.metric_graph.md)
  or
  [`rspde.spacetime()`](https://davidbolin.github.io/rSPDE/reference/rspde.spacetime.md)
  function.

- name:

  A character string with the base name of the effect.

- repl:

  Which replicates? If there is no replicates, one can set `repl` to
  `NULL`. If one wants all replicates, then one sets to `repl` to
  `.all`.

- repl_col:

  Which "column" of the data contains the replicate variable?

- group:

  Which groups? If there is no groups, one can set `group` to `NULL`. If
  one wants all groups, then one sets to `group` to `.all`.

- group_col:

  Which "column" of the data contains the group variable?

- only_pred:

  Should only return the `data.frame` to the prediction data?

- time:

  Column containing times for space time models. Not needed when using
  inlabru. Only for INLA implementation of space time model.

- bru:

  Should the data be processed for `inlabru`?

- tibble:

  Should the data be returned as a
  [`tidyr::tibble`](https://tibble.tidyverse.org/reference/tibble.html)?

- drop_na:

  Should the rows with at least one NA for one of the columns be
  removed? DEFAULT is `FALSE`. This option is turned to `FALSE` if
  `only_pred` is `TRUE`.

- drop_all_na:

  Should the rows with all variables being NA be removed? DEFAULT is
  `TRUE`. This option is turned to `FALSE` if `only_pred` is `TRUE`.

## Value

An 'INLA' and 'inlabru' friendly list with the data.
