# Observation matrix for space-time models

Observation matrix for space-time models

## Usage

``` r
rSPDE.Ast(
  mesh_space = NULL,
  space_loc = NULL,
  mesh_time = NULL,
  time_loc = NULL,
  graph = NULL,
  obs.s = NULL,
  obs.t = NULL,
  ind.field = NULL
)
```

## Arguments

- mesh_space:

  mesh object for models on 1d or 2d domains

- space_loc:

  mesh locations for models on 1d domains

- mesh_time:

  mesh object for time discretization

- time_loc:

  mesh locations for time discretization

- graph:

  MetricGraph object for models on metric graphs

- obs.s:

  spatial locations of observations

- obs.t:

  time points for observations

- ind.field:

  indices for field (if Dirichlet conditions are used)

## Value

Observation matrix linking observation locations to mesh nodes

## Examples

``` r
s <- seq(from = 0, to = 20, length.out = 11)
t <- seq(from = 0, to = 20, length.out = 5)
n.obs <- 10
obs.loc <- data.frame(x = max(s)*runif(n.obs), 
t = max(t)*runif(n.obs))
A <- rSPDE.Ast(space_loc = s,time_loc = t, 
               obs.s = obs.loc$x, obs.t = obs.loc$t)
```
