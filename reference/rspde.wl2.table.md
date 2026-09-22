# Weighted-L2 coefficients for a given spectral interval

Builds a table of weighted-L2 rational coefficients for a specific lower
end of the spectral interval, so that the approximation is optimal on
the interval the discrete operator actually has rather than on
\\(0,1\]\\.

The mesh-free tables that
[`matern.operators()`](https://davidbolin.github.io/rSPDE/reference/matern.operators.md)
uses by default do not depend on the mesh or on \\\kappa\\ and are
stored in the package, which is why they cost nothing at set-up. They
are, however, fitted on all of \\(0,1\]\\, while a discretised operator
only ever sees \\x \in (x\_{\min}, 1\]\\ with \\x\_{\min}\\ set by the
mesh resolution and by \\\kappa\\. Fitting on that shorter interval
spends the same number of terms where they matter and is appreciably
more accurate; the price is that the table depends on the mesh and
therefore has to be computed, which takes of the order of a second.

With the cache turned on (see
[`rspde.cache()`](https://davidbolin.github.io/rSPDE/reference/rspde.cache.md))
there is usually nothing to do by hand: a model asked for with
`kappa_ref` computes its table once and finds it again in later sessions
without being told where it is. The two routes share one store, so a
table built here is what a later `matern.operators(kappa_ref = ...)`
picks up, and the other way round.

This function is for what the cache does not cover: holding
\\x\_{\min}\\ fixed across several meshes or models rather than letting
each derive its own, inspecting the coefficients and the errors of the
fit, or working without writing to a cache at all. Pass the result to
[`matern.operators()`](https://davidbolin.github.io/rSPDE/reference/matern.operators.md)
as `wl2_table`; supplying it overrides `x_min` and `kappa_ref`. The
result is an ordinary data frame with attributes, so
[`saveRDS()`](https://rdrr.io/r/base/readRDS.html) is also an option if
you would rather manage it yourself than enable the cache.

Since \\x\_{\min}\\ grows with \\\kappa\\, a *lower bound* for
\\\kappa\\ gives the shortest interval that is certainly long enough,
and hence a table that stays valid while \\\kappa\\ is estimated. That
is what `kappa_ref` is.

## Usage

``` r
rspde.wl2.table(
  m,
  d = NULL,
  nu = NULL,
  alpha = NULL,
  kappa_ref = NULL,
  x_min = NULL,
  C = NULL,
  G = NULL,
  mesh = NULL,
  loc_mesh = NULL,
  type = c("covariance", "operator"),
  eigenvalue = c("bound", "exact"),
  s = 0,
  k_term = FALSE,
  ...
)
```

## Arguments

- m:

  The order of the rational approximation.

- d:

  The dimension of the domain. Taken from `mesh` when possible.

- nu, alpha:

  The smoothness. Give one of them; `alpha = nu + d/2`. Only the integer
  part of `alpha` is used, since one table covers a whole unit interval
  of `alpha`.

- kappa_ref:

  A lower bound for `kappa`, from which `x_min` is computed.

- x_min:

  The lower end of the spectral interval, if it is known; then neither
  `kappa_ref` nor the mesh is needed.

- C, G, mesh, loc_mesh:

  The finite element matrices or the mesh, used with `kappa_ref` to find
  `x_min`. See
  [`rspde.xmin()`](https://davidbolin.github.io/rSPDE/reference/rspde.xmin.md).

- type:

  Either `"covariance"` or `"operator"`, matching the model the table is
  for.

- eigenvalue:

  Passed to
  [`rspde.xmin()`](https://davidbolin.github.io/rSPDE/reference/rspde.xmin.md).

- s, k_term:

  The weight exponent and the constant term of the fit; see
  [`rational.coefficients.wl2()`](https://davidbolin.github.io/rSPDE/reference/rational.coefficients.wl2.md).
  These change what is being approximated and are for studying the
  approximation rather than for building a model: a constant term has no
  block in the covariance-based construction, and a table carrying one
  is refused by
  [`matern.operators()`](https://davidbolin.github.io/rSPDE/reference/matern.operators.md).

- ...:

  Passed to the fit, for instance `n_starts` or `by`.

## Value

A table of coefficients, with the attributes `type`, `d`, `m`,
`m_alpha`, `x_min` and `kind`, for use as the `wl2_table` argument of
[`matern.operators()`](https://davidbolin.github.io/rSPDE/reference/matern.operators.md).

## See also

[`rspde.xmin()`](https://davidbolin.github.io/rSPDE/reference/rspde.xmin.md),
[`matern.operators()`](https://davidbolin.github.io/rSPDE/reference/matern.operators.md),
[`rational.coefficients.wl2()`](https://davidbolin.github.io/rSPDE/reference/rational.coefficients.wl2.md)

## Examples

``` r
# \donttest{
mesh <- fmesher::fm_mesh_1d(seq(0, 1, length.out = 101))
tab <- rspde.wl2.table(m = 2, d = 1, nu = 0.8, kappa_ref = 5, mesh = mesh)
attr(tab, "x_min")
#> [1] 0.0006246096
# }
```
